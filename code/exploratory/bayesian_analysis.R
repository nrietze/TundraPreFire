library(terra)
library(tidyterra)
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)
library(broom)
library(tictoc)
library(future)
library(gt)
set.seed(10)

# 0. Set up functions ----
# ========================.
load_data <- function(fire_attrs,burn_severity_index,frac_int,
                      return_df_only = FALSE){
  # Fire perimeter attributes
  UTM_TILE_ID <- fire_attrs$opt_UTM_tile
  year <- fire_attrs$tst_year
  FIRE_ID <- fire_attrs$fireid
  
  # Load filtered and formatted dataframe
  fn_filtered_df <- paste0(TABLE_DIR,sprintf(
    "model_dataframes/%spct/%s_model_dataframe.csv",
    frac_int,FIRE_ID
  ))
  
  if (!file.exists(fn_filtered_df)){
    return(NULL)
  }
  
  df_filtered <- read_csv2(fn_filtered_df, col_names = TRUE,show_col_types = FALSE) %>% 
    as_tibble() %>% 
    mutate(date = as.Date(date),
           burn_date = as.Date(burn_date),
           doy = yday(burn_date),
           fireid = FIRE_ID)
  
  if (return_df_only) {return(df_filtered)}
  
  # Load optimal burn severity raster
  fname_optimal_severity_raster <- optimality_lut %>% 
    filter(fireid == FIRE_ID,
           severity_index == burn_severity_index) %>% 
    pull(fname_severity_raster)
  
  severity_raster <- rast(fname_optimal_severity_raster)
  
  # Subset to single perimeter
  selected_fire_perimeter <- fire_perimeters %>% 
    filter(fireid  == FIRE_ID) %>%
    project(crs(severity_raster))
  
  # Buffer the selected fire perimeter
  fire_perimeter_buffered <- buffer(selected_fire_perimeter, 1200)
  
  cropped_severity_raster <- severity_raster %>% 
    mask(fire_perimeter_buffered, updatevalue = NA) %>% 
    crop(fire_perimeter_buffered)
  
  # Load sample points (with associated burn dates)
  fname_sample_points <-  paste0(
    sprintf("%s/feature_layers/%spct/%s_sample_points_%spct_burn_date.gpkg",
            DATA_DIR,frac_int,FIRE_ID,frac_int)
  )
  sample_points <- vect(fname_sample_points) %>% 
    project(crs(severity_raster)) %>% 
    mutate(ObservationID = 1:nrow(.))
  
  # Extract data at random points 
  severity_raster_sample <- terra::extract(severity_raster, sample_points,
                                           ID = FALSE, xy = TRUE)
  
  return(list(df_filtered,
              severity_raster_sample,
              sample_points,
              cropped_severity_raster,
              selected_fire_perimeter))
}

# 1. Configure and load stuff ----
# ================================.
# Config, loading and preparing data
OS <- Sys.info()[['sysname']]
if(OS == "Linux"){
  
  HLS_DIR <- "~/scratch/raster/hls/"
  TABLE_DIR <- "~/data/tables/"
  DATA_DIR <- "~/data/"
} else {
  
  HLS_DIR <- "data/raster/hls/"
  TABLE_DIR <- "data/tables/"
  DATA_DIR <- "data/"
}

index_name <- "NDMI"
burn_severity_index <- "dNBR"
pct_cutoff <- 0.5
SAVE_FIGURES <- TRUE

frac_to_sample <- 0.01
frac_int <- frac_to_sample *100
SCALE_FACTOR <- ifelse(burn_severity_index == "dNBR", 1000, 1)

# Load lookup tables
processing_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017)

optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

# Load fire perimeters
fire_perimeters <- vect(
  sprintf("%s/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg",DATA_DIR)
)

## Build global model data table ----
topN_fires <- fire_perimeters %>%
  filter(tst_year >= 2017) %>% 
  arrange(desc(farea)) %>% 
  slice_head(n = 25) 

TEST_ID <- topN_fires$fireid

if (length(TEST_ID)>0){
  subset_lut <- filter(processing_lut, fireid %in% TEST_ID)
}

# Load model dataframe
fname_model_data <-  paste0(
  TABLE_DIR,sprintf("model_dataframes/%spct/final_model_dataframe.csv",frac_int))

if (file.exists(fname_model_data)){
  
  final_df <- read_csv2(fname_model_data)
  
} else {
  data_list <- subset_lut %>%
    split(seq(nrow(.))) %>%
    map(~ load_data(.x, burn_severity_index, frac_int,
                    return_df_only = TRUE)) %>%
    compact()
  
  bind_fill <- function(dfs) {
    all_cols <- unique(unlist(lapply(dfs, colnames)))  # all unique column names
    dfs_filled <- lapply(dfs, function(df) {
      missing <- setdiff(all_cols, names(df))       # columns that are missing
      df[missing] <- NA                              # add missing columns filled with NA
      df[all_cols]                                    # reorder to consistent column order
    })
    do.call(rbind, dfs_filled)
  }
  
  # Combine dataframes
  final_df <- bind_fill(data_list)
  final_df$burn_doy <- yday(final_df$burn_date)
  
  write_csv2(final_df,fname_model_data )
}

# subset of upper 75th percentile
df_subset <- final_df %>% 
  group_by(fireid) %>% 
  # filter(dnbr >= 0.245) %>%  # (threshold from Kolden et al. 2012)
  filter(quantile(dnbr, 0.75,na.rm=T)<dnbr) %>% 
  mutate(fireid = as.factor(fireid),
         # scale predictors
         burn_doy = burn_doy / 366,
         elevation = scale(elevation)[,1],
         across(contains("lst"), ~ as.numeric(scale(.)))) %>% 
  ungroup()

# 4. Run brms model ----
# ========================.

y_var <- "dnbr"

TEST_RUN_SINGELDAY <- FALSE

# Test run single-day model
if (TEST_RUN_SINGELDAY){
  day <- 10
  ndvi_var <- paste0("NDVI.d_prefire_", day)
  ndmi_var <- paste0("NDMI.d_prefire_", day)
  lst_var <- paste0("cumsum_lst_d_prefire_", day)
  
  predictors <- c(ndvi_var, ndmi_var,lst_var, "elevation", "slope", "northness", "eastness", "doy")
  fixed_effects <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))
  
  df_d10_subset <- final_df %>% 
    group_by(fireid) %>% 
    sample_frac(.01)
  
  tic()
  model <- brm(
    bf(formula),
    family = "gaussian",
    data = df_d10_subset,
    prior=set_prior("normal(0,1)", class="b"), # we expect small coefficients for slopes (b)
    chains = 4,
    iter = 10000,
    thin = 10,
    cores = 4,
    seed = 1234)
  toc()
  
  
} else{
  
  predictors <- c("NDVI", "NDMI","LST", "elevation", "slope", "northness", "eastness", "doy")
  fixed_effects <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))
  
  df_test_subset <- df_subset
  # df_test_subset <- df_subset %>% 
  #   group_by(fireid) %>% 
  #   sample_frac(.01)
  
  n_days <- 30
  
  df_list_by_day <- list()
  
  for (day in 1:n_days) {
    ndvi_col <- paste0("NDVI.d_prefire_", day)
    ndmi_col <- paste0("NDMI.d_prefire_", day)
    lst_col  <- paste0("cumsum_lst_d_prefire_", day)
    
    df_day <- df_test_subset[, c(ndvi_col, ndmi_col, lst_col,
                                 "elevation", "slope", "northness", "eastness",
                                 "doy","fireid","dnbr")]
    names(df_day) <- c(predictors,"fireid","dnbr")
    
    # Add to list
    df_list_by_day[[day]] <- df_day
  }
  
  plan(multisession, workers = 8)
  
  tic()
  all_models <- brm_multiple(
    bf(formula),
    family = "gaussian",
    data = df_list_by_day,
    prior=set_prior("normal(0,1)", class="b"), # we expect small coefficients for slopes (b)
    chains = 4,
    iter = 10000,
    thin = 10,
    seed = 1234,
    combine = FALSE
  )
  toc()
  saveRDS(all_models, file="data/tables/all_models_fuldata.RData")
  
  # tic()
  # model <- brm(
  #   bf(formula),
  #   family = "gaussian",
  #   data = df_d10_subset,
  #   prior=set_prior("normal(0,1)", class="b"),
  #   chains = 4,
  #   iter = 10000,
  #   thin = 10,
  #   cores = 4, 
  #   seed = 1234)
  # toc()
}



# 3. Model assessment ----

ggplot(all_models[[1]]) +
  geom_histogram(aes(x=dnbr), color = "tomato",fill = "salmon") +
  labs(x = "Burn severity (dNBR)",
       title = "Distribution of dNBR for benchmark subset") +
  theme_cowplot()

ggsave("figures/hist_dNBR_benchmark.png",bg = "white",
       width = 8, height = 8)


## a. Report posterior checks ----

model_summaries <- lapply(all_models, function(model) {
  list(
    ess = bayestestR::effective_sample(model),
    r2  = performance::r2_bayes(model)
  )
})

# View effective sample sizes from model 1
model_summaries[[1]]$ess %>% gt()

# View R² for model 1
r2s <- unlist(lapply(model_summaries, function(x) {x$r2$R2_Bayes_marginal}))

plot(r2s)

model <- all_models[[1]]

# Report effective sample size
bayestestR::effective_sample(model) %>% gt()

# Get pseudo R2
performance::r2_bayes(model)

# model sumamry
summary(model)

# plots
plot(model)

pp_check(model)

pp_check(model, type='error_scatter_avg')
pp_check(model, type='error_hist')

mod_labs <-c(elevation = 'Elevation',
             slope = 'Slope',
             northness = 'Northness',
             eastness = 'Eastness',
             doy = "Burn timing (Day of Year)",
             cumsum_lst_d_prefire_10 = "Land surface temperature \n(LST Landsat)",
             NDVI.d_prefire_10 = "Greenness (NDVI HLS)", 
             NDMI.d_prefire_10 = "Vegetation moisture content \n(NDMI HLS)")

## b. Plot fixed effects ----
(p <- model %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   mutate(`.variable` = gsub("b_", "", `.variable`) ) %>%
   mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>%
   # build plot 
   ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.95),size = 4,fill = "#E8CEB6") +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   scale_y_discrete(labels = mod_labs) +
   labs(x = "Coefficient estimate \n(scaled)",
        y = "") +
   theme_cowplot(18)
)


# Get random effects
ranef(model)$fireid %>% 
  as_tibble(rownames = "fireid")

## c. Plot random effects ----
(p_re <- model %>%
   gather_draws(`r_fireid\\[.*`, regex = TRUE) %>%          # filter estimated effects
   mutate(`.variable` = factor(`.variable`)) %>% 
   ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.05,0.95),size = 0.5,fill = '#b5ccb9') +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   scale_y_discrete(labels = mod_labs) +
   theme_minimal_hgrid() + 
   theme(panel.grid = element_blank(),
         legend.position = "none") +
   labs(x = "Posterior estimates (standardized)", y = ""))

## d. Plot all mode effects over time ----

# Function to extract posteriors andfor each day and plot per predictor
plot_predictor_effects <- function(models, days, predictors, palette = NULL) {
  
  predictor_regex <- paste0("b_", predictors, collapse = "|")
  
  # Extract draws for all predictors
  all_effects <- map2_dfr(models, days, ~{
    gather_draws(.x, !!sym(predictor_regex), regex = TRUE) %>%
      mutate(day_before_fire = .y)
  })
  
  # Clean predictor names
  all_effects <- all_effects %>%
    mutate(predictor = gsub("^b_", "", .variable))
  
  # Add significance flag (does 95% CI exclude 0?)
  significance_flags <- all_effects %>%
    group_by(day_before_fire, predictor) %>%
    summarise(
      lower = quantile(.value, 0.025),
      upper = quantile(.value, 0.975),
      sig = !(lower < 0 & upper > 0),
      .groups = "drop"
    )
  
  # Join with full draws
  all_effects <- left_join(all_effects, significance_flags, by = c("day_before_fire", "predictor"))
  
  # Optional color palette
  if (is.null(palette)) {
    palette <- scales::hue_pal()(length(predictors))
    names(palette) <- predictors
  }
  
  # Plot with colored predictors
  ggplot(all_effects, aes(x = day_before_fire, y = .value,
                          fill = predictor, alpha = sig)) +
    stat_halfeye(.width = 0.95, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3)) +
    scale_fill_manual(values = palette) +
    labs(
      x = "Days before fire",
      y = "Posterior estimate",
      fill = "Predictor",
      alpha = "Significant"
    ) +
    scale_x_reverse() +
    facet_wrap(~ predictor, scales = "free_y") +
    theme_cowplot(14)
}

plot_predictor_effects(models = all_models, days = 1:30, predictors = predictors)

ggsave("figures/dnbr_effect_sizes_brms_benchmark_10pct.png",bg = "white",
       height = 10, width = 14)
