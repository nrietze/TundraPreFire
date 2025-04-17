library(terra)
library(tidyterra)
library(spatstat.utils)
library(tidyverse)
library(cowplot)
library(patchwork)
library(lubridate)
library(lme4)
library(brms)
library(tidybayes)
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
    sprintf("%s/feature_layers/%s_sample_points_%spct_burn_date.gpkg",
            DATA_DIR,FIRE_ID,frac_int)
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

PlotSeverityMapHist <- function(cropped_severity_raster){
  burn_severity_index <- names(cropped_severity_raster)
  SCALE_FACTOR <- ifelse(burn_severity_index == "dNBR", 1000, 1)
  
  binwidth <- .02 * SCALE_FACTOR
  
  # Plot dNBR map
  p1 <- ggplot() +
    geom_spatraster(data = cropped_severity_raster*SCALE_FACTOR) +
    scale_fill_viridis_c(option = "inferno",
                         na.value = "white",
                         name = burn_severity_index) +
    theme_cowplot()
  
  # Plot histogram for this map
  xvar <- as.symbol(burn_severity_index)
  p2 <- ggplot(data = cropped_severity_raster, aes(x = !!xvar)) +
    geom_histogram(binwidth = binwidth, color = "black", fill = "indianred1") +
    xlim(c(-.5*SCALE_FACTOR,1*SCALE_FACTOR)) +
    labs(
      title = sprintf("Histogram of %s",burn_severity_index),
      x = sprintf("%s Values",burn_severity_index),
      y = "Frequency"
    ) +
    theme_cowplot()
  
  pg <- p1 + p2
  return(pg)
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
top20_fires <- fire_perimeters %>%
  arrange(desc(farea)) %>% 
  slice_head(n = 25) 

TEST_ID <- top20_fires$fireid
# TEST_ID <- c(14211,14664,10792,17548)

if (length(TEST_ID)>0){
  subset_lut <- filter(processing_lut, fireid %in% TEST_ID)
}

data_list <- subset_lut %>%
  split(seq(nrow(.))) %>%
  map(~ load_data(.x, burn_severity_index, frac_int,
                  return_df_only = TRUE)) %>% 
  compact() 

final_df <- do.call(rbind, data_list)

# Filter out data for burned points (threshold from Kolden et al. 2012)
final_df_subset_burnclass <- final_df %>% 
  mutate(dnbr = dnbr) %>% 
  filter(dnbr >= 0.245) %>% 
  mutate(fireid = as.factor(fireid))

# Load individual fire data
FIRE_ID <- 14664
if (!is.na(FIRE_ID)){
  subset_lut <- filter(processing_lut, fireid == FIRE_ID)
  
  data_list <- load_data(subset_lut, burn_severity_index = "dNBR",
                         frac_int)
  
  final_df_14664 <- data_list[[1]]
  severity_raster_sample <- data_list[[2]]
  sample_points <- data_list[[3]]
  cropped_severity_raster <- data_list[[4]]
  selected_fire_perimeter <- data_list[[5]]
  
  final_df_14664_burnclass <- final_df_14664 %>% 
    filter(dnbr >= 0.1)
  
  rm(data_list)
}

# 2. Descriptive statistics ----
# ==============================.

### Map & Histogram of burn severity data ----
pg <- PlotSeverityMapHist(cropped_severity_raster)
if (SAVE_FIGURES){
  ggsave2(pg,filename = sprintf("figures/%s_%s_map_hist.png",burn_severity_index,FIRE_ID),
          width = 8,height = 6,bg = "white")
}

### Plot distribution of total raster vs. sample points ----
xvar <- as.symbol(burn_severity_index)

p2 <- ggplot(data = cropped_severity_raster, aes(x = !!xvar)) +
  geom_histogram(binwidth = binwidth, color = "black", fill = "indianred1") +
  xlim(c(-.5*SCALE_FACTOR,1*SCALE_FACTOR)) +
  labs(
    title = sprintf("Histogram of %s",burn_severity_index),
    x = sprintf("%s Values",burn_severity_index),
    y = "Frequency"
  ) +
  theme_cowplot()

p3 <- ggplot(data = df_all_subset, aes(x = burn_severity * SCALE_FACTOR)) +
  geom_histogram(binwidth = binwidth, color = "black", fill = "orange") +
  xlim(c(-500,1000)) +
  labs(
    title = sprintf("Histogram of %s (%s%% sample)",burn_severity_index,frac*100),
    x = sprintf("%s Values",burn_severity_index),
    y = "Frequency"
  ) +
  theme_cowplot()

pg <- p2 / p3

if (SAVE_FIGURES){
  ggsave2(pg, 
          filename = sprintf("figures/%s_hist_vs_%spct_sampled_%s.png",
                             burn_severity_index,frac*100,FIRE_ID),
          width = 8, height = 6,bg = "white")
}

### Scatterplot NDMI vs. LST ----
ggplot(df_all_subset) +
  geom_point(aes(x = DailyMeanNDMI,y = LST,color = burn_severity),
             alpha = .5) +
  labs(x = "Canopy moisture (Daily averaged NDMI)",
       y = "Land surface temperature (Daily averaged LST)",
       color = sprintf("Burn severity (%s)",burn_severity_index)) +
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(filename = sprintf("figures/Scatter_LST_NDMI_%s_%s.png",
                             burn_severity_index,FIRE_ID),
          width = 8,height = 6,bg = "white")
}

### dNBR sensitivity over all dates (for largest scar) ----
UTM_TILE_ID <- "54WXE"
year <- 2020

severity_raster_list <- list.files(path = paste0(HLS_DIR,"severity_rasters"),
                                   pattern = sprintf("^%s_%s_%s.*\\.tif$",
                                                     burn_severity_index,UTM_TILE_ID,year),
                                   full.names = TRUE)

all_dnbr_rast_list <- lapply(severity_raster_list, rast)

# Extract dates (or full names) from filenames
raster_names <- basename(severity_raster_list) |> 
  tools::file_path_sans_ext()  # remove ".tif" if present

# Set the names of each raster before stacking
for (i in seq_along(all_dnbr_rast_list)) {
  names(all_dnbr_rast_list[[i]]) <- raster_names[i]
}

# Combine into one multi-layer raster
all_dnbr_rast <- rast(all_dnbr_rast_list)

# Convert raster to data frame
dnbr_df <- as.data.frame(all_dnbr_rast, xy = FALSE, na.rm = TRUE)

# Reshape to long format
dnbr_long <- dnbr_df %>%
  pivot_longer(cols = everything(), names_to = "date", values_to = "dNBR")

# Optional: Extract just the date part (if name is like "dNBR_54WXE_2020-09-04")
dnbr_long <- dnbr_long %>%
  mutate(date = str_extract(date, "\\d{4}-\\d{2}-\\d{2}"))

ggplot(dnbr_long, aes(x = dNBR)) +
  geom_histogram(bins = 50, fill = "tomato", color = "white") +
  facet_wrap(~ date, scales = "free") +
  labs(title = "dNBR composite histograms by date", 
       x = "dNBR", y = "Frequency") +
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(filename = sprintf("figures/%s_%s_histograms_all_dates.png",
                             burn_severity_index,FIRE_ID),
          width = 10,height = 10,bg = "white")
}



## Histogram of dNBR in sampled data ----
burn_severity_index <- "dnbr"

pdnbr <- ggplot(data = final_df) +
  geom_histogram(aes(x = dnbr, 
                     fill = ifelse(dnbr >= 0.245,
                                   "burned", "unburned")),
                 binwidth = .005) +
  facet_wrap(~fireid) +
  labs(x = "Burn Severity (dNBR)",
       fill = "",
       title = "Histograms per fire scar") +
  theme_cowplot()

pdgemi <- ggplot(data = final_df) +
  geom_histogram(aes(x = dgemi, 
                     fill = ifelse(dnbr >= 0.245,
                                   "burned", "unburned")),
                 binwidth = .005) +
  facet_wrap(~fireid) +
  labs(x = "Burn Severity (dgemi)",
       fill = "") +
  theme_cowplot()

fig <- pdnbr / pdgemi +
  plot_layout(axis_titles = "collect")

if (SAVE_FIGURES){
  ggsave2(fig,filename = "figures/Histograms_dnbrvs_dgemi_all_scars.png",
          width = 8,height = 6,bg = "white")
}


### Check spatial autocorrelation ----
# library(spatstat)
# library(spdep)
# 
# df_autocorr <- sample_points_subset %>% 
#   terra::as.data.frame(geom = "XY") %>% 
#   inner_join(df_all_subset %>% 
#                select(ObservationID, burn_severity) %>% 
#                distinct(),
#             by = "ObservationID")
# 
# knn <- knn2nb(knearneigh(df_autocorr[c("x","y")], k = 5))
# 
# # Convert neighbors to weights
# weights <- nb2listw(knn, style = "W")
# 
# # Compute Moran's I
# moran_i <- moran.test(df_autocorr$burn_severity, listw = weights)
# 
# (ggplot() +
#     geom_spatraster(data = mask(cropped_severity_raster,
#                                 selected_fire_perimeter)) +
#     scale_fill_viridis_c(option = "inferno",
#                          na.value = "white",
#                          name = burn_severity_index) +
#     geom_spatvector(data = sample_points_subset) +
#     labs(title = sprintf("Randomly samlped points\nMoran's I: %.2f (p-value: %.4f)", 
#                          moran_i$estimate[1],moran_i$p.value)) +                   
#     theme_cowplot())
# 
# ggsave2(sprintf("figures/%s_stratified_sample_points_moranI.png",FIRE_ID),
#         bg = "white",width = 10, height = 10)


# 3. Plot time series data ----
# ==============================.
# Plot time-integrated NDMI 7 days before fire vs. dNBR
xvar <- "NDMI.d_prefire_7"

df_all_subset <- final_df %>% 
  mutate(severity_group = if_else(burn_severity >= 0.245,
                                  "burned (≥ 0.245)",
                                  "unburned (< 0.245)"))

fig <- ggplot(df_all_subset, aes(x = !!sym(xvar),
                                 y = burn_severity,
                                 color = severity_group)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(x = sprintf("Cumulative %s (sum over 1 to 7 days pre-fire)", index_name),
       y = burn_severity_index,
       color = "Burn severity class \n (threshold from Kolden et al.(2012))") +
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(fig,
          filename = sprintf("figures/%s_%s_%s.png",
                             burn_severity_index,xvar,FIRE_ID),
          width = 8,height = 6,bg = "white")
}

# Plot time-integrated LST 7 days before fire vs. dNBR
ggplot(df_all, aes(x = cumsum_lst_d_prefire_7 ,
                   y = burn_severity,
                   color = descals_burn_class)) +
  geom_point(alpha = .2) +
  geom_smooth(method='lm') +
  labs(x = "Cumulative LST (sum over 1 to 7 days pre-fire)",
       y = burn_severity_index,
       color = "Burn class (Descals et al., 2022)") +
  theme_cowplot()


# Plot single-day predicted LST 7 days before fire vs. dNBR
ggplot(df_all, aes(x = lst_pred_d_prefire_7 ,
                   y = burn_severity, color = descals_burn_class)) +
  geom_point(alpha = .2) +
  geom_smooth(method='lm') +
  labs(x = "LST (predicted 7 days pre-fire)",
       y = burn_severity_index,
       color = "Burn class (Descals et al., 2022)") +
  theme_cowplot()


# 4. Run linear model ----
# ========================.
burn_severity_index <- "dnbr"

## a. mixed model ----
ndmi_prefire_cols <- grep("^NDMI.d_prefire_", names(final_df), value = TRUE)
ndvi_prefire_cols <- grep("^NDVI.d_prefire_", names(final_df), value = TRUE)

preds <- c(ndmi_prefire_cols, ndvi_prefire_cols)

# Function to run lm and extract statistics
run_lm_day <- function(day,y_var, data) {
  ndvi_var <- paste0("NDVI.d_prefire_", day)
  ndmi_var <- paste0("NDMI.d_prefire_", day)
  
  predictors <- c(ndvi_var, ndmi_var, "elevation", "slope", "northness", "eastness", "doy")
  fixed_effects <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))
  
  model <- lmer(formula, data = data, REML = FALSE)
  model_summary <- summary(model)
  
  # Extract fixed effect estimates
  coefs <- fixef(model)  # Named vector
  # pvals <- coef(summary(model))[, "Pr(>|t|)"]
  
  # R-squared using performance package (optional but helpful)
  r2_vals <- performance::r2(model)
  
  tibble(
    days_before_fire = day,
    r_squared = r2_vals$R2_marginal,
    adj_r_squared = NA  # Not directly available in lmer
  ) %>%
    bind_cols(
      as_tibble_row(
        setNames(as.list(coefs),
                 c("Intercept","NDVI", "NDMI",
                   "elevation", "slope", "northness", "eastness",
                   "doy")))
    )
}

# only model dNBR above burned threshold
model_results_by_day <- map_dfr(c(1:40), ~ run_lm_day(y_var = burn_severity_index,
                                                      day = .x,
                                                      data = final_df_subset_burnclass))
results <- model_results_by_day %>% arrange(days_before_fire)

fig <- ggplot(results) + 
  geom_point(aes(x = days_before_fire, y = r_squared)) +
  labs(x = "Days before fire", y = "R-squared",
       title = sprintf("Rsquared for %s",burn_severity_index)) + 
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(fig,filename = sprintf("figures/%s_%s_subset_R2_fullmod_burnedonly.png",burn_severity_index,FIRE_ID),
          width = 8,height = 6,bg = "white")
}


# Plot effect sizes over time
results %>% 
  pivot_longer(
    cols = -c(days_before_fire, r_squared, adj_r_squared),
    names_to = "coefficient", values_to = "values") %>% 
  ggplot() + 
  geom_point(aes(x = days_before_fire, y = values)) + 
  labs(x = "Days before fire", y = "Effect size") +
  facet_wrap(~coefficient,scales="free",ncol = 2) + 
  theme_cowplot()

## b. Run linear model on individual site  ----
run_lm <- function(day,y_var, data) {
  ndvi_var <- paste0("NDVI.d_prefire_", day)
  ndmi_var <- paste0("NDMI.d_prefire_", day)
  
  predictors <- c(ndvi_var, ndmi_var, "elevation", "slope", "northness", "eastness", "doy")
  fixed_effects <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(y_var, "~", fixed_effects))
  
  model <- lm(formula, data = data)
  model_summary <- summary(model)
  
  coefs <- coef(model)
  
  tibble(
    days_before_fire = day,
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared
  ) %>%
    bind_cols(
      as_tibble_row(
        setNames(as.list(coefs),
                 c("Intercept","NDVI", "NDMI",
                   "elevation", "slope", "northness", "eastness",
                   "doy")))
    )
}

# Run lm for each d_prefire column and combine results
burn_severity_index <- "dnbr"

results_all_burn_classes <- map_dfr(c(1:40),
                                    ~run_lm(burn_severity_index, day = .x,data = final_df_14664))
results_all_burn_classes <- results_all_burn_classes %>% arrange(days_before_fire)

# Plot R2
p1 <- ggplot(results_all_burn_classes) + 
  geom_point(aes(x = days_before_fire, y = r_squared)) + 
  ylim(c(0.15,0.47)) +
  labs(x = "Days before fire", y = "R-squared", 
       title = sprintf("%s ~ NDVI.d_prefire_x + NDMI.d_prefire_x + elevation + 
         slope + northness + eastness + doy",burn_severity_index)) +
  theme_cowplot()

results_burned <- map_dfr(c(1:40),
                          ~run_lm(burn_severity_index, day = .x,data = final_df_14664_burnclass))
results_burned <- results_burned %>% arrange(days_before_fire)

p2 <- ggplot(results_burned) + 
  geom_point(aes(x = days_before_fire, y = r_squared)) + 
  ylim(c(0.15,0.47)) +
  labs(x = "Days before fire", y = "R-squared") +
  theme_cowplot()

fig <- p1 + p2 +
  plot_layout(axis_titles = "collect")

if (SAVE_FIGURES){
  ggsave2(fig,filename = sprintf("figures/%s_%s_R2_comparison.png",burn_severity_index,FIRE_ID),
          width = 8,height = 6,bg = "white")
}


# Plot effect sizes over time
results_burned %>% 
  pivot_longer(
    cols = -c(days_before_fire, r_squared, adj_r_squared),
    names_to = "coefficient", values_to = "values") %>% 
  ggplot() + 
  geom_point(aes(x = days_before_fire, y = values)) + 
  labs(x = "Days before fire", y = "Effect size") +
  facet_wrap(~coefficient,scales="free",ncol = 2) + 
  theme_cowplot()

## check multicollinearity ----
mod1 <- lm(burn_severity ~ NDMI.d_prefire_10 + NDVI.d_prefire_10 + 
             elevation + slope + northness + eastness,
           data=final_df)
summary(mod1)
mc_result <- performance::multicollinearity(mod1)
print(mc_result)


## c. Bayesian implementation ----
n_iter <- 4000
n_thin <- ifelse(n_iter > 5000, 10, 1)

run_brm <- function(day,y_var, data) {
  ndvi_var <- paste0("NDVI.d_prefire_", day)
  ndmi_var <- paste0("NDMI.d_prefire_", day)
  
  predictors <- c(ndvi_var, ndmi_var, "elevation", "slope", "northness", "eastness", "doy")
  fixed_effects <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))
  
  model <- brm(
    bf(formula),
    family = gaussian(),
    data = data,
    chains = 4, 
    iter = n_iter, 
    thin = n_thin,
    cores = 4, seed = 1234,
    file = paste0(burn_severity_index,
                  "full_model_d_prefire",day))
  model_summary <- summary(model)
  
  coefs <- coef(model)
  
  tibble(
    days_before_fire = day,
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared
  ) %>%
    bind_cols(
      as_tibble_row(
        setNames(as.list(coefs),
                 c("Intercept","NDVI", "NDMI",
                   "elevation", "slope", "northness", "eastness",
                   "doy")))
    )
}



## d. Random forest ----
y_var <- "dnbr"
day <- 10
ndvi_var <- paste0("NDVI.d_prefire_", day)
ndmi_var <- paste0("NDMI.d_prefire_", day)

predictors <- c(ndvi_var, ndmi_var, "elevation", "slope", "northness", "eastness", "doy")
fixed_effects <- paste(predictors, collapse = " + ")
formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))

bay_qr_day10 <- brm(
  bf(formula, quantile = 0.75),
  data = final_df,
  family = asym_laplace()
)

summary(bay_qr_day10)