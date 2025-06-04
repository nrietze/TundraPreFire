library(terra)
library(tidyterra)
library(spatstat.utils)
library(tidyverse)
library(cowplot)
library(patchwork)
library(colorspace)
library(lubridate)
library(pbapply)
library(lme4)
library(lmerTest)
library(performance)
library(sjPlot)
library(tictoc)
library(mclust)
library(gt)
library(latex2exp)
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
           descals_burn_class = case_when(
             descals_burn_class == 99 ~ "nodata",
             TRUE ~ as.character(descals_burn_class)
           ),
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

rowwise_cumsum <- function(df, pattern) {
  prefix <- "cum_"
  # Select NDMI/NDVI column names
  cols <- df %>% 
    select(matches(pattern)) %>%
    select(order(as.integer(str_extract(names(.), "\\d+")))) %>%
    names()
  
  # Compute row-wise cumulative sums
  cumsum_df <- df %>%
    select(all_of(cols)) %>%
    t() %>%
    apply(2, cumsum) %>%
    t() %>%
    as.data.frame()
  
  # Assign new column names
  names(cumsum_df) <- paste0(prefix, cols)
  
  return(cumsum_df)
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
FONT_SIZE <- 18

frac_to_sample <- 0.01
frac_int <- frac_to_sample *100
SCALE_FACTOR <- ifelse(burn_severity_index == "dNBR", 1000, 1)

# Load lookup tables
processing_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017)

optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

df_optimal_raster_doy <- optimality_lut %>% 
  mutate(fireid = as.factor(fireid),
         raster_doy = yday(date_raster)) %>% 
  mutate(raster_doy_scaled = raster_doy / 366,
         raster_doy = as.factor(raster_doy)) %>% 
  select(raster_doy,raster_doy_scaled,fireid)

# Load fire perimeters
fire_perimeters <- vect(
  sprintf("%s/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg",DATA_DIR)
)

df_all_fire_perims <- fire_perimeters %>% as.data.frame() %>% 
  filter(tst_year >= 2017)

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
  
  df_allfires <- read_csv2(fname_model_data)
  
} else {
  df_allfires_fromfile <- subset_lut %>%
    split(seq(nrow(.))) %>%
    map_dfr(~ load_data(.x, burn_severity_index, frac_int,
                        return_df_only = TRUE)) %>%
    compact()
  
  df_allfires_fromfile$burn_doy <- yday(df_allfires_fromfile$burn_date)
  
  ndmi_cumsum <- rowwise_cumsum(df_allfires_fromfile, "^NDMI\\.d_prefire_\\d+$")
  ndvi_cumsum <- rowwise_cumsum(df_allfires_fromfile, "^NDVI\\.d_prefire_\\d+$")
  df_allfires <- bind_cols(df_allfires_fromfile, ndmi_cumsum, ndvi_cumsum)
  
  write_csv2(df_allfires,fname_model_data )
}

# scale variables
df_allfires <- df_allfires %>% 
  mutate(fireid = as.factor(fireid),
         doy = yday(date),
         # scale predictors
         burn_doy_scaled = burn_doy / 366,
         elevation = scale(elevation)[,1],
         across(contains("lst"), ~ as.numeric(scale(.)))) %>% 
  left_join(df_optimal_raster_doy,
            by = "fireid"
  )

# 2. Descriptive statistics ----
# ==============================.

## a. Summary tables of studied fire scars ----

### Figure & Table S1: Fire size ----
table_S1 <- topN_fires %>% as.data.frame() %>% 
  mutate(farea = 100 * farea, #convert to ha
         burn_date = ymd(paste(tst_year,tst_month,tst_day,sep="-"))) %>% 
  select(fireid,farea, burn_date) %>%
  gt() %>% 
  fmt_number(
    columns = farea,
    decimals = 0,
    sep_mark ="'"
  ) %>% 
  cols_label(
    fireid = "Fire Event ID",
    farea = "Burned area (ha)",
    burn_date = "Earliest burn date"
  ) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels());table_S1

if (SAVE_FIGURES){
  gtsave(table_S1,filename = "data/tables/Table_S1.html")
}

# Plot bar chart of fire sizes
burned_area_all <- sum(df_all_fire_perims$farea)
burned_area_top25 <- sum(topN_fires$farea)

pct_studied <- (burned_area_top25 / burned_area_all) * 100

s <- sprintf("The 25 studied fire scars made up %.1f%% \n of the total burned area in the study region since 2017.",round(pct_studied,1))

annot_xy <- df_all_fire_perims %>% 
  arrange(farea) %>% 
  filter(fireid %in% topN_fires$fireid) %>% 
  first() %>% 
  select(fireid, farea)

(p1 <- ggplot(df_all_fire_perims,aes(y = reorder(fireid,-farea), x = farea *100,
                                     fill = ifelse(fireid %in% topN_fires$fireid,"studied","not studied"))) +
    geom_bar(stat = 'identity') +
    geom_vline(aes(xintercept = mean(topN_fires$farea*100)),color = "grey40", lty = "dashed") +
    annotate("text", y = factor(annot_xy$fireid), x = annot_xy$farea + 100 *200,
             label = s,
             size = 7, color = "grey40") +
    # geom_curve(aes(x = log(1e4), y = .3 * max(counts),
    #                xend = log(900), yend = .03 * max(counts)),
    #            arrow = arrow(length = unit(0.08, "inch")), linewidth = 1,
    #            color = "grey40", curvature = -0.3) +
    labs(y = "Fire Event ID",x = "Burned area (ha)",
         fill = "Fire scar included in study?") +
    scale_y_discrete(labels = NULL, breaks = NULL) +
    theme_cowplot(FONT_SIZE) +
    theme(legend.position = "bottom"))

(p2 <- ggplot(topN_fires,aes(x = tst_year, y = farea *100)) +
    geom_bar(stat = 'identity', fill = "#730000") +
    labs(x = "Fire Year",
         y = "Total annual burned area of studied fires (ha)") +
    theme_cowplot(FONT_SIZE) +
    theme(legend.position = "bottom"))

if (SAVE_FIGURES){
  ggsave2(p2,filename = "figures/Figure_1b.png", bg ="white",
          width = 6, height = 7)
}  

fig_s1 <- p1 + p2

if (SAVE_FIGURES){
  ggsave2(fig_s1,filename = "figures/Figure_S1.png", bg ="white",
          width = 10, height = 8)
}

## b. Number of sample points per fire scar ----
df_subset %>% 
  group_by(fireid) %>% 
  distinct(ObservationID,.keep_all = TRUE) %>% 
  tally() %>% 
  ggplot() + 
  geom_bar(aes(x = reorder(fireid,-n),y = n),stat = "identity") + 
  labs(y = "Number of valid sample points", x = "Fire ID") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

if (SAVE_FIGURES){
  ggsave2(filename = "figures/N_sample_points_per_scar.png", bg ="white",
          width = 10, height = 6)
}

## c. Histogram of dNBR in sampled data ----
PlotSeverityHistogram <- function(df,burn_severity_index,
                                  burn_severity_index_label, 
                                  facet = TRUE, SAVEFIG =FALSE) {
  p <- df %>% 
    group_by(fireid) %>% 
    distinct(!!sym(burn_severity_index), .keep_all = TRUE) %>% 
    ggplot() +
    geom_histogram(aes(x = !!sym(burn_severity_index)),fill = "salmon") +
    labs(x = sprintf("Burn Severity (%s)",burn_severity_index_label),
         fill = "") +
    theme_cowplot()
  
  if (facet){
    p + facet_wrap(~fireid, scales = "free_y")
  }  
  
  if (SAVEFIG == TRUE){
    ggsave2(p,filename = sprintf("figures/Histograms_%s_top25_scars.png",burn_severity_index),
            width = 12,height = 10,bg = "white")
  }
  
  return(p)
}

burn_severity_index <- "dnbr"
(hist_dnbr <- PlotSeverityHistogram(df_subset,burn_severity_index,
                                    burn_severity_index_label = "dNBR", SAVEFIG = SAVE_FIGURES))

burn_severity_index <- "dnbr_corr"
(hist_dnbr_corr <- PlotSeverityHistogram(df_subset,burn_severity_index,
                                         burn_severity_index_label = "dNBR_corr", SAVEFIG = F))

burn_severity_index <- "rdnbr"
(hist_rdnbr <- PlotSeverityHistogram(df_subset,burn_severity_index,
                                burn_severity_index_label = "RdNBR", SAVEFIG = SAVE_FIGURES))

burn_severity_index <- "rdnbr_corr"
(hist_rdnbr_corr <- PlotSeverityHistogram(df_subset,burn_severity_index,
                                burn_severity_index_label = "RdNBR_corr", SAVEFIG = SAVE_FIGURES))
burn_severity_index <- "rbr"
(hist_rbr <- PlotSeverityHistogram(df_subset,burn_severity_index,
                                burn_severity_index_label = "RBR", SAVEFIG = SAVE_FIGURES))

burn_severity_index <- "dgemi"
(hist_dgemi <- PlotSeverityHistogram(df_subset,burn_severity_index,
                                burn_severity_index_label = "dGEMI", 
                                facet = F, SAVEFIG = F))


### i. Map & Histogram of burn severity data ----
pg <- PlotSeverityMapHist(cropped_severity_raster)
if (SAVE_FIGURES){
  ggsave2(pg,filename = sprintf("figures/%s_%s_map_hist.png",burn_severity_index,FIRE_ID),
          width = 8,height = 6,bg = "white")
}

### ii. Plot distribution of total raster vs. sample points ----
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

## d. Plot NDMI, NDVI & LST time series ----

### i. NDMI obs. & fitted  ----
obsid_per_fire <- df_subset %>%
  group_by(fireid) %>%
  summarize(random_obsid = sample(unique(ObservationID), 1), .groups = "drop")

df_subset_filtered <- df_subset %>% 
  inner_join(obsid_per_fire, by = "fireid") %>%
  filter(ObservationID == random_obsid) # picks a random ObservationID per fireid

df_subset_filtered_transformed <- df_subset_filtered %>% 
  pivot_longer(
    cols = starts_with("NDMI.d_prefire_"),
    names_to = "d_prefire_name",
    values_to = "NDMI_fit"
  ) %>%
  mutate(d_prefire = as.integer(str_extract(d_prefire_name, "\\d+"))) %>%
  filter(d_prefire >= 0 & d_prefire <= 30)

df_subset_filtered_obs <- df_subset_filtered %>%
  mutate(d_prefire = burn_doy - doy) %>%
  filter(d_prefire >= 0 & d_prefire <= 30)

ggplot() +
  geom_point(data = df_subset_filtered_obs, aes(x = d_prefire, y = DailyMeanNDMI), color = "blue") +
  geom_line(data = df_subset_filtered_transformed, aes(x = d_prefire, y = NDMI_fit), color = "red") +
  facet_wrap(~ fireid) +
  labs(
    x = "Days Before Fire", 
    y = "NDMI", 
    title = "NDMI Observations vs. Fitted"
  ) +
  scale_x_reverse() +
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(filename = "figures/sample_NDMI_fitted_and_observed.png",
          width = 8,height = 6,bg = "white")
}

### ii. Lineplots NDMI time series ----
obsid_per_fire <- df_subset %>%
  group_by(fireid) %>%
  summarize(random_obsid = sample(unique(ObservationID), 10), .groups = "drop")

df_subset_filtered <- df_subset %>% 
  inner_join(obsid_per_fire, by = "fireid") %>%
  filter(ObservationID == random_obsid)

df_subset_transformed <- df_subset_filtered %>% 
  pivot_longer(
    cols = starts_with("NDMI.d_prefire_"),
    names_to = "d_prefire_name",
    values_to = "NDMI_fit"
  ) %>%
  mutate(d_prefire = as.integer(str_extract(d_prefire_name, "\\d+"))) %>%
  filter(d_prefire >= 0 & d_prefire <= 30)


df_subset_transformed <- df_subset_filtered %>%
  pivot_longer(
    cols = matches("^(NDMI\\.d_prefire_|NDVI\\.d_prefire_|lst_pred_d_prefire_)"),
    names_to = c(".value", "d_prefire"),
    names_pattern = "(NDMI\\.d_prefire_|NDVI\\.d_prefire_|lst_pred_d_prefire_)(\\d+)"
  ) %>%
  mutate(d_prefire = as.integer(d_prefire)) %>%
  filter(d_prefire >= 0 & d_prefire <= 30) %>% 
  rename(
    NDMI_fit = NDMI.d_prefire_,
    NDVI_fit = NDVI.d_prefire_,
    LST_fit = lst_pred_d_prefire_
  )

ggplot(data = df_subset_transformed) +
  geom_line(aes(x = d_prefire, y = NDMI_fit,group = ObservationID), alpha = .1) +
  geom_vline(xintercept = 0,lty = "dashed") +
  facet_wrap(~ fireid) +
  labs(
    x = "Days Before Fire", 
    y = "NDMI", 
    title = "NDMI fits"
  ) +
  scale_x_reverse() +
  theme_cowplot()

ggplot(data = df_subset_transformed) +
  geom_line(aes(x = d_prefire, y = LST_fit,group = ObservationID), alpha = .1) +
  geom_vline(xintercept = 0,lty = "dashed") +
  facet_wrap(~ fireid) +
  labs(
    x = "Days Before Fire", 
    y = "LST (scaled)", 
    title = "LST fits"
  ) +
  scale_x_reverse() +
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(filename = "figures/sample_NDMI_timeseries.png",
          width = 8,height = 8,bg = "white")
}


### iii. Scatterplot NDMI vs. LST ----
ggplot(df_subset) +
  geom_point(aes(x = DailyMeanNDMI,y = LST,color = !!sym(burn_severity_index)),
             alpha = .5) +
  labs(x = "Daily mean NDMI",
       y = "Daily mean LST (scaled)",
       color = sprintf("Burn severity (%s)",burn_severity_index)) +
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(filename = sprintf("figures/Scatter_LST_NDMI_%s_%s.png",
                             burn_severity_index,FIRE_ID),
          width = 8,height = 6,bg = "white")
}

## e. Histograms of all predictors ----
predictors <- c("NDVI.d_prefire_", "NDMI.d_prefire_","cumsum_lst_d_prefire_", 
                "elevation", "slope", "northness", "eastness", "burn_doy_scaled")

predictor_pattern <- paste(predictors, collapse = "|")

# Omit columns where prefire day > 30
selected_cols <- df_subset %>%
  select(matches(predictor_pattern)) %>%
  select(where(~ {
    colname <- deparse(substitute(.))
    if (str_detect(colname, "_\\d+$")) {
      num <- as.integer(str_extract(colname, "\\d+$"))
      !is.na(num) && num <= 30
    } else {
      TRUE
    }
  }))

df_long <- selected_cols %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

p_hist <- ggplot(df_long, aes(x = value, fill = variable)) +
  geom_histogram(bins = 30, alpha = 0.7, color = "black") +
  facet_wrap(~ variable, scales = "free") +
  theme_cowplot() +
  labs(title = "Histograms of Predictor Variables", x = "Value", y = "Count") +
  theme(legend.position = "none")

if (SAVE_FIGURES){
  ggsave2(p_hist,filename = "figures/Histogram_predictors.png",
          width = 12,height = 12,bg = "white")
}


## e. dNBR sensitivity over all dates (for largest scar) ----
UTM_TILE_ID <- "54WXE"
year <- 2020
selected_perim <- fire_perimeters %>% filter(fireid == 14211)

severity_raster_list <- list.files(path = paste0(HLS_DIR,"severity_rasters"),
                                   pattern = sprintf("^%s_%s_%s.*\\.tif$",
                                                     burn_severity_index,UTM_TILE_ID,year),
                                   full.names = TRUE)

all_dnbr_rast <- rast(severity_raster_list)

all_dnbr_masked <- mask(all_dnbr_rast,selected_perim)

names(all_dnbr_rast) <- raster_names

# Convert raster to data frame
dnbr_df <- as.data.frame(all_dnbr_rast, xy = FALSE, na.rm = TRUE)

# Reshape to long format
dnbr_long <- dnbr_df %>%
  pivot_longer(cols = everything(), names_to = "date", values_to = "dNBR") %>% 
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

## f. Distributions of cumulative NDMI and NDVI ----
plot_distributions <- function(data, varname, cmap){
  df_plot <- data %>%
    select(contains(paste0(varname,"."))) %>%
    pivot_longer(cols = everything(), names_to = "var", values_to = "value") %>% 
    mutate(var_day = as.integer(str_extract(var, "\\d+$"))) %>%
    mutate(var = fct_reorder(var, var_day))
  
  xlims <- if(varname == "NDMI") c(-5,10) else c(0,15)
  
  ggplot(df_plot,aes(x = value, group = var_day)) +
    geom_density(aes(fill = var_day),alpha = 0.2, size = 0.1) +  
    scale_fill_continuous_sequential(palette = cmap) +
    xlim(xlims) +
    theme_cowplot() +
    labs(x = sprintf("Cumulative %s",varname), y = "Density",
         fill = "Sum period (days before fire)",
         title = sprintf("Distribution of %s for each day pre-fire",varname))
}

(fig <- plot_distributions(df_subset,"NDMI",cmap = "Mako"))

if (SAVE_FIGURES){
  ggsave2(fig ,filename = "figures/cumNDMI_distribution_top25_scars.png",
          width = 6,height = 8,bg = "white")
}

(fig <- plot_distributions(df_subset,"NDVI",cmap = "GreenYellow"))

if (SAVE_FIGURES){
  ggsave2(fig ,filename = "figures/cumNDVI_distribution_top25_scars.png",
          width = 6,height = 8,bg = "white")
}

## g. Check spatial autocorrelation ----
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


## h. cluster plots ----
burn_severity_index <- "dnbr"

df_allfires %>% 
  mutate(avgNDMI = rowMeans(across(c(NDMI.d_prefire_1, NDMI.d_prefire_2, NDMI.d_prefire_3)), na.rm = TRUE),
         avgNDVI = rowMeans(across(c(NDVI.d_prefire_1, NDVI.d_prefire_2, NDVI.d_prefire_3)), na.rm = TRUE)) %>% 
  ggplot() + 
  geom_point(aes(x = elevation,y =  avgNDMI, color = !!sym(burn_severity_index)),alpha = .2) + 
  theme_cowplot()

ggsave2(sprintf("figures/cluster_ndmi_slope_%s.png",burn_severity_index),
        bg = "white",width = 10, height = 10)

## i. Plot time series data ----

# Plot daily NMDI sooths over days before fire
df_transformed <- df_allfires %>%
  pivot_longer(
    cols = matches("NDMI_prefire_d\\.|NDVI_prefire_d\\.|cumsum_lst_d_prefire_"),
    names_to = c(".value", "days_before_fire"),
    names_pattern = "(NDMI_prefire_d|NDVI_prefire_d|cumsum_lst_d_prefire)\\.?([0-9]+)",
    values_drop_na = TRUE
  ) %>%
  rename(
    NDMI = NDMI_prefire_d,
    NDVI = NDVI_prefire_d,
    LST = cumsum_lst_d_prefire
  ) %>% 
  mutate(
    # Convert dates to day-of-year (doy)
    date_doy = yday(date),
    burn_doy = yday(burn_date),
    
    # Calculate days before fire (negative = pre-fire)
    days_before_fire = burn_doy - date_doy
  ) %>%
  # Keep only pre-fire observations (optional)
  filter(days_before_fire > 0)

fig <- ggplot(df_transformed,aes(x = days_before_fire, y = DailyMeanNDMI)) +
  stat_smooth(aes(group = ObservationID),geom="line",
              alpha=0.1, linewidth=.5, span=0.5) +
  # geom_point(alpha=0.1) +
  facet_wrap(~fireid) +
  labs(
    x = "Days Before Fire",
    y = "Daily Mean NDMI",
    title = "NDMI Trend Before Fire"
  ) +
  lims(y = c(-1.2,1.2)) +
  scale_x_reverse(limits = c(30,0)) +
  theme_cowplot()

ggsave2(fig,"figures/ndmi_points_per_fire_q75_20pct.png",
        bg = "white",width = 10, height = 10)

# LST smooths pre-fire
ggplot(df_transformed, aes(x = days_before_fire, y = LST)) +
  stat_smooth(aes(group = ObservationID),geom="line", 
              alpha=0.1, linewidth=.5, span=0.5) +
  facet_wrap(~fireid) +
  labs(
    x = "Days Before Fire",
    y = "Daily LST",
    title = "LST Before Fire"
  ) +
  scale_x_reverse() +
  theme_cowplot()

ggsave2("figures/lst_curves_per_fire_alldnbr.png",
        bg = "white",width = 10, height = 10)

# Plot nr. of non-NA NDMI observations for each day in 25 largest fires
df_allfires %>% 
  group_by(date,fireid) %>% 
  summarise(non_na_count = sum(!is.na(DailyMeanNDMI)),
            burn_date = first(burn_date)) %>% 
  ggplot(aes(x = yday(date), y = non_na_count)) +
  geom_col(fill = "steelblue", width = 0.7) + 
  geom_vline(aes(xintercept = yday(burn_date))) +
  labs(
    x = "Date",
    y = "Non-NA Observations",
    title = "Daily Mean NDMI Data Availability"
  ) +
  facet_wrap(~fireid,ncol = 2,scales = "free") +
  theme_cowplot()

ggsave2("figures/n_obs_ndmi_per_fire_scar.png",
        bg = "white",width = 8, height = 16)

# Plot nr. of non-NA LST observations for each day in 25 largest fires 
df_allfires %>% 
  group_by(date,fireid) %>% 
  summarise(non_na_count = sum(!is.na(LST)),
            burn_date = first(burn_date)) %>% 
  ggplot(aes(x = yday(date), y = non_na_count)) +
  geom_col(fill = "steelblue", width = 0.7) + 
  geom_vline(aes(xintercept = yday(burn_date))) +
  labs(
    x = "Date",
    y = "Non-NA Observations",
    title = "Daily LST Data Availability"
  ) +
  facet_wrap(~fireid,ncol = 2,scales = "free") +
  theme_cowplot()

ggsave2("figures/n_obs_lst_per_fire_scar.png",
        bg = "white",width = 8, height = 16)

# Plot all cumulated NDMI before fire
df_long <- df_subset %>%
  sample_frac(.01) %>% 
  select(
    starts_with("NDMI.d_prefire_"),
    any_of(c("day", "dnbr", "rdnbr", "dgemi", "rbr"))
  ) %>% 
  pivot_longer(
    cols = starts_with("NDMI.d_prefire_"),
    names_to = "day",
    names_prefix = "NDMI.d_prefire_",
    values_to = "NDMI_value"
  ) %>%
  mutate(day = as.numeric(day))

ggplot(df_long, aes(x = day, y = NDMI_value)) +
  geom_point(aes(color = dnbr),alpha = 0.02) +
  theme_cowplot()

# Plot time-integrated NDMI 7 days before fire vs. dNBR
xvar <- "NDMI.d_prefire_7"

df_all_subset <- df_allfires %>% 
  mutate(severity_group = if_else(!!sym(burn_severity_index) >= 0.245,
                                  "burned (≥ 0.245)",
                                  "unburned (< 0.245)"))

(fig <- ggplot(df_all_subset, aes(x = !!sym(xvar),
                                  y = !!sym(burn_severity_index),
                                  color = severity_group)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = 'lm', se = FALSE) +
    labs(x = sprintf("Cumulative %s (sum over 1 to 7 days pre-fire)", index_name),
         y = burn_severity_index,
         color = "Burn severity class \n (threshold from Kolden et al.(2012))") +
    theme_cowplot())

if (SAVE_FIGURES){
  ggsave2(fig,
          filename = sprintf("figures/%s_%s_%s.png",
                             burn_severity_index,xvar,FIRE_ID),
          width = 8,height = 6,bg = "white")
}

# Plot time-integrated LST 7 days before fire vs. dNBR
df_all_subset %>% sample_n(1e5) %>% 
  filter(descals_burn_class != "nodata") %>% 
  ggplot( aes(x = cumsum_lst_d_prefire_7 ,
              y = !!sym(burn_severity_index),
              color = descals_burn_class)) +
  geom_point(alpha = .2) +
  geom_smooth(method='lm') +
  labs(x = "Cumulative LST (sum over 1 to 7 days pre-fire)",
       y = burn_severity_index,
       color = "Burn class (Descals et al., 2022)") +
  theme_cowplot()

# 3. Cluster data ----
# ====================.
burn_severity_index <- "dnbr_corr"

df_clust_subset <- df_allfires %>% 
  sample_frac(.05) %>% 
  drop_na(!!sym(burn_severity_index))

ggplot(df_clust_subset) +
  geom_histogram(aes(x = !!sym(burn_severity_index)),fill = "salmon") +
  theme_cowplot()

ggplot(df_clust_subset) +
  geom_boxplot(aes(y = !!sym(burn_severity_index), x = as.factor(raster_doy))) +
  theme_cowplot()

ggplot(df_clust_subset) +
  geom_boxplot(aes(y = !!sym(burn_severity_index), x = as.factor(raster_doy))) +
  theme_cowplot()

df_clust_subset %>% 
  ggplot() +
  geom_point(aes(x = burn_doy, y = raster_doy)) +
  lims(x = c(170,270), y = c(170,270)) +
  theme_cowplot()

# mcl <- Mclust(df_clust_subset[c('NDMI.d_prefire_30','elevation','burn_doy')])
mcl <- Mclust(df_allfires %>%
                drop_na(!!sym(burn_severity_index)) %>% 
                select(burn_severity_index,raster_doy))

summary(mcl)
plot(mcl,  what = "classification")
plot(mcl,  what = "uncertainty")

df_mcl <- data.frame(mcl$BIC[], G = 1:nrow(mcl$BIC)) %>% 
  pivot_longer(cols = 1:14, names_to = "Model", values_to = "BIC") %>% 
  mutate(Model = factor(Model, levels = mclust.options("emModelNames")))

ggplot(df_mcl, aes(x = G, y = BIC, colour = Model, shape = Model)) +
  geom_point() + 
  geom_line() +
  scale_shape_manual(values = mclust.options("bicPlotSymbols")) +
  scale_color_manual(values = mclust.options("bicPlotColors")) +
  scale_x_continuous(breaks = unique(df_mcl$G)) +
  xlab("Number of mixture components") +
  guides(shape = guide_legend(ncol=2)) +
  theme_cowplot()

ggsave("figures/clustering/BIC_5pct_NDMI30_elev_doyburn.png",bg = 'white',
       height = 6, width = 10)


# report uncertainty statistics for clusters
df_uncertainty <- data.frame(cluster = mcl$classification, 
                             uncertainty = mcl$uncertainty)

# Summarize stats per cluster
summary_table <- df_uncertainty %>%
  group_by(cluster) %>%
  summarise(
    mean = mean(uncertainty),
    min = min(uncertainty),
    max = max(uncertainty),
    p95 = quantile(uncertainty, 0.95)
  ) %>%
  pivot_longer(-cluster, names_to = "statistic", values_to = "value") %>%
  pivot_wider(names_from = cluster, names_prefix = "Cluster_", values_from = value) %>%
  relocate(statistic)

# Convert to gt table
gt_table <- summary_table %>%
  gt() %>%
  tab_header(
    title = "Uncertainty Summary Statistics per Cluster"
  ) %>%
  fmt_number(
    columns = starts_with("Cluster_"),
    decimals = 3
  )

gt_table

# 4. Run linear model ----
# ========================.
## a. mixed model ----

# Function to run lm and extract statistics
run_lm_day <- function(day,y_var, data) {
  ndvi_var <- paste0("NDVI.d_prefire_", day)
  ndmi_var <- paste0("NDMI.d_prefire_", day)
  lst_var <- paste0("lst_pred_d_prefire_", day)
  
  predictors <- c(ndvi_var, ndmi_var,lst_var, 
                  "elevation", "slope", "northness", "eastness",
                  "burn_doy_scaled")
  fixed_effects <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)+ (1 | raster_doy)"))
  
  model <- lmerTest::lmer(formula, data = data, REML = FALSE)
  model_summary <- summary(model)
  
  # Extract fixed effect estimates
  coefs <- fixef(model)  # Named vector
  std_err <- coef(summary(model))[, "Std. Error"]
  t_vals <- coef(summary(model))[, "t value"]
  p_vals <- coef(summary(model))[, "Pr(>|t|)"]
  
  # R-squared using performance package (optional but helpful)
  r2_vals <- performance::r2(model)
  
  coefs_tbl <- tibble(
    predictor = names(coefs),
    coef = unlist(coefs),
    std_err = unlist(std_err),
    t_val = unlist(t_vals),
    p_val = unlist(p_vals)
  ) %>%
    mutate(
      days_before_fire = day,
      r_squared_marg = r2_vals$R2_marginal,
      r_squared_cond = r2_vals$R2_conditional,
    )
  
  # Random effects
  ranef_vals <- ranef(model, condVar = TRUE)
  re_list <- lapply(names(ranef_vals), function(grp) {
    df <- as.data.frame(ranef_vals[[grp]])
    df$level <- rownames(ranef_vals[[grp]])
    
    # Extract conditional variances for CI
    post_var <- attr(ranef_vals[[grp]], "postVar")[1, 1, ]
    post_sd <- sqrt(post_var)
    df$lower <- df$`(Intercept)` - 1.96 * post_sd
    df$upper <- df$`(Intercept)` + 1.96 * post_sd
    df$significant <- df$lower * df$upper > 0  # True if CI does not include 0
    
    df <- df %>%
      select(level, `(Intercept)`, lower, upper, significant) %>%
      rename(!!paste0("day_", day) := `(Intercept)`)
    
    return(df)
  })
  names(re_list) <- names(ranef_vals)
  
  return(list(
    fixed_effects = coefs_tbl,
    random_effects = re_list
  ))
}

burn_severity_index <- "dgemi"

# subset of upper 75th percentile
df_subset <- df_allfires %>% 
  group_by(fireid) %>% 
  # filter(dnbr >= 0.245) %>%  # (threshold from Kolden et al. 2012)
  filter(quantile(!!sym(burn_severity_index), 0.75,na.rm=T)<!!sym(burn_severity_index)) %>% 
  ungroup()

# only model dNBR above burned threshold
model_results_by_day <- pblapply(1:30, function(day) {
  run_lm_day(y_var = burn_severity_index,
             day = day,
             data = df_subset)
})

# combine fixed effects
fixed_effects_all <- bind_rows(lapply(model_results_by_day, `[[`, "fixed_effects"))

fname_results_rds <- sprintf("data/models/results_lmer_%s_nocumsum_fireid_rastdoy.rds",burn_severity_index)

if (file.exists(fname_results_rds)){
  results <- readRDS(fname_results_rds)
} else{
  results <- fixed_effects_all %>% arrange(days_before_fire)
  saveRDS(results,fname_results_rds)
}

results_clean <- results %>%
  mutate(
    # Extract base predictor: NDVI, NDMI, LST (or keep as-is for Intercept, etc.)
    base_predictor = case_when(
      str_detect(predictor, "^NDVI\\.d_prefire_") ~ "NDVI",
      str_detect(predictor, "^NDMI\\.d_prefire_") ~ "NDMI",
      str_detect(predictor, "^lst_pred_d_prefire_") ~ "LST",
      TRUE ~ predictor  # Keep Intercept or others as-is
    ),
    
    # Extract numeric day from "d_prefire_" part
    day_from_name = as.integer(str_extract(predictor, "(?<=d_prefire_)[0-9]+")),
    
    # Replace with parsed day if available
    days_before_fire = coalesce(day_from_name, days_before_fire)
  ) %>%
  select(-day_from_name) 

# Combine random effects separately by group
combine_random_effects <- function(group) {
  re_group_list <- lapply(model_results_by_day, function(res) res$random_effects[[group]])
  full_df <- reduce(re_group_list, full_join, by = "level")
  return(full_df)
}

fireid_re_df <- combine_random_effects("fireid")
raster_doy_re_df <- combine_random_effects("raster_doy")

## b. Plot model results ----

### i. R2 curves ----
max_r2_day <- results_clean %>%
  filter(r_squared_marg == max(r_squared_marg, na.rm = TRUE)) %>%
  select(days_before_fire, r_squared_marg) %>% 
  distinct()

max_day <- max_r2_day$days_before_fire
max_r2 <- max_r2_day$r_squared_marg

s <- sprintf("Maximum marginal R² is \n %s days before the fire", max_day)

ylims <- c(0, 0.25)

(p_r2_marg <- ggplot(results) + 
    geom_point(aes(x = days_before_fire, y = r_squared_marg)) +
    geom_vline(aes(xintercept = max_day), lty = "dashed") +
    annotate("text", x = max_day + 11, y = ylims[2] - .1,
             label = s,
             size = 5, color = "grey40") +
    geom_curve(aes(x = max_day + 13, y = ylims[2] - .11, 
                   xend = max_day, yend = max_r2 -0.005),
               arrow = arrow(length = unit(0.08, "inch")), linewidth = 1,
               color = "grey40", curvature = 0.3) +
    scale_x_reverse() +
    scale_y_continuous(breaks = seq(ylims[1],ylims[2],.05),
                       labels = seq(ylims[1],ylims[2],.05),
                       limits = ylims
    ) +
    labs(x = "Days before fire", y = expression("marginal"~R^2)) + 
    theme_cowplot(FONT_SIZE) )

if (SAVE_FIGURES){
  ggsave2(p_r2_marg,filename = sprintf("figures/Figure_3a_%s.png",burn_severity_index),
          width = 8,height = 8,bg = "white")
}

max_r2_cond_day <- results_clean %>%
  filter(r_squared_cond == max(r_squared_cond, na.rm = TRUE)) %>%
  select(days_before_fire, r_squared_cond) %>% 
  distinct()

max_cond_day <- max_r2_cond_day$days_before_fire
max_cond_r2 <- max_r2_cond_day$r_squared_cond

s_cond <- sprintf("Maximum conditional R² is \n %s days before the fire", max_cond_day)

ylims <- c(0.6, 0.8)

(p_r2_cond <- ggplot(results) + 
    geom_point(aes(x = days_before_fire, y = r_squared_cond)) +
    geom_vline(aes(xintercept = max_cond_day), lty = "dashed") +
    annotate("text", x = max_cond_day + 8, y = max_cond_r2 + .05,
             label = s_cond,
             size = 5, color = "grey40") +
    geom_curve(aes(x = max_cond_day + 8, y = max_cond_r2 + .035, 
                   xend = max_cond_day, yend = max_cond_r2 + 0.005),
               arrow = arrow(length = unit(0.08, "inch")), linewidth = 1,
               color = "grey40", curvature = 0.3) +
    scale_x_reverse() +
    labs(x = "Days before fire", y = "conditional R²") + 
    scale_y_continuous(breaks = seq(ylims[1],ylims[2],.05),
                       labels = seq(ylims[1],ylims[2],.05),
                       limits = ylims
                       ) +
    theme_cowplot(FONT_SIZE) )

if (SAVE_FIGURES){
  ggsave2(p_r2_cond,filename = sprintf("figures/Figure_3b_%s.png",burn_severity_index),
          width = 8,height = 8,bg = "white")
}

p_r2_marg / p_r2_cond + 
  plot_annotation(tag_levels = "a", tag_suffix = ')')

if (SAVE_FIGURES){
  ggsave2(filename = sprintf("figures/Figure_3_%s.png",burn_severity_index),
          width = 8,height = 8,bg = "white")
}

# find maxima of each predictor
ndmi_max_d <- results_clean %>%
  filter(base_predictor == "NDMI") %>%
  filter(coef == max(coef, na.rm = TRUE)) %>%
  select(base_predictor, days_before_fire, coef) %>%
  ungroup()
ndvi_max_d <- results_clean %>%
  filter(base_predictor == "NDVI") %>%
  filter(coef == max(coef, na.rm = TRUE)) %>%
  select(base_predictor, days_before_fire, coef) %>%
  ungroup()

### ii.  Plot effect sizes over time ----
mod_labs <- c(
  burn_doy_scaled = "Date of burn \n(Day of year)",
  NDVI = "Vegetation greenness \n(NDVI HLS)", 
  NDMI = "Vegetation moisture \n(NDMI HLS)",
  LST = "Land surface temperature \n(LST Landsat-8)",
  raster_doy_scaled = "Season of burn severity raster \n(Day of year)",
  elevation = "Elevation",
  slope = "Terrain slope",
  northness = "Northness",
  eastness = "Eastness"
)

results_clean %>%
  mutate(
    base_predictor = factor(base_predictor, levels = names(mod_labs)),
    significant = p_val < 0.05,
    alpha_val = ifelse(significant, 1, 0.5)
  ) %>%
  filter(base_predictor %in% c("NDVI", "NDMI", "LST",
                               "elevation", "slope", "northness", "eastness", 
                               "burn_doy_scaled")) %>%
  ggplot(aes(
    x = days_before_fire,
    y = coef,
    color = base_predictor,
    alpha = alpha_val,
    group = base_predictor
  )) +
  geom_point() +
  geom_errorbar(aes(
    ymin = coef - std_err,
    ymax = coef + std_err
  ), width = 0.5) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_x_reverse() +
  facet_wrap(~ base_predictor,
             scales = "free_y", 
             labeller = as_labeller(mod_labs),
             nrow = 2) +
  labs(x = "Days before fire", y = "Effect size (coefficient)") +
  scale_alpha_identity() +
  theme_cowplot(FONT_SIZE) +
  theme(legend.position = "none")

if (SAVE_FIGURES){
  ggsave2(filename = sprintf("figures/Figure_2_%s.png",burn_severity_index),
          width = 14,height = 8,bg = "white")
}

### Random effects ----
cols <- colnames(fireid_re_df)

# Find positions of 'day_*' and 'significant_*' columns based on structure
effect_cols <- cols[seq(2, ncol(fireid_re_df), by = 4)]
signif_cols <- cols[seq(5, ncol(fireid_re_df), by = 4)]

# Get the days from the names (e.g., from day_1, day_2)
day_nums <- as.integer(gsub("day_", "", effect_cols))

# Pivot effects
effects_long <- fireid_re_df %>%
  select(level, all_of(effect_cols)) %>%
  pivot_longer(
    cols = -level,
    names_to = "day_col",
    values_to = "effect"
  ) %>%
  mutate(day = as.integer(gsub("day_", "", day_col))) %>%
  select(-day_col)

# Pivot significance
signif_long <- fireid_re_df %>%
  select(level, all_of(signif_cols)) %>%
  pivot_longer(
    cols = -level,
    names_to = "signif_col",
    values_to = "significant"
  ) %>%
  mutate(day = rep(day_nums,25)) %>%
  select(-signif_col)

# Join the two long tables
fireid_plot_df <- left_join(effects_long, signif_long, by = c("level", "day"))

cols <- colnames(raster_doy_re_df)

# Find positions of 'day_*' and 'significant_*' columns based on structure
effect_cols <- cols[seq(2, ncol(raster_doy_re_df), by = 4)]
signif_cols <- cols[seq(5, ncol(raster_doy_re_df), by = 4)]

# Get the days from the names (e.g., from day_1, day_2)
day_nums <- as.integer(gsub("day_", "", effect_cols))

# Pivot effects
effects_long <- raster_doy_re_df %>%
  select(level, all_of(effect_cols)) %>%
  pivot_longer(
    cols = -level,
    names_to = "day_col",
    values_to = "effect"
  ) %>%
  mutate(day = as.integer(gsub("day_", "", day_col))) %>%
  select(-day_col)

# Pivot significance
signif_long <- raster_doy_re_df %>%
  select(level, all_of(signif_cols)) %>%
  pivot_longer(
    cols = -level,
    names_to = "signif_col",
    values_to = "significant"
  ) %>%
  mutate(day = rep(day_nums,20)) %>%
  select(-signif_col)

# Join the two long tables
raster_doy_plot_df <- left_join(effects_long, signif_long, by = c("level", "day"))

plot_random_effects <- function(plot_df, title) {
  ggplot(plot_df, aes(x = day, y = level)) +
    # Plot non-significant points in gray and small
    geom_point(
      data = filter(plot_df, significant == FALSE),
      color = "gray70",
      size = 1
    ) +
    # Plot significant points with effect-size coloring and larger size
    geom_point(
      data = filter(plot_df, significant == TRUE),
      aes(color = effect),
      size = 3
    ) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = "Effect size"
    ) +
    scale_x_reverse() +
    labs(
      title = title,
      x = "Days Before Fire",
      y = "Random Effect Level"
    ) +
    theme_cowplot() +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    )
}

plot_random_effects(fireid_plot_df, "Random Effect (fireid) Over Time")
if (SAVE_FIGURES){
  ggsave2(filename = "figures/Figure_S55.png",
          width = 8,height = 8,bg = "white")
}

plot_random_effects(raster_doy_plot_df, "Random Effect (raster_doy) Over Time")
if (SAVE_FIGURES){
  ggsave2(filename = "figures/Figure_S56.png",
          width = 8,height = 8,bg = "white")
}

## c. Write results to tables ----
make_predictor_table <- function(predictor_name) {
  results_clean %>%
    filter(base_predictor == predictor_name) %>%
    dplyr::select(days_before_fire, coef, std_err, t_val, p_val) %>%
    gt() %>%
    tab_header(
      title = md(glue::glue("**Results for Predictor: {predictor_name}**"))
    ) %>%
    cols_label(
      days_before_fire = md("**Days Before Fire**"),
      coef = md("**Coefficient**"),
      std_err = md("**Std. Error**"),
      t_val = md("**t Value**"),
      p_val = md("**p Value**")
    ) %>% 
    fmt_number(
      columns = c(coef, std_err, t_val, p_val),
      n_sigfig = 2
    ) %>%
    fmt_number(
      columns = days_before_fire,
      decimals = 0
    ) %>% 
    fmt(
      columns = p_val,
      fns = function(x) {
        ifelse(x < 0.05, "< 0.05", formatC(x, format = "f", digits = 2))
      }
    )
}

make_predictor_table("burn_doy_scaled")

unique_preds <- unique(results_clean$base_predictor)

tables_list <- lapply(unique_preds, make_predictor_table)
names(tables_list) <- unique_preds

lapply(seq_along(tables_list), function(i) {
  gtsave(
    tables_list[[i]],
    filename = paste0(burn_severity_index,"_predictor_table_", unique_preds[i], ".html"),
    path = "data/tables/"
  )
})

# cretae R2 table
r2_table <- results_clean %>%
  dplyr::select(base_predictor, days_before_fire, r_squared_marg, r_squared_cond) %>%
  distinct() %>%
  arrange(base_predictor, days_before_fire) %>%
  gt(groupname_col = "base_predictor") %>%
  cols_label(
    days_before_fire = md("**Days Before Fire**"),
    r_squared_marg = md("**Marginal R²**"),
    r_squared_cond = md("**Conditional R²**")
  ) %>%
  tab_header(title = md("**R² Summary by Predictor and Days Before Fire**")) %>% 
  fmt_number(decimals = 2, )


