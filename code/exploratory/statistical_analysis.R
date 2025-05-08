library(terra)
library(tidyterra)
library(spatstat.utils)
library(tidyverse)
library(cowplot)
library(patchwork)
library(colorspace)
library(lubridate)
library(pbapply)
library(INLA)
library(lme4)
library(brms)
library(tidybayes)
library(tictoc)
library(mclust)
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
  final_df <- subset_lut %>%
    split(seq(nrow(.))) %>%
    map_dfr(~ load_data(.x, burn_severity_index, frac_int,
                    return_df_only = TRUE)) %>%
    compact()
  
  final_df$burn_doy <- yday(final_df$burn_date)
  
  write_csv2(final_df,fname_model_data )
}

# scale variables
final_df <- final_df %>% 
  mutate(fireid = as.factor(fireid),
         # scale predictors
         burn_doy = burn_doy / 366,
         elevation = scale(elevation)[,1],
         across(contains("lst"), ~ as.numeric(scale(.))))

# subset of upper 75th percentile
df_subset <- final_df %>% 
  group_by(fireid) %>% 
  # filter(dnbr >= 0.245) %>%  # (threshold from Kolden et al. 2012)
  filter(quantile(dnbr, 0.75,na.rm=T)<dnbr) %>% 
  ungroup()

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
burn_severity_index <- "rdnbr_corr"

ggplot(data = df_subset) +
    geom_histogram(aes(x = !!sym(burn_severity_index)),binwidth = .005) +
    labs(x = sprintf("Burn Severity (%s)",burn_severity_index),
         fill = "",
         title = "Histograms per fire scar") +
    theme_cowplot()

(pdnbr <- ggplot(data = df_subset) +
  geom_histogram(aes(x = dnbr),binwidth = .005) +
  facet_wrap(~fireid) +
  labs(x = "Burn Severity (dNBR)",
       fill = "",
       title = "Histograms per fire scar") +
  theme_cowplot())

ggsave2(filename = "figures/Histograms_dnbr_top25_scars.png",
        width = 12,height = 10,bg = "white")


(prbr <- ggplot(data = final_df) +
    geom_histogram(aes(x = rbr),binwidth = .005) +
    facet_wrap(~fireid) +
    labs(x = "Burn Severity (RBR)",
         fill = "",
         title = "Histograms per fire scar") +
    theme_cowplot())

ggsave2(filename = "figures/Histograms_rbr_top25_scars.png",
        width = 12,height = 10,bg = "white")


pdgemi <- ggplot(data = final_df) +
  geom_histogram(aes(x = dgemi),binwidth = .005) +
  facet_wrap(~fireid) +
  labs(x = "Burn Severity (dgemi)",
       title = "Histograms per fire scar",
       fill = "") +
  theme_cowplot()

ggsave2(pdgemi,filename = "figures/Histograms_dgemi_top25_scars.png",
        width = 12,height = 10,bg = "white")

prdnbr <- ggplot(data = final_df) +
  geom_histogram(aes(x = rdnbr),binwidth = .005) +
  facet_wrap(~fireid) +
  labs(x = "Burn Severity (RdNBR)",
       title = "Histograms per fire scar",
       fill = "") +
  theme_cowplot()

ggsave2(prdnbr,filename = "figures/Histograms_rdnbr_top25_scars.png",
        width = 12,height = 10,bg = "white")


(fig <- pdnbr / pdgemi +
  plot_layout(axis_titles = "collect"))

if (SAVE_FIGURES){
  ggsave2(fig,filename = "figures/Histograms_dnbrvs_dgemi_top25_scars.png",
          width = 12,height = 10,bg = "white")
}

# overall distribution
ggplot(data = final_df) +
  geom_density(aes(x = dnbr),color = 'tomato', fill = 'salmon',alpha = 0.7) +
  labs(x = "Burn Severity (dNBR)",
       title = "dNBR distribution of the 25 largest tundra fires") +
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(filename = "figures/dNBR_distribution_top25_scars.png",
          width = 8,height = 6,bg = "white")
}

## Distributions of cumulative NDMI and NDVI ----
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


## cluster plots ----
burn_severity_index <- "dnbr"

final_df %>% 
  mutate(avgNDMI = rowMeans(across(c(NDMI.d_prefire_1, NDMI.d_prefire_2, NDMI.d_prefire_3)), na.rm = TRUE),
         avgNDVI = rowMeans(across(c(NDVI.d_prefire_1, NDVI.d_prefire_2, NDVI.d_prefire_3)), na.rm = TRUE)) %>% 
  ggplot() + 
  geom_point(aes(x = elevation,y =  avgNDMI, color = !!sym(burn_severity_index)),alpha = .2) + 
  theme_cowplot()

ggsave2(sprintf("figures/cluster_ndmi_slope_%s.png",burn_severity_index),
        bg = "white",width = 10, height = 10)

## b. Plot time series data ----

# Plot daily NMDI sooths over days before fire
df_transformed <- final_df %>%
  mutate(
    # Convert dates to day-of-year (doy)
    date_doy = yday(date),
    burn_doy = yday(burn_date),
    
    # Calculate days before fire (negative = pre-fire)
    days_before_fire = burn_doy - date_doy
  ) %>%
  # Keep only pre-fire observations (optional)
  filter(days_before_fire > 0)

fig <- ggplot(df_transformed, aes(x = days_before_fire, y = DailyMeanNDMI)) +
  stat_smooth(aes(group = ObservationID),geom="line", 
               alpha=0.1, linewidth=.5, span=0.5) +
  facet_wrap(~fireid) +
  labs(
    x = "Days Before Fire",
    y = "Daily Mean NDMI",
    title = "NDMI Trend Before Fire"
  ) +
  ylim(c(-1.2,1.2)) +
  scale_x_reverse() +
  theme_cowplot()

ggsave2(fig,"figures/ndmi_curves_per_fire_alldnbr_v2.png",
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
final_df %>% 
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

ggsave2("figures/n_obs_ndmi_per_fire_scar_v3.png",
        bg = "white",width = 8, height = 16)

# Plot nr. of non-NA LST observations for each day in 25 largest fires 
final_df %>% 
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

ggsave2("figures/n_obs_lst_per_fire_scar_v3.png",
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

df_all_subset <- final_df %>% 
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

# 3. Cluster data ----
# ====================.
df_clust_subset <- final_df %>% 
  sample_frac(.05)

mcl <- Mclust(df_clust_subset[c('NDMI.d_prefire_30','elevation','burn_doy')])

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
burn_severity_index <- "rdnbr_corr"

## a. mixed model ----

# Function to run lm and extract statistics
run_lm_day <- function(day,y_var, data) {
  ndvi_var <- paste0("NDVI.d_prefire_", day)
  ndmi_var <- paste0("NDMI.d_prefire_", day)
  lst_var <- paste0("cumsum_lst_d_prefire_", day)
  
  predictors <- c(ndvi_var, ndmi_var,lst_var, "elevation", "slope", "northness", "eastness", "burn_doy")
  fixed_effects <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))
  
  model <- lmer(formula, data = data, REML = FALSE)
  model_summary <- summary(model)
  
  # reduced_formula <- as.formula(paste(y_var, "~", 1, "+ (1 | fireid)"))
  # reduced.lmer <- lmer(reduced_formula,data = data, REML = FALSE)
  # 
  # anova(reduced.lmer, model)  # the two models are not significantly different
  
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
                 c("Intercept","NDVI", "NDMI","LST",
                   "elevation", "slope", "northness", "eastness",
                   "doy")))
    )
}

# only model dNBR above burned threshold
model_results_by_day <- pblapply(1:30, function(day) {
  run_lm_day(y_var = burn_severity_index,
             day = day,
             data = df_subset)
}) %>%
  bind_rows()

# model_results_by_day <- map_dfr(c(1:40), ~ run_lm_day(y_var = burn_severity_index,
#                                                       day = .x,
#                                                       data = df_subset))
results <- model_results_by_day %>% arrange(days_before_fire)

(fig <- ggplot(results) + 
  geom_point(aes(x = days_before_fire, y = r_squared)) +
  scale_x_reverse() +
  labs(x = "Days before fire", y = "R-squared",
       title = sprintf("Rsquared for %s",burn_severity_index)) + 
  theme_cowplot() )

if (SAVE_FIGURES){
  ggsave2(fig,filename = sprintf("figures/%s_top25_R2_fullmod.png",burn_severity_index),
          width = 8,height = 6,bg = "white")
}


# Plot effect sizes over time
results %>% 
  pivot_longer(
    cols = -c(days_before_fire, r_squared, adj_r_squared),
    names_to = "coefficient", values_to = "values") %>% 
  ggplot() + 
  geom_point(aes(x = days_before_fire, y = values)) + 
  scale_x_reverse() +
  labs(x = "Days before fire", y = "Effect size") +
  facet_wrap(~coefficient,scales="free",nrow = 2) + 
  theme_cowplot()

if (SAVE_FIGURES){
  ggsave2(filename = sprintf("figures/%s_effect_sizes_top25.png",
                             burn_severity_index),
          width = 10,height = 8,bg = "white")
}

## b. Run linear model on individual site  ----
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


## c. INLA ----
run_inla <- function(day,y_var, data) {
  ndvi_var <- paste0("NDVI.d_prefire_", day)
  ndmi_var <- paste0("NDMI.d_prefire_", day)
  lst_var <- paste0("cumsum_lst_d_prefire_", day)
  
  predictors <- c(ndvi_var, ndmi_var,lst_var, "elevation", "slope", "northness", "eastness", "burn_doy")
  fixed_effects <- paste(predictors, collapse = " + ")
  # formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))
  formula <- as.formula(paste0(y_var, " ~ ",fixed_effects, 
                               " +  f(fireid, model = 'iid')")) 
  
  model  <- inla(formula, 
                 family = "gaussian",
                 data = data)
  
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

library(ggregplot)
Efxplot(model) + 
  theme_cowplot()
ggsave("figures/inla_10d_efxplot.png",width = 8, height = 8, bg = "white")
# dNBR_model_sel <- INLAModelSel(y_var, predictors, "fireid", "iid", "gaussian", df_subset)
# Finalcovar <- dNBR_model_sel$Removed[[length(dNBR_model_sel$Removed)]]


## d. brms ----
y_var <- "dnbr"
day <- 10
ndvi_var <- paste0("NDVI.d_prefire_", day)
ndmi_var <- paste0("NDMI.d_prefire_", day)
lst_var <- paste0("cumsum_lst_d_prefire_", day)

predictors <- c(ndvi_var, ndmi_var, "elevation", "slope", "northness", "eastness", "doy")
fixed_effects <- paste(predictors, collapse = " + ")
formula <- as.formula(paste(y_var, "~", fixed_effects, "+ (1 | fireid)"))

df_d10_subset <- df_subset %>% 
  group_by(fireid) %>% 
  sample_frac(.05)

# tic()
# bay_qr_day10 <- brm(
#   bf(formula, quantile = 0.75),
#   data = df_d10_subset,
#   family = asym_laplace(),
#   chains = 4,
#   iter = 10000,
#   thin = 10,
#   cores = 4, seed = 1234
# )
# toc()

tic()
model <- brm(
  bf(formula),
  family = "gaussian",
  data = df_d10_subset,
  chains = 4,
  iter = 10000,
  thin = 10,
  cores = 4, 
  seed = 1234)
toc()

summary(model)