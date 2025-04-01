library(terra)
library(tidyterra)
library(spatstat.utils)
library(tidyverse)
library(cowplot)
library(patchwork)
library(lubridate)
library(quantreg)
set.seed(10)

# 0. Set up functions ----
# ========================.
load_data <- function(final_lut,severity_index){
  # Fire perimeter attributes
  UTM_TILE_ID <- final_lut$opt_UTM_tile
  year <- final_lut$tst_year
  
  SCALE_FACTOR <- ifelse(severity_index == "dNBR", 1000, 1)
  
  # Load dNBR raster
  severity_raster_list <- list.files(path = paste0(HLS_DIR,"severity_rasters"),
                                     pattern = sprintf("^%s_%s_%s.*\\.tif$",
                                                       severity_index,UTM_TILE_ID,year),
                                     full.names = TRUE)
  severity_raster <- rast(severity_raster_list[1]) * SCALE_FACTOR
  
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
  fname_sample_points <- sprintf("~/data/feature_layers/%s_sample_points_burn_date.gpkg",
                                 FIRE_ID)
  sample_points <- vect(fname_sample_points) %>% 
    project(crs(severity_raster)) %>% 
    mutate(ObservationID = 1:nrow(.))
  
  # Extract data at random points 
  severity_raster_sample <- terra::extract(severity_raster, sample_points,
                                           ID = FALSE, xy = TRUE)
  
  # Load filtered and formatted dataframe
  fn_filtered_df <- sprintf(
    "~/data/tables/sampled_data/%s_%s_merged_filtered_%sth_pctile.csv",
    FIRE_ID, year,pct_cutoff*100
  )
  
  df_filtered <- read.csv2(fn_filtered_df) %>% 
    mutate(date = as.Date(date),
           burn_date = as.Date(burn_date),
           doy = yday(date))
  
  return(list(df_filtered,
              severity_raster_sample,
              sample_points,
              cropped_severity_raster,
              selected_fire_perimeter))
}

model_fit_smoothedspline <- function(x,y,spar = 0.5) {
  # Using a spline smoother
  smooth.spline(x = x, y = y, spar = spar)
}

# Restructured function to return vector of predictions
model_index_smoothedspline <- function(x,y,full_data, spar = 0.5) {
  
  # Use function to fit model
  model <- model_fit_smoothedspline(x,y,spar)
  
  # Generate predictions for curve plotting (for time-period doy 130-300) 
  pred <- data.frame(predict(model, data.frame(doy = 130:300))) %>% 
    rename(doy = doy, index = doy.1)
  
  # Get DOY of burn
  doy_burn <- yday(full_data$burn_date[1])
  
  pred_before_burn <- pred %>% 
    filter(doy < doy_burn) %>% 
    mutate(TI_index = revcumsum(index),
           days_before_fire = doy_burn - doy)
  
  # Prepare output vector of 40-day predictions
  out_data <- pred_before_burn %>% 
    filter(days_before_fire <= 40) %>% 
    select(days_before_fire, TI_index) %>% 
    pivot_wider(
      names_from = days_before_fire,
      values_from = TI_index,
      names_prefix = "d_prefire_"
    ) 
  
  return(out_data)
}

# Restructured function to return vector of predictions
model_LST_polynomial <- function(x,y,full_data) {
  
  # Use function to fit model
  model <- lm(y ~ poly(x, 2, raw=TRUE))
  
  # Generate predictions for curve plotting (for time-period doy 130-300) 
  pred <- data.frame(preds = predict(model, data.frame(x = 130:300))) %>% 
    mutate(doy = 130:300)
  
  # Get DOY of burn
  doy_burn <- yday(full_data$burn_date[1])
  
  pred_before_burn <- pred %>% 
    filter(doy < doy_burn) %>% 
    mutate(cumsum_lst = revcumsum(preds),
           lst_pred = preds,
           days_before_fire = doy_burn - doy)
  
  # Prepare output vector of 40-day predictions
  out_data <- pred_before_burn %>% 
    filter(days_before_fire <= 40) %>% 
    select(days_before_fire, cumsum_lst,lst_pred) %>% 
    pivot_wider(
      names_from = days_before_fire,
      values_from = c(cumsum_lst,lst_pred),
      names_prefix = "d_prefire_"
    ) 
  
  return(out_data)
}

# 1. Configure and load stuff ----
# ================================.
# Config, loading and preparing data
HLS_DIR <- "~/data/raster/hls/"
TABLE_DIR <- "~/data/tables/"

FIRE_IDs <- c(14664,10792,17548,14211)
index_name <- "NDMI"
severity_index <- "dNBR"
pct_cutoff <- 0.5
SAVE_FIGURES <- TRUE
OVERWRITE_DATA <- TRUE

# Load lookup table
final_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017)

# Load fire perimeters
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)

for (FIRE_ID in FIRE_IDs){
  final_lut <- filter(final_lut, fireid == FIRE_ID)
  
  # Load all data
  data_list <- load_data(final_lut, severity_index)
  
  df_filtered <- data_list[[1]]
  severity_raster_sample <- data_list[[2]]
  sample_points <- data_list[[3]]
  cropped_severity_raster <- data_list[[4]]
  selected_fire_perimeter <- data_list[[5]]
  
  rm(data_list)
  
  # 2. Fit splines to NDMI & NDVI ----
  # ==================================.
  spar <- 0.5
  
  # Apply spline to each time series point
  ndmi_smooth <- df_filtered %>%
    filter(!is.na(DailyMeanNDMI)) %>% 
    group_by(ObservationID) %>% 
    group_modify(
      ~model_index_smoothedspline(.x$doy,.x$DailyMeanNDMI,.x,spar = spar)
      ) %>% 
    rename_with(~ paste0("NDMI.", .), -ObservationID)
  
  ndvi_smooth <- df_filtered %>%
    filter(!is.na(DailyMeanNDVI)) %>% 
    group_by(ObservationID) %>% 
    group_modify(
      ~model_index_smoothedspline(.x$doy,.x$DailyMeanNDVI,.x,spar = spar)
      ) %>% 
    rename_with(~ paste0("NDVI.", .), -ObservationID)
  
  # Merge dataframes
  data_index <- df_filtered %>% 
    right_join(ndmi_smooth, by = "ObservationID") %>% 
    right_join(ndvi_smooth, by = "ObservationID")
  
  if (any(df_filtered$descals_burn_class != 99,na.rm = T)){
    data_index <- data_index %>% 
      mutate(descals_burn_class = factor(descals_burn_class,
                                         labels = c("unburned","burned")))
  }
  
  # 3. Fit splines to LST ----
  # ==========================.
  
  # Apply spline to each time series point
  lst_smooth <- df_filtered %>%
    filter(!is.na(doy)) %>% 
    group_by(ObservationID) %>% 
    group_modify(~model_LST_polynomial(.x$doy,.x$LST,.x))
  
  # Merge dataframes
  data_lst <- df_filtered %>% 
    right_join(lst_smooth, by = "ObservationID")
  
  if (any(df_filtered$descals_burn_class != 99,na.rm = T)){
    data_lst <- data_lst %>% 
      mutate(descals_burn_class = factor(descals_burn_class,
                                       labels = c("unburned","burned")))
  }
  
  # 4. Build full dataframe for export ----
  # =======================================.
  data_all <- data_lst %>% 
    right_join(data_index,by = c("ObservationID","date")) 
  
  if (any(df_filtered$descals_burn_class != 99,na.rm = T)){
    data_all <- data_all %>% 
      mutate(descals_burn_class = factor(descals_burn_class,
                                         labels = c("unburned","burned")))
  }
  
  # remove unnecessary or duplicate columns
  data_reduced <- data_all %>% 
    select(-contains(".y")) %>% 
    rename_with(~ sub(".x$", "", .x), everything()) %>%
    distinct()
  
  # Export as table
  filename <- sprintf("model_dataframes/%s_%s_model_dataframe.csv",
                      FIRE_ID,severity_index)
  
  if (!file.exists(paste0(TABLE_DIR,filename)) || OVERWRITE_DATA){
    write.csv2(data_reduced,paste0(TABLE_DIR,filename))
  }
}