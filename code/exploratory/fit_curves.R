library(terra)
library(tidyterra)
library(readr)
library(spatstat.utils)
library(tidyverse)
library(data.table)
library(pbapply)
library(lubridate)
library(tictoc)
set.seed(10)

# 0. Set up functions ----
# ========================.
load_data <- function(fire_attrs,severity_index,frac_int){
  # Fire perimeter attributes
  UTM_TILE_ID <- fire_attrs$opt_UTM_tile
  year <- fire_attrs$tst_year
  
  SCALE_FACTOR <- ifelse(severity_index == "dNBR", 1000, 1)
  
  # Load filtered and formatted dataframe
  fn_filtered_df <- paste0(TABLE_DIR,sprintf(
    "sampled_data/%spct/%s_%s_merged_filtered_%sth_pctile.csv",
    frac_int,FIRE_ID, year,pct_cutoff*100
  ))
  
  if (!file.exists(fn_filtered_df)){
    return(NULL)
  }
  
  df_filtered <- read_csv2(fn_filtered_df, col_names = TRUE,show_col_types = FALSE) %>% 
    as_tibble() %>% 
    mutate(date = as.Date(date),
           burn_date = as.Date(burn_date),
           doy = yday(date))
  
  # Load optimal burn severity raster
  fname_optimal_severity_raster <- optimality_lut %>% 
    filter(fireid == FIRE_ID,
           severity_index == severity_index) %>% 
    pull(fname_severity_raster)
  
  severity_raster <- rast(fname_optimal_severity_raster)  * SCALE_FACTOR
  
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
    DATA_DIR,
    sprintf("feature_layers/%spct/%s_sample_points_%spct_burn_date.gpkg",
            frac_int,FIRE_ID,frac_int)
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

model_fit_smoothedspline <- function(x,y,spar = 0.5) {
  # Using a spline smoother
  smooth.spline(x = x, y = y, spar = spar)
}

# Restructured function to return vector of predictions
model_index_smoothedspline <- function(x,y,full_data, spar = 0.5) {
  
  # Use function to fit model
  model <- model_fit_smoothedspline(x,y,spar)
  
  # Generate predictions for curve plotting (for time-period doy 130-300) 
  pred <- predict(model, data.frame(doy = 130:300))
  pred_df <- data.table(doy = 130:300, spectral_index = pred$y)
  colnames(pred_df)[2] = "spectral_index"
  
  # Get DOY of burn
  doy_burn <- yday(full_data$burn_date[1])
  
  pred_before_burn <- pred_df[doy < doy_burn]
  pred_before_burn[, `:=`(
    # TI_index = revcumsum(spectral_index),
    TI_index = spectral_index,
    days_before_fire = doy_burn - doy
  )]
  
  # Prepare output vector of 40-day predictions
  out_data <- pred_before_burn[days_before_fire <= 40, .(ObservationID = full_data$ObservationID[1], 
                                                         TI_index, 
                                                         days_before_fire)] 
  
  # Reformat as wide datafrmae
  out_data_wide <- dcast(out_data, ObservationID ~ days_before_fire,
                         value.var = "TI_index", 
                         fun.aggregate = list)
  
  return(out_data_wide)
}

#  function to process each ObservationID
process_group <- function(dt, index_name) {
  x <- dt$doy
  y <- dt[[index_name]]
  
  if ("LST" %in% index_name){
    result <- model_LST_polynomial(x, y, dt)
  } else {
    result <- model_index_smoothedspline(x, y, dt, spar = 0.5)
  }
  
  return(result)
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
    ) %>% 
    mutate(ObservationID = full_data$ObservationID[1])
  
  return(out_data)
}

# 1. Configure and load stuff ----
# ================================.
# Config, loading and preparing data
OS <- Sys.info()[['sysname']]
DATA_DIR <- ifelse(OS == "Linux","~/data/","data/")

HLS_DIR <- "~/scratch/raster/hls/"
TABLE_DIR <- paste0(DATA_DIR,"tables/")

severity_index <- "dNBR"
pct_cutoff <- 0.5
OVERWRITE_DATA <- FALSE

frac_to_sample <- 0.01
frac_int <- frac_to_sample *100

# Load fire perimeters
fire_perimeters <- vect(
  paste0(DATA_DIR,"feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")
)

topN_fires <- fire_perimeters %>%
  filter(tst_year >= 2017) %>% 
  arrange(desc(farea)) %>% 
  slice_head(n = 25) 

TEST_ID <- topN_fires$fireid

# Load lookup table
final_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017)

optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

if (length(TEST_ID) > 0){final_lut <- filter(final_lut,fireid %in% TEST_ID)}

for (FIRE_ID in final_lut$fireid){
  cat(sprintf("Processing data from fire: %s \n",FIRE_ID))
  
  filename <- sprintf("model_dataframes/%spct/%s_model_dataframe.csv",
                      frac_int,FIRE_ID)
  
  if (!file.exists(paste0(TABLE_DIR,filename)) || OVERWRITE_DATA){
  
    fire_attrs <- filter(final_lut, fireid == FIRE_ID)
    
    # Load all data
    cat("Loading data...\n")
    data_list <- load_data(fire_attrs, severity_index,frac_int)
    
    if (is.null(data_list)){
      cat("Sampled data table does not exist. Skipping this scar.\n")
      next
    }
    
    df_filtered <- data_list[[1]]
    severity_raster_sample <- data_list[[2]]
    sample_points <- data_list[[3]]
    cropped_severity_raster <- data_list[[4]]
    selected_fire_perimeter <- data_list[[5]]
    
    rm(data_list)
    
    # 2. Fit splines to NDMI & NDVI ----
    # ==================================.
    spar <- 0.5
    
    cat("Fitting smoothing splines for NDMI... \n")
    
    # convert to data table for faster processing
    setDT(df_filtered) 
    
    # Filter observations in pre-fire + 14 d window & points with > 4 observations
    df_filtered_nonan <- df_filtered[
      date < (ymd(burn_date) + 14) & !is.na(DailyMeanNDMI),
      .SD[.N >= 4],
      by = ObservationID
    ]
      
    # Group by ObservationID and apply function in parallel
    df_list <- split(df_filtered_nonan, df_filtered_nonan$ObservationID)
    results <- pblapply(df_list, process_group,cl = 30,  index_name = "DailyMeanNDMI")
    
    # unlist the results sublists to merge into data.table
    results_clean <- lapply(results, function(dt) {
      as.data.table(lapply(dt, function(col) {
        if (is.list(col)) unlist(col, recursive = FALSE) else col
      }))
    })
    
    ndmi_smooth <- rbindlist(results_clean, fill = TRUE) %>% 
      rename_with(~ paste0("NDMI.d_prefire_", .), -ObservationID)
    
    cat("Fitting smoothing splines for NDVI... \n")
    
    # Apply spline to each NDMI time series point (filter out NA and points with < 4 obs.)
    results <- pblapply(df_list, process_group,cl = 30,  index_name = "DailyMeanNDVI")
    # df_filtered_nonan <- df_filtered[!is.na(DailyMeanNDVI)][,
    #                                                         if (.N >= 4) .SD,
    #                                                         by = ObservationID] 
    # df_list <- split(df_filtered_nonan, df_filtered_nonan$ObservationID)

    results_clean <- lapply(results, function(dt) {
      as.data.table(lapply(dt, function(col) {
        if (is.list(col)) unlist(col, recursive = FALSE) else col
      }))
    })
    
    # merge results into  data.table
    ndvi_smooth <- rbindlist(results_clean, fill = TRUE) %>% 
      rename_with(~ paste0("NDVI.d_prefire_", .), -ObservationID)
    
    # Merge dataframes
    data_index <- df_filtered %>% 
      right_join(ndmi_smooth, by = "ObservationID") %>% 
      right_join(ndvi_smooth, by = "ObservationID")
    
    if (any(df_filtered$descals_burn_class != 99,na.rm = T)){
      data_index <- data_index %>% 
        mutate(descals_burn_class = factor(descals_burn_class,
                                           levels = c(0,1),
                                           labels = c("unburned","burned")))
    }
    
    # 3. Fit splines to LST ----
    # ==========================.
    cat("Fitting polynomials to LST... \n")
    
    # Apply spline to each time series point
    df_filtered_nonan <- df_filtered[!is.na(LST) & .N >= 4, ]
    
    # Group by ObservationID and apply function in parallel
    df_list <- split(df_filtered_nonan, df_filtered_nonan$ObservationID)
    results <- pblapply(df_list, process_group,cl = 30,  index_name = "LST")
    
    # unlist the results sublists to merge into data.table
    results_clean <- lapply(results, function(dt) {
      as.data.table(lapply(dt, function(col) {
        if (is.list(col)) unlist(col, recursive = FALSE) else col
      }))
    })
    
    # merge results into  data.table
    lst_smooth <- rbindlist(results_clean, fill = TRUE)
    
    # Merge dataframes
    data_lst <- df_filtered %>% 
      right_join(lst_smooth, by = "ObservationID")
    
    # 4. Build full dataframe for export ----
    # =======================================.
    data_all <- data_lst %>% 
      right_join(data_index,by = c("ObservationID","date")) 
    
    # remove unnecessary or duplicate columns
    data_reduced <- data_all %>% 
      select(-contains(".y")) %>% 
      rename_with(~ sub(".x$", "", .x), everything()) %>%
      distinct()
    
    if (any(df_filtered$descals_burn_class != 99,na.rm = T)){
      data_reduced <- data_reduced %>% 
        mutate(descals_burn_class = factor(descals_burn_class,
                                           levels = c(0,1),
                                           labels = c("unburned","burned")))
    }
    
    # Export as table
    write_csv2(data_reduced,paste0(TABLE_DIR,filename))
    }
}


# 5. Check data fitting ----
# data_reduced <- read_csv2(paste0(TABLE_DIR,filename))
# data_reduced_before <- read_csv2(paste0(TABLE_DIR,"model_dataframes/1pct/14211_model_dataframe_copy.csv"))
# 
# data_reduced_long <- data_reduced %>%
#   select(-contains("lst")) %>%
#   select(-contains("ndvi")) %>%
#   pivot_longer(
#     cols = contains("NDMI.d_prefire_"),
#     names_to = c(".value", "days_before_fire_from_colname"),
#     names_pattern = "(NDMI)\\.d_prefire_([0-9]+)",
#     values_drop_na = TRUE
#   )%>% mutate(
#     days_before_fire = yday(burn_date) - doy,
#     days_before_fire_from_colname = as.integer(days_before_fire_from_colname)) %>% 
#   group_by(ObservationID) %>%
#   filter(days_before_fire_from_colname == days_before_fire) %>%
#   ungroup() 
# 
# data_reduced_before_long <- data_reduced_before %>%
#   select(-contains("lst")) %>%
#   select(-contains("ndvi")) %>%
#   pivot_longer(
#     cols = contains("NDMI.d_prefire_"),
#     names_to = c(".value", "days_before_fire_from_colname"),
#     names_pattern = "(NDMI)\\.d_prefire_([0-9]+)",
#     values_drop_na = TRUE
#   ) %>% 
#   mutate(
#     days_before_fire = yday(burn_date) - doy,
#     days_before_fire_from_colname = as.integer(days_before_fire_from_colname)) %>% 
#   group_by(ObservationID) %>%
#   filter(days_before_fire_from_colname == days_before_fire) %>%
#   ungroup()
# 
# 
# ggplot(data_reduced_long) + 
#   geom_point(aes(x = days_before_fire,y = DailyMeanNDMI,
#                 group = ObservationID),
#             color = "red",alpha = .5) +
#   geom_line(aes(x = days_before_fire,y = NDMI,
#                 group = ObservationID),
#             color = "blue",alpha = .5) +
#   geom_vline(aes(xintercept = 0)) + 
#   scale_x_reverse() +
#   theme_cowplot()
# 
# ggsave(sprintf("figures/NDMI_splines_perpoint_shorter_period_%s.png",FIRE_ID),
#        bg = "white",height = 12,width = 12)


# name <- "dgemi"
# q <- .5
# n_points <- length(unique(data_reduced %>% 
#                             filter(quantile(!!sym(name), q,na.rm=T)<!!sym(name)) %>%
#                             pull(ObservationID)))
# data_reduced %>% distinct(!!sym(name),.keep_all = T) %>% 
#   ggplot(aes(x = !!sym(name))) + 
#   stat_ecdf(color = "salmon") + 
#   geom_hline(aes(yintercept = q)) +
#   geom_vline(aes(xintercept = quantile(!!sym(name), q,na.rm=T))) +
#   theme_cowplot() +
#   labs(x = name,
#        title = sprintf("Number of points with \n %s >= %s decile: %s",name, q,n_points)  )
# 
# ggsave(sprintf("figures/%s_ecdf_%s_cutoff_%s.png",name,q,FIRE_ID),
#        bg = "white",height = 12,width = 12)
