library(terra)
library(sf)
library(tidyterra)
library(rts)
library(readr)
# library(gt)
library(tidyverse)
library(cowplot)
library(patchwork)
library(lubridate)
print(getwd())

set.seed(10)

# 0. Configure functions ----
# ===========================.

sample_dnbr_points <- function(raster, sample_pct = 0.10) {
  # Frequency table for raster bins
  freq_table <- freq(raster, digits = 0)
  
  # Calculate the number of points per bin
  sample_size_per_bin <- ceiling(sample_pct * freq_table[, "count"])
  
  # Initialize empty SpatVector for sampled points
  sample_points <- vect()
  
  # Efficient sampling using vectorized approach
  for (i in seq_len(nrow(freq_table))) {
    category_value <- freq_table[i, "value"]
    
    # Skip if no points to sample
    if (sample_size_per_bin[i] == 0) next
    
    # Mask for the current category
    category_mask <- raster == category_value
    
    # Sample random points
    points <- spatSample(category_mask, size = sample_size_per_bin[i],
                         method = "random", na.rm = TRUE, as.points = TRUE)
    
    # Append sampled points if valid
    if (!is.null(points)) {
      sample_points <- rbind(sample_points, points)
    }
  }
  
  return(sample_points)
}

assign_rast_time <- function(lyr) {
  # Extract the layer's timestamp
  timestamp <- sub(
    ".*\\.(\\d{7}T\\d{6}).*\\.tif$", "\\1",
    basename(sources(lyr))
  )
  
  # Convert string to datetime
  datetime <- as.POSIXct(timestamp, format = "%Y%jT%H%M%S", tz = "UTC")
  
  # Assign time to raster layer
  terra::time(lyr, "seconds") <- datetime
  # names(lyr) <- datetime
  
  return(lyr)
}

stack_time_series <- function(hls_dir,utm_tile_id, index_name){
  # List spectral_index COGs
  spectral_index_files_s30 <- list.files(hls_dir,
                                         pattern = sprintf("HLS.S30.T%s.*%s\\.tif$",
                                                           utm_tile_id,index_name))
  spectral_index_files_l30 <- list.files(hls_dir,
                                         pattern = sprintf("HLS.L30.T%s.*%s\\.tif$",
                                                           utm_tile_id,index_name))
  
  cat("Loading rasters...\n")
  
  # Load rasters to list
  spectral_index_s30_list <- lapply(paste0(hls_dir, spectral_index_files_s30),rast)
  spectral_index_l30_list <- lapply(paste0(hls_dir, spectral_index_files_l30),rast)
  
  cat("Creating image stacks...\n")
  
  # Concatenate as single raster and assign time dimension
  spectral_index_s30 <- rast(spectral_index_s30_list)
  spectral_index_s30 <- sapp(spectral_index_s30, assign_rast_time)
  
  spectral_index_l30 <- rast(spectral_index_l30_list)
  spectral_index_l30 <- sapp(spectral_index_l30, assign_rast_time)
  
  # Merge S30 and L30 stack
  spectral_index <- c(spectral_index_s30, spectral_index_l30)
  
  # Shift NDMI by 1
  if (index_name == "NDMI"){
    spectral_index <- spectral_index + 1
  }
  
  cat("done.\n")
  return(spectral_index)
}

get_ls_datetime <- function(raster){
  datenum <- global(raster,"max",na.rm=T)
  date <- as.Date(as.numeric(datenum), origin = "1970-01-01")
  return(date)
}

# 1. Loading and preparing data: ----
# ===================================.

# Set name of severity index to extract data for
severity_index <- "dNBR"

# TRUE to overwrite existing data form time series extraction
OVERWRITE_DATA <- FALSE

TEST_ID <- c(14664,17548,10792) # fire ID for part of the large fire scar

# Define percentile for sample cutoff
pct_cutoff <- 0.5

OS <- Sys.info()[['sysname']]

# Output directory for sample tables
TABLE_DIR <- ifelse(OS == "Linux","~/data/tables/","data/tables/")
OUT_DIR <- paste0(TABLE_DIR,"sampled_data/")

# Load lookup tables
final_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017) %>% 
  filter(fireid %in% TEST_ID)

dem_lut <- read.csv(paste0(TABLE_DIR,"dem_fire_perim_intersect.csv")) # DEM tiles

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)

# Run data_extraction_part1.R first

# 2. Run burn date assignment in python ----
# ==========================================.
# Run on python instance on S3IT

# ...or use below if on local machine

# give permission: chmod +x code/data_processing/run_burn_date_extraction.sh
# system("code/data_processing/run_burn_date_extraction.sh")

# 3. Stack rasters and extract data----
# =====================================.
for(i in 1:nrow(final_lut)) {
  row <- final_lut[i, ]
  
  # extract fire perimeter attributes
  UTM_TILE_ID <- row$opt_UTM_tile
  year <- row$tst_year
  FIRE_ID <- row$fireid
  
  cat(sprintf("Extracting data for fire %s in UTM tile: %s \n",
              FIRE_ID,UTM_TILE_ID))
  
  # Load burn severity rasters
  severity_rasters <- list.files(path = "~/data/raster/hls/severity_rasters",
                                 pattern = sprintf("^%s_%s_%s.*\\.tif$",
                                                   severity_index,UTM_TILE_ID,year),
                                 full.names = TRUE)
  rast_burn_severity <- rast(severity_rasters[1])
  
  # Load sample points with burn dates
  fname_sample_points <- sprintf("~/data/feature_layers/%s_sample_points_burn_date.gpkg",
                                 FIRE_ID)
  sample_points <- vect(fname_sample_points) %>% 
    project(crs(rast_burn_severity)) %>% 
    mutate(ObservationID = 1:nrow(.))
  
  # Extract dNBR at random points 
  dnbr_sample <- terra::extract(rast_burn_severity, sample_points,
                                ID = FALSE, xy = TRUE)
  
  ## a. HLS spectral indices ----
  
  for (index_name in c("NDMI","NDVI")){
    # Check if filename exists and reprocessing is activated
    filename <- sprintf("%s_sampled_%s_%s.csv",index_name,FIRE_ID,year)
    
    fmt  <- ifelse(index_name == "NDMI","%Y-%m-%d %H:%M:%S","X%Y.%m.%d.%H.%M.%S")
    
    # only stack rasters if needed (takes a while)
    if (!file.exists(paste0(OUT_DIR,filename)) || OVERWRITE_DATA){
      cat("Stacking and extracting spectral index data... \n")
      
      HLS_DIR <- "~/scratch/raster/hls/processed/"
      image_stack <- stack_time_series(HLS_DIR,UTM_TILE_ID, index_name)
      
      # create raster time series object
      rt <- rts(image_stack, time(image_stack))
      
      # Extract spectral index time series at each sample point
      df_spectral_index <- terra::extract(rt, sample_points) %>%
        t() %>%
        as_tibble() 
      
      df_long <- df_spectral_index %>% 
        mutate(ObservationID = 1:nrow(.),
               burn_severity = dnbr_sample[,1],
               burn_date = sample_points$burn_date,
               descals_burn_class = sample_points$descals_burned,
               .before = 1) %>% 
        gather(key = "Time", value = !!sym(index_name), 
               -c(ObservationID,burn_severity,burn_date,descals_burn_class)) %>% 
        mutate(Time = ymd_hms(Time))
      
      # Export as table
      cat("Writing data to csv...\n")
      write_csv2(df_long,paste0(OUT_DIR,filename))
      
    } 
  }
  
  # Load  NDMI time series at points 
  cat("Loading NDMI samples...\n")
  filename <- sprintf("NDMI_sampled_%s_%s.csv",FIRE_ID,year)
  
  df_ndmi <- read_csv2(paste0(OUT_DIR,filename), col_names = TRUE) %>% 
    as_tibble()
  df_ndmi <- df_ndmi %>% 
    mutate(Time = ymd_hms(Time))
  
  # Load NDVI time series at points
  cat("Loading NDVI samples...\n")
  filename <- sprintf("NDVI_sampled_%s_%s.csv",FIRE_ID,year)
  
  df_ndvi <- read_csv2(paste0(OUT_DIR,filename), col_names = TRUE) %>% 
    as_tibble()
  df_ndvi <- df_ndvi %>% 
    mutate(Time = ymd_hms(Time))
  
  ## b. Landsat-8 LST ----
  
  # Check if filename exists and reprocessing is activated
  filename <- sprintf("LST_sampled_%s_%s.csv",FIRE_ID,year)
  
  # only stack rasters if needed (takes a while)
  if (!file.exists(paste0(OUT_DIR,filename)) || OVERWRITE_DATA){
    cat("Stacking and extracting LST data... \n")
    
    # List LST tiffs (except "tile" files)
    LS8_DIR <- "~/data/raster/landsat8"
    lst_files <- list.files(LS8_DIR, 
                            pattern = paste0("^", FIRE_ID, ".*\\.tif$"),
                            full.names = TRUE)
    
    lst_files <- lst_files[!grepl("tile", lst_files)]
    
    # Load LST rasters to list
    lst_list <- lapply(lst_files, function(fn){
      r <- rast(fn)[[1]]
      names(r) <- "TIR"
      return(r)})
    lst_time_list <- lapply(lst_files, function(fn){
      r <- rast(fn)[[2]]
      names(r) <- "TIME"
      return(r)})
    
    # Concatenate as single raster 
    lst <- rast(lst_list) %>% project(crs(rast_burn_severity))
    
    # Extract dates and format to date
    lst_dates <- lapply(lst_time_list,get_ls_datetime)
    lst_dates_df <- data.frame(datenum = unlist(lst_dates)) %>% 
      mutate(date = as.Date(datenum, origin = "1970-01-01"))
    
    # assign time dimension to raster stack
    terra::time(lst, "days") <- lst_dates_df$date
    
    # Remove duplicate LST rasters
    n <- time(lst)
    dup.idx <- which(duplicated(n))
    
    if (length(dup.idx) > 0) {
      lst <- lst[[-dup.idx]]
    }
    
    # Convert raster stack to time series stack
    rt_lst <- rts(lst, time(lst))
    
    # Extract LST time series at each sample point
    fmt = "%Y-%m-%d"
    
    df_lst <- terra::extract(rt_lst, sample_points) %>%
      t() %>%
      as_tibble() %>% 
      mutate(ObservationID = 1:nrow(.),
             burn_severity = dnbr_sample[,1],
             burn_date = sample_points$burn_date,
             descals_burn_class = sample_points$descals_burned,
             .before = 1) %>% 
      gather(key = "Time", value = "LST", 
             -c(ObservationID,burn_severity,burn_date,descals_burn_class)) %>% 
      mutate(Time = as_datetime(Time, format = fmt))
    
    # Export as table
    cat("Writing data to csv...\n")
    write_csv2(df_lst,paste0(OUT_DIR,filename))
  } else {
    
    cat("Loading LST samples...\n")
    df_lst <- read_csv2(paste0(OUT_DIR,filename), col_names = TRUE) %>% 
      as_tibble()
    
    # df_lst <- df_lst %>%
    #   mutate(Time = as_datetime(Time, format = fmt))
  }
  
  # 5. Filter and format sampled data ----
  
  # Set name of filtered and formatted dataframe
  fn_filtered_df <- sprintf(
    "%s_%s_merged_filtered_%sth_pctile.csv",
    FIRE_ID, year,pct_cutoff*100
    )
  
  # Check if data frame exists or needs to be overwritten
  if (!file.exists(paste0(OUT_DIR,fn_filtered_df)) || OVERWRITE_DATA){
    cat("Filtering and formatting raster data... \n")
    
    # Concatenate the data frames
    df_hls <- df_ndmi %>%
      left_join(df_ndvi %>% select(ObservationID, Time, NDVI), 
                by = c("ObservationID", "Time")) %>% 
      mutate(date = as.Date(Time))
    
    # Compute daily means
    df_daily_spectral_index <- df_hls %>% 
      group_by(date, ObservationID) %>% 
      summarise(DailyMeanNDMI = mean(NDMI, na.rm = TRUE),
                DailyMeanNDVI = mean(NDVI, na.rm = TRUE),
                burn_severity = first(burn_severity),
                burn_date = first(burn_date),
                descals_burn_class = first(descals_burn_class))
    
    # get nr. of valid observations
    valid_counts_by_id <- df_daily_spectral_index %>%
      group_by(ObservationID) %>%
      summarise(valid_count = sum(!is.na(DailyMeanNDMI)))
    
    thr_nobs <- quantile(valid_counts_by_id$valid_count,pct_cutoff , na.rm = TRUE)
    
    # Filter out observations with fewer than X observations
    df_filtered <- df_daily_spectral_index %>%
      inner_join(valid_counts_by_id, by = "ObservationID") %>%
      filter(valid_count >= thr_nobs,
             year(date) == year) %>% 
      select(-valid_count) %>% 
      # Flag pre- & post-fire observations
      mutate(
        burn_date = ymd_hms(burn_date),
        BeforeBurnDate = date < burn_date
      )
    
    # Extract topography at random points
    fname_dem <- sprintf("~/data/raster/arcticDEM/%s_dem_30m.tif",FIRE_ID)
    dem_30m <- rast(fname_dem)
    names(dem_30m) <- "elevation"
    
    elevation_sample <- terra::extract(dem_30m, sample_points, bind = TRUE) %>% 
      as.data.frame() %>% select(c(ObservationID,elevation))
    
    # Slope
    slope <- terrain(dem_30m,'slope')
    slope_sample <- terra::extract(slope,sample_points, bind = TRUE) %>% 
      as.data.frame() %>% select(c(ObservationID,slope))
    
    # Aspect
    aspect <- terrain(dem_30m,'aspect')
    northness <- cos(aspect * pi/180)
    names(northness) <- 'northness'
    
    eastness <- sin(aspect * pi/180)
    names(eastness) <- 'eastness'
    
    northness_sample <- terra::extract(northness, sample_points, bind = TRUE) %>% 
      as.data.frame() %>% select(c(ObservationID,northness))
    eastness_sample <- terra::extract(eastness, sample_points, bind = TRUE) %>% 
      as.data.frame() %>% select(c(ObservationID,eastness))
    
    # combine all DEM data
    df_list <- list(elevation_sample, slope_sample, 
                    northness_sample,eastness_sample)
    
    # merge all data frames in list
    dem_sample <- df_list %>% reduce(full_join, by='ObservationID')
    
    df_filtered <- full_join(df_filtered, dem_sample,by = "ObservationID")
    
    # Merge with LST data
    df_filtered <- df_filtered %>% 
      left_join(df_lst %>% select(ObservationID, Time, LST), 
                by = c("ObservationID" = "ObservationID", 
                       "date" = "Time"))
    
    # Write filtered data frame to CSV
    write_csv2(df_filtered,paste0(OUT_DIR,fn_filtered_df))
  }
}
