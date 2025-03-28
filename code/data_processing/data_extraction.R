library(terra)
library(sf)
library(tidyterra)
library(rts)
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
  
  cat("Loading rasters...")
  
  # Load rasters to list
  spectral_index_s30_list <- lapply(paste0(hls_dir, spectral_index_files_s30),rast)
  spectral_index_l30_list <- lapply(paste0(hls_dir, spectral_index_files_l30),rast)
  
  cat("Creating image stacks...")
  
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
  
  cat("done.")
  return(spectral_index)
}

get_ls_datetime <- function(raster){
  datenum <- global(raster,"max",na.rm=T)
  date <- as.Date(as.numeric(datenum), origin = "1970-01-01")
  return(date)
}

read_hls_data_frames <- function(index_name, UTM_TILE_ID,year,
                                 sample_points, dnbr_sample){
  filename <- sprintf("%s_sampled_%s_%s.csv",index_name,UTM_TILE_ID,year)
  
  df <- read.csv2(paste0(TABLE_DIR,filename)) %>% 
    select(-1) %>% 
    mutate(burn_severity = dnbr_sample$burn_severity) %>% 
    as_tibble() %>% 
    mutate(burn_date = sample_points$burn_date,
           descals_burn_class = sample_points$descals_burned)
  
  if (index_name == "LST"){
    fmt = "X%Y.%m.%d"
  } else{
    fmt = "X%Y.%m.%d.%H.%M.%S"
  }
  
  df_long <- df %>%
    mutate(ObservationID = 1:nrow(df)) %>%  # Add a column for observation IDs (row numbers)
    gather(key = "Time", value = index_name, 
           -c(ObservationID,burn_severity,burn_date,descals_burn_class)) %>% 
    mutate(Time = as.POSIXct(Time, format = fmt)) 
  
  return(df_long)
}


# 1. Loading and preparing data: ----
# ===================================.

# Set name of severity index to extract data for
severity_index <- "dNBR"

# Set name of spectral index to extract data for
index_name <- "NDMI"

# TRUE to overwrite existing data form time series extraction
OVERWRITE_DATA <- TRUE

TEST_ID <- c(17548,10792,14664) # fire ID for part of the large fire scar

# Define percentile for sample cutoff
pct_cutoff <- 0.5

OS <- Sys.info()[['sysname']]

# Output directory for sample tables
TABLE_DIR <- ifelse(OS == "Linux", 
                    "~/data/tables/","data/tables/")

# Load lookup tables
final_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017) %>% 
  filter(fireid %in% TEST_ID)

dem_lut <- read.csv(paste0(TABLE_DIR,"dem_fire_perim_intersect.csv")) # DEM tiles

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)

# 2. Execute spatial sampling ----
#=================================.

# Apply data extraction to all fire perimeters
for(i in 1:nrow(final_lut)) {
  row <- final_lut[i, ]
  
  # extract fire peimeter attributes
  UTM_TILE_ID <- row$opt_UTM_tile
  year <- row$tst_year
  FIRE_ID <- row$fireid
  
  print(sprintf("Extracting data for fire %s in UTM tile: %s",FIRE_ID,UTM_TILE_ID))
  
  # Load burn severity rasters
  severity_rasters <- list.files(path = "~/data/raster/hls/severity_rasters",
                                 pattern = sprintf("^%s_%s_%s.*\\.tif$",
                                                   severity_index,UTM_TILE_ID,year),
                                 full.names = TRUE
                                 )
  rast_burn_severity <- rast(severity_rasters[1])
  
  # get fire perimter
  selected_fire_perimeter <- fire_perimeters %>% 
    filter(fireid  == FIRE_ID) %>%
    project(crs(rast_burn_severity))
  
  # Check if fire extent is above Arctic circle
  min_lat <- ext(project(selected_fire_perimeter, "EPSG:4326"))[3]
  is_above_arctic <- min_lat > 66.56
  
  # Prepare Descals et al. (2022) burned area maps if available
  if (is_above_arctic){
    
    # Load descals burned area raster for the selected fire
    descals_tilename <- final_lut %>% 
      filter(fireid  == FIRE_ID) %>% 
      select(descals_file)
    
    # read raster
    descals_rast <- rast(paste0("~/data/raster/burned_area_descals/",descals_tilename)) 
    
    # Reclassify to binary burned area raster where 1 is burned that year (of selected fire perimeter), 0 not
    year_value_descals <- year - 1990 # value 30 is for burned area in 2020
    
    # ... and crop with reprojected fire perimeter (and 1km buffer)
    descals_rast_bin <- terra::as.int(descals_rast == year_value_descals) %>% 
      crop(buffer(
        project(selected_fire_perimeter,crs(descals_rast))
        ,1e3))
    
    # Reproject raster if necessary
    if (crs(descals_rast_bin) != crs(selected_fire_perimeter)){
      descals_rast_bin <- project(descals_rast_bin,
                                  crs(selected_fire_perimeter),
                                  method = "near")
    }
  }
  
  # Load DEM tiles for this fire perimeter
  dem_tiles <- dem_lut %>% 
    filter(fireid == FIRE_ID) %>% 
    mutate(filename = paste0(dem_id,"_dem.tif")) %>% 
    pull(filename)
  
  dem_list <- lapply(paste0("~/scratch/raster/arcticDEM/",dem_tiles), rast,
                     drivers="GTiff")
  
  selected_fire_perimeter_stereo <- project(selected_fire_perimeter,
                                            crs(dem_list[[1]])) %>% 
    buffer(1000)
  
  cropped_dems <- lapply(dem_list,
                         function(x) crop(x,ext(selected_fire_perimeter_stereo)))
  
  # Mosaic DEMS
  rsrc <- sprc(cropped_dems)
  dem_mos <- mosaic(rsrc)
  
  # Resample to 30m HLS UTM
  raster_grid_template <- rast(dem_mos)
  res(raster_grid_template) <- 30
  
  # Resample to 30m resolution using bilinear interpolation
  dem_30m <- resample(dem_mos, raster_grid_template, method="cubicspline") 
  names(dem_30m) <- 'elevation'
  
  # delete 2m DEM mosaic to free up space
  rm(cropped_dems,dem_mos)
  
  # reproject
  dem_30m_utm <- project(dem_30m,crs(selected_fire_perimeter))
  
  # export resampled DEM
  fname_dem_out <- paste0("~/data/raster/arcticDEM/",
                          FIRE_ID,
                          "_dem_30m.tif")
  writeRaster(dem_30m_utm,filename = fname_dem_out,overwrite = T)
  
  rm(dem_30m,dem_30m_utm)
  
  # Buffer the selected fire perimeter
  fire_perimeter_buffered <- buffer(selected_fire_perimeter, 1200)
  
  dnbr_in_perimeter <- rast_burn_severity %>% 
    mask(fire_perimeter_buffered, updatevalue = NA) %>% 
    crop(fire_perimeter_buffered) * 1000
  
  fname_sample_points <- sprintf("~/data/feature_layers/%s_sample_points.gpkg",
                                 FIRE_ID)
  
  # Sample points only when gpkg file doens't exist
  if (!file.exists(fname_sample_points)){
    print("Creating random point sample... \n")
    
    # Set dNBR bin size for sampling (20 if dNBR * 1000)
    binwidth <- 20
    
    # get bins for the raster
    bins <- seq(minmax(dnbr_in_perimeter)[1],
                minmax(dnbr_in_perimeter)[2],
                binwidth)
    
    # Create binned raster for sampling
    rast_binned <- classify(dnbr_in_perimeter,
                            bins, 
                            include.lowest=TRUE,
                            brackets=TRUE)
    
    # Set proportion of sampled values per bin
    frac_to_sample <- 0.01
    
    # Create spatial point sample
    sample_points <- sample_dnbr_points(rast_binned, sample_pct = frac_to_sample)
    
    if (is_above_arctic){
      # Extract Descals et al. (2022) burn class to point
      sample_points$descals_burned <- terra::extract(descals_rast_bin, 
                                                     sample_points, ID = F)
    } else {
      sample_points$descals_burned <- 99
    }
    
    # Convert to sf
    sample_points_sf <- st_as_sf(sample_points)
    
    # Export points as gpkg
    st_write(sample_points_sf, fname_sample_points,
             layer = "sample_points", delete_layer = TRUE)
  }
  
  print("Random point GPKG already exists, continue. \n")
  
}

# 3. Run burn date assignment in python ----
# ==========================================.
# Run on python instance on S3IT

# ...or use below if on local machine

# give permission: chmod +x code/data_processing/run_burn_date_extraction.sh
# system("code/data_processing/run_burn_date_extraction.sh")

# 4. Stack rasters and extract data----
# =====================================.
for(i in 1:nrow(final_lut)) {
  row <- final_lut[i, ]
  
  # extract fire perimeter attributes
  UTM_TILE_ID <- row$opt_UTM_tile
  year <- row$tst_year
  FIRE_ID <- row$fireid
  
  print(sprintf("Extracting data for fire %s in UTM tile: %s",FIRE_ID,UTM_TILE_ID))
  
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
  
  ## a. HLS spectral indices ----
  
  # Check if filename exists and reprocessing is activated
  filename <- sprintf("%s_sampled_%s_%s.csv",index_name,FIRE_ID,year)
  
  # only stack rasters if needed (takes a while)
  if (!file.exists(paste0(TABLE_DIR,filename)) || OVERWRITE_DATA){
    print("Stacking and extracting spectral index data... \n")
    
    HLS_DIR <- "~/scratch/raster/hls/processed/"
    image_stack <- stack_time_series(HLS_DIR,UTM_TILE_ID, index_name)
    
    # create raster time series object
    rt <- rts(image_stack, time(image_stack))
    
    # Extract spectral index time series at each sample point
    df_spectral_index <- terra::extract(rt, sample_points) %>%
      t() %>%
      as_tibble()
    
    # Export as table
    write.csv2(df_spectral_index,paste0(TABLE_DIR,filename))
  }
  
  ## b. Landsat-8 LST ----
  
  # Check if filename exists and reprocessing is activated
  filename <- sprintf("LST_sampled_%s_%s.csv",UTM_TILE_ID,year)
  
  # only stack rasters if needed (takes a while)
  if (!file.exists(paste0(TABLE_DIR,filename)) || OVERWRITE_DATA){
    print("Stacking and extracting LST data... \n")
    
    # List LST tiffs
    LS8_DIR <- "~/scratch/raster/landsat8/"
    lst_files <- list.files(LS8_DIR,pattern = paste0("^", FIRE_ID, ".*\\.tif$"))
    
    # Load LST rasters to list
    lst_list <- lapply(paste0(LS8_DIR, lst_files), function(fn) rast(fn)$TIR)
    lst_time_list <- lapply(paste0(LS8_DIR, lst_files), function(fn) rast(fn)$TIME)
    
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
    lst <- lst[[-dup.idx]]
    
    # Convert raster stack to time series stack
    rt_lst <- rts(lst, time(lst))
    
    # Extract LST time series at each sample point
    df_lst <- terra::extract(rt_lst, sample_points) %>%
      t() %>%
      as_tibble()
    
    # Export as table
    write.csv2(df_lst,paste0(TABLE_DIR,filename))
  }
  
  # 5. Filter and format sampled data ----
  
  # Set name of filtered and formatted dataframe
  fn_filtered_df <- sprintf("data/tables/%s_filtered_%sth_pctile.csv",
                            index_name,pct_cutoff*100)
  
  # Check if data frame exists or needs to be overwritten
  if (!file.exists(fn_filtered_df) || OVERWRITE_DATA){
    print("Filtering and formatting raster data... \n")
    
    # Extract dNBR at random points 
    dnbr_sample <- terra::extract(rast_burn_severity, sample_points,
                                  ID = FALSE, xy = TRUE)
    
    # Load NDVI and NDMI time series at points 
    df_ndmi <- read_hls_data_frames("NDMI", UTM_TILE_ID,year,sample_points, dnbr_sample)
    df_ndvi <- read_hls_data_frames("NDVI",UTM_TILE_ID,year,sample_points, dnbr_sample)
    df_lst <- read_hls_data_frames("LST",UTM_TILE_ID,year,sample_points, dnbr_sample)
    
    df_lst <- df_lst %>%
      mutate(Time = as.Date(Time))
    
    # Concatenate the data frames
    df_hls <- df_ndmi %>% 
      mutate(NDVI = df_ndvi$index_name) %>% 
      rename(NDMI = index_name) %>% 
      mutate(date = as.Date(Time))
    
    # Compute daily means
    df_daily_spectral_index <- df_hls %>% 
      group_by(date, ObservationID) %>% 
      summarise(DailyMeanNDMI = mean(NDMI, na.rm = TRUE),
                DailyMeanNDVI = mean(NDVI, na.rm = TRUE),
                burn_severity = first(rast_burn_severity),
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
             year(date) == 2020) %>% 
      select(-valid_count) %>% 
      # Flag pre- & post-fire observations
      mutate(
        burn_date = ymd_hms(burn_date),
        BeforeBurnDate = date < burn_date
      )
    
    # Extract topography at random points 
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
    df_list <- list(elevation_sample, slope_sample, northness_sample,eastness_sample)
    
    # merge all data frames in list
    dem_sample <- df_list %>% reduce(full_join, by='ObservationID')
    
    df_filtered <- full_join(df_filtered, dem_sample,by = "ObservationID")
    
    # Merge with LST data
    df_filtered <- df_filtered %>% 
      left_join(df_lst %>% select(ObservationID, Time, index_name), 
                by = c("ObservationID" = "ObservationID", 
                       "date" = "Time")) %>%
      rename(LST = index_name) # Rename index_name to LST
    
    # Write filtered data frame to CSV
    write.csv2(df_filtered,
               sprintf("data/tables/merged_filtered_%sth_pctile.csv",
                       pct_cutoff*100))
  }
}
