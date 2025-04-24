library(terra)
library(sf)
library(tidyterra)
library(rts)
library(readr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(lubridate)
library(tictoc)

set.seed(10)

# 0. Configure functions ----
# ===========================.
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

get_ls_datetime <- function(raster){
  datenum <- global(raster,"max",na.rm=T)
  date <- as.Date(as.numeric(datenum), origin = "1970-01-01")
  return(date)
}

# 1. Loading and preparing data: ----
# ===================================.

# TRUE to overwrite existing data form time series extraction
OVERWRITE_DATA <- TRUE

# Define percentile for sample cutoff
pct_cutoff <- 0.5

frac_to_sample <- 0.01
frac_int <- frac_to_sample *100

OS <- Sys.info()[['sysname']]

# Output directory for sample tables
TABLE_DIR <- ifelse(OS == "Linux","~/data/tables/","data/tables/")
OUT_DIR <- paste0(TABLE_DIR,"sampled_data/",frac_int,"pct/")
HLS_DIR <- normalizePath("~/scratch/raster/hls/processed/")

dir.create(OUT_DIR, showWarnings = FALSE)

# Load lookup tables
dem_lut <- read.csv(paste0(TABLE_DIR,"dem_fire_perim_intersect.csv")) # DEM tiles

optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

processing_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017) %>% 
  arrange(by = opt_UTM_tile)

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)

FALSE_FIRES_ID <- c(23633,21461,15231,15970,17473,13223,
                    14071,12145,10168,24037,13712)

topN_fires <- fire_perimeters %>%
  filter(tst_year >= 2017) %>% 
  arrange(desc(farea)) %>% 
  slice_head(n = 25) 

# TEST_ID <- c(14664,10792,17548) 
TEST_ID <- topN_fires$fireid

if (length(TEST_ID) > 0){final_lut <- filter(processing_lut,fireid %in% TEST_ID)}


# Run data_extraction_part1.R first

# 2. Run burn date assignment in python ----
# ==========================================.
# Run on python instance on S3IT

# ...or use below if on local machine

# give permission: chmod +x code/data_processing/run_burn_date_extraction.sh
# system("code/data_processing/run_burn_date_extraction.sh")

# 3. Stack rasters and extract data----
# =====================================.
UTM_TILE_ID_old <- 9999 # dummy ID to start with

for(i in 22:nrow(final_lut)) {
  row <- final_lut[i, ]
  
  # extract fire perimeter attributes
  UTM_TILE_ID <- row$opt_UTM_tile
  year <- row$tst_year
  FIRE_ID <- row$fireid
  
  if (FIRE_ID %in% FALSE_FIRES_ID){
    next
  }
  
  cat(sprintf("Extracting data for fire %s in UTM tile: %s \n",
              FIRE_ID,UTM_TILE_ID))
  
  # Load optimal burn severity raster
  fname_optimal_dnbr_raster <- optimality_lut %>% 
    filter(fireid == FIRE_ID,severity_index == 'dNBR') %>% 
    pull(fname_severity_raster)
  
  rdnbr_path <- gsub("dNBR", "RdNBR", fname_optimal_dnbr_raster)
  rbr_path <- gsub("dNBR", "RBR", fname_optimal_dnbr_raster)
  
  if (length(rbr_path) == 0){
    cat("missing burn severity files. Skipping data extraction. \n")
    next
  }
  
  fname_optimal_dgemi_raster <- optimality_lut %>% 
    filter(fireid == FIRE_ID,severity_index == 'dGEMI') %>% 
    pull(fname_severity_raster)
    
  rast_dnbr <- rast(fname_optimal_dnbr_raster)
  rast_rdnbr <- rast(rdnbr_path)
  rast_rbr <- rast(rbr_path)
  rast_dgemi <- rast(fname_optimal_dgemi_raster)
  
  # Load sample points with burn dates
  fname_sample_points <- sprintf("~/data/feature_layers/%spct/%s_sample_points_%spct_burn_date.gpkg",
                                 frac_int,FIRE_ID,frac_int)
  sample_points <- vect(fname_sample_points) %>% 
    project(crs(rast_dnbr)) %>% 
    mutate(ObservationID = 1:nrow(.))
  
  # Extract dNBR at random points 
  dnbr_sample <- terra::extract(rast_dnbr, sample_points,ID = FALSE, xy = TRUE)
  rdnbr_sample <- terra::extract(rast_rdnbr, sample_points,ID = FALSE, xy = TRUE)
  rbr_sample <- terra::extract(rast_rbr, sample_points,ID = FALSE, xy = TRUE)
  dgemi_sample <- terra::extract(rast_dgemi, sample_points,ID = FALSE, xy = TRUE)
  
  ## a. HLS spectral indices ----
  
  for (index_name in c("NDMI","NDVI")){
    # Check if filename exists and reprocessing is activated
    filename <- sprintf("%s_sampled_%s_%s.csv",index_name,FIRE_ID,year)
    
    fmt  <- ifelse(index_name == "NDMI","%Y-%m-%d %H:%M:%S","X%Y.%m.%d.%H.%M.%S")
    
    # only stack rasters if needed (takes a while)
    if (!file.exists(paste0(OUT_DIR,filename)) || OVERWRITE_DATA){
      
      # use previous stacks, if same UTM tile
      if ((UTM_TILE_ID == UTM_TILE_ID_old) &
          exists("ndvi_rt") &
          exists("ndmi_rt") ){
        # If still the same UTM tile, call old image stacks
        if ((index_name == "NDMI") & (crs(ndmi_rt[[1]])==crs(rast_dnbr))){
          rast_series <- ndmi_rt
        } else if ((index_name == "NDVI") & (crs(ndvi_rt[[1]])==crs(rast_dnbr))){
          rast_series <- ndvi_rt
        }
      } else{
        cat(sprintf("Stacking and extracting %s data... \n",index_name))

        # List spectral_index COGs
        spectral_index_files_s30 <- list.files(HLS_DIR,
                                               pattern = sprintf("HLS.S30.T%s.%s.*%s\\.tif$",
                                                                 UTM_TILE_ID,year,index_name),
                                               full.names = TRUE)
        spectral_index_files_l30 <- list.files(HLS_DIR,
                                               pattern = sprintf("HLS.L30.T%s.%s.*%s\\.tif$",
                                                                 UTM_TILE_ID,year,index_name),
                                               full.names = TRUE)
        
        # Concatenate as single raster and assign time dimension
        spectral_index_s30 <- rast(spectral_index_files_s30)
        spectral_index_s30 <- sapp(spectral_index_s30, assign_rast_time)
        
        spectral_index_l30 <- rast(spectral_index_files_l30)
        spectral_index_l30 <- sapp(spectral_index_l30, assign_rast_time)
        
        # Merge S30 and L30 stack
        spectral_index <- c(spectral_index_s30, spectral_index_l30)
        
        # create raster time series object
        cat("Converting to time series.")
        rast_series <- rts(spectral_index, time(spectral_index))

        if (index_name == "NDMI"){
          ndmi_rt <- rast_series
        } else if (index_name == "NDVI"){
          ndvi_rt <- rast_series
        }
      }

      # Extract spectral index time series at each sample point
      tic()
      df_spectral_index <- terra::extract(rast_series, sample_points) %>%
        t() %>%
        as_tibble()
      toc()
      
      df_long <- df_spectral_index %>% 
        mutate(ObservationID = 1:nrow(.),
               dnbr = dnbr_sample[,1],
               rdnbr = rdnbr_sample[,1],
               rbr = rbr_sample[,1],
               dgemi = dgemi_sample[,1],
               burn_date = sample_points$burn_date,
               descals_burn_class = sample_points$descals_burned,
               .before = 1) %>% 
        gather(key = "Time", value = !!sym(index_name), 
               -c(ObservationID,dnbr,rdnbr,rbr,dgemi,burn_date,descals_burn_class)) %>% 
        mutate(Time = ymd_hms(Time))
      
      # Export as table
      cat("Writing data to csv...\n")
      write_csv2(df_long,paste0(OUT_DIR,filename))
      
    } 
  }
  
  # Load  NDMI time series at points 
  cat("Loading NDMI samples...\n")
  filename <- sprintf("NDMI_sampled_%s_%s.csv",FIRE_ID,year)
  
  df_ndmi <- read_csv2(paste0(OUT_DIR,filename), col_names = TRUE,show_col_types = FALSE) %>% 
    as_tibble()
  df_ndmi <- df_ndmi %>% 
    mutate(Time = ymd_hms(Time))
  
  # Load NDVI time series at points
  cat("Loading NDVI samples...\n")
  filename <- sprintf("NDVI_sampled_%s_%s.csv",FIRE_ID,year)
  
  df_ndvi <- read_csv2(paste0(OUT_DIR,filename), col_names = TRUE,show_col_types = FALSE) %>% 
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
                            pattern = paste0("^", FIRE_ID,"_",UTM_TILE_ID, ".*\\.tif$"),
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
    lst <- rast(lst_list) %>% project(crs(rast_dnbr))
    
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
             burn_date = sample_points$burn_date,
             descals_burn_class = sample_points$descals_burned,
             .before = 1) %>% 
      gather(key = "Time", value = "LST", 
             -c(ObservationID,burn_date,descals_burn_class)) %>% 
      mutate(Time = as_datetime(Time, format = fmt))
    
    # Export as table
    cat("Writing data to csv...\n")
    write_csv2(df_lst,paste0(OUT_DIR,filename))
  } else {
    
    cat("Loading LST samples...\n")
    df_lst <- read_csv2(paste0(OUT_DIR,filename), col_names = TRUE,show_col_types = FALSE) %>% 
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
                dnbr = first(dnbr),
                rdnbr = first(rdnbr),
                rbr = first(rbr),
                dgemi = first(dgemi),
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
  UTM_TILE_ID_old <- UTM_TILE_ID
}
