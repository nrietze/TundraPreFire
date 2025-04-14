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


# 1. Loading and preparing data: ----
# ===================================.

# Set name of severity index to extract data for
burn_severity_index <- "dNBR"

# TRUE to overwrite existing data form time series extraction
OVERWRITE_DATA <- TRUE

# TRUE to create random point samples on dNBR maps
SAMPLE_FROM_DNBR_RASTER <- TRUE

# Define percentile for sample cutoff
pct_cutoff <- 0.5

# Set proportion of sampled values per bin
frac_to_sample <- 0.01
frac_int <- frac_to_sample *100

OS <- Sys.info()[['sysname']]

# Output directory for sample tables
TABLE_DIR <- ifelse(OS == "Linux", 
                    "~/data/tables/","data/tables/")

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)

top20_fires <- fire_perimeters %>%
  arrange(desc(farea)) %>% 
  slice_head(n = 20) 

# TEST_ID <- top20_fires$fireid
TEST_ID <- c(14211,14664,10792,17548) # fire ID for part of the large fire scar

# Load lookup tables
final_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017) 

optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

if (length(TEST_ID) > 0){final_lut <- filter(final_lut,fireid %in% TEST_ID)}

dem_lut <- read.csv(paste0(TABLE_DIR,"dem_fire_perim_intersect.csv")) # DEM tiles

# 2. Execute spatial sampling ----
# =================================.

# Apply data extraction to all fire perimeters
for(i in 1:nrow(final_lut)) {
  row <- final_lut[i, ]
  
  # extract fire peimeter attributes
  UTM_TILE_ID <- row$opt_UTM_tile
  year <- row$tst_year
  FIRE_ID <- row$fireid
  
  print(sprintf("Extracting data for fire %s in UTM tile: %s",FIRE_ID,UTM_TILE_ID))
  
  # Load optimal burn severity raster
  fname_optimal_severity_raster <- optimality_lut %>% 
    filter(fireid == FIRE_ID,
           severity_index == burn_severity_index) %>% 
    pull(fname_severity_raster)
  rast_burn_severity <- rast(fname_optimal_severity_raster)
  
  # get fire perimeter
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
  
  fname_dem_out <- paste0("~/data/raster/arcticDEM/",
                          FIRE_ID,
                          "_dem_30m.tif")
  
  if (!file.exists(fname_dem_out)){
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
    writeRaster(dem_30m_utm,filename = fname_dem_out,overwrite = T)
    
    rm(dem_30m,dem_30m_utm)
  }
  
  # Buffer the selected fire perimeter
  fire_perimeter_buffered <- buffer(selected_fire_perimeter, 1200)
  
  dnbr_in_perimeter <- rast_burn_severity %>% 
    mask(fire_perimeter_buffered, updatevalue = NA) %>% 
    crop(fire_perimeter_buffered) * 1000
  
  fname_sample_points <- sprintf("~/data/feature_layers/%s_sample_points_%spct.gpkg",
                                 FIRE_ID,frac_int)
  
  # Sample points only when gpkg file doens't exist
  if (!file.exists(fname_sample_points) || OVERWRITE_DATA){
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
    
    # Create spatial point sample
    if (SAMPLE_FROM_DNBR_RASTER){
      sample_points <- sample_dnbr_points(rast_binned, sample_pct = frac_to_sample)
    } else {
      # Weigh nr. of sampled points by area burned
      n_points <- round(100 * fire_perimeter_buffered$farea)
      
      sample_points <- spatSample(fire_perimeter_buffered,
                                  size = n_points,
                                  method = "random")
    }
    
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
# run in part2 of the script