library(terra)
library(tidyterra)

# Create list of dNBR files for the test tile
UTM_TILE_ID <- "60VWR"
year <- 2020
severity_index <- "dNBR"

HLS_DIR <- "~/data/raster/hls/severity_rasters/"

burn_severity_filelist <- list.files(HLS_DIR,
                                     pattern = sprintf("^%s_%s_%s-.*\\.tif$",
                                                       severity_index,UTM_TILE_ID,year))

# Load and stack all dNBR rasters for that tile & year
severity_raster_list <- lapply(paste0(HLS_DIR, burn_severity_filelist),rast)

severity_stack <- rast(severity_raster_list)
severity_stack

# Load lookup table
final_lut <- read.csv("data/tables/processing_LUT.csv")

TEST_ID <- final_lut %>% 
  filter(opt_UTM_tile==UTM_TILE_ID,
         ted_year == year) %>% 
  pull(fireid)

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)
selected_fire_perimeter <- fire_perimeters %>% 
  filter(fireid  == TEST_ID[2]) %>%
  project(crs(severity_stack))

plot(severity_stack[[14]] %>% crop(selected_fire_perimeter))
plot(selected_fire_perimeter, add = T)
