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
read_hls_data_frames <- function(index_name, FIRE_ID,year,
                                 sample_points, dnbr_sample){
  filename <- sprintf("%s_sampled_%s_%s.csv",index_name,FIRE_ID,year)
  
  df <- read.csv2(paste0(TABLE_DIR,filename)) %>% 
    select(-1) %>% 
    mutate(dNBR = dnbr_sample$dNBR) %>% 
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
           -c(ObservationID,dNBR,burn_date,descals_burn_class)) %>% 
    mutate(Time = as.POSIXct(Time, format = fmt)) 
  
  return(df_long)
}

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
    "~/data/tables/%s_%s_merged_filtered_%sth_pctile.csv",
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

PlotSeverityMapHist <- function(cropped_severity_raster){
  severity_index <- names(cropped_severity_raster)
  SCALE_FACTOR <- ifelse(severity_index == "dNBR", 1000, 1)
  
  binwidth <- .02 * SCALE_FACTOR
  
  # Plot dNBR map
  p1 <- ggplot() +
    geom_spatraster(data = cropped_severity_raster) +
    scale_fill_viridis_c(option = "inferno",
                         na.value = "white",
                         name = severity_index) +
    theme_cowplot()
  
  # Plot histogram for this map
  xvar <- as.symbol(severity_index)
  p2 <- ggplot(data = cropped_severity_raster, aes(x = !!xvar)) +
    geom_histogram(binwidth = binwidth, color = "black", fill = "indianred1") +
    xlim(c(-.5*SCALE_FACTOR,1*SCALE_FACTOR)) +
    labs(
      title = sprintf("Histogram of %s",severity_index),
      x = sprintf("%s Values",severity_index),
      y = "Frequency"
    ) +
    theme_cowplot()
  
  pg <- p1 + p2
  return(pg)
}

# 1. Configure and load stuff ----
# ================================.
# Config, loading and preparing data
HLS_DIR <- "~/data/raster/hls/"
TABLE_DIR <- "~/data/tables/"

FIRE_ID <- 14664
index_name <- "NDMI"
severity_index <- "dNBR"
pct_cutoff <- 0.5
SAVE_FIGURES <- TRUE

# Load lookup table
final_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017)

# Load fire perimeters
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)

if (!is.na(FIRE_ID)){
  final_lut <- filter(final_lut, fireid == FIRE_ID)
}

data_list <- load_data(final_lut, severity_index)

df_filtered <- data_list[[1]]
severity_raster_sample <- data_list[[2]]
sample_points <- data_list[[3]]
cropped_severity_raster <- data_list[[4]]
selected_fire_perimeter <- data_list[[5]]

rm(data_list)

# 2. Descriptive statistics ----
# ==============================.

# Map & Histogram of burn severity data
pg <- PlotSeverityMapHist(cropped_severity_raster)
if (SAVE_FIGURES){
  ggsave2(pg,filename = sprintf("figures/%s_%s_map_hist.png",severity_index,FIRE_ID),
          width = 8,height = 6,bg = "white")
}

