library(terra)
library(sf)
library(spatstat)
library(tidyterra)
library(rts)
library(spdep)
library(gt)
library(tidyverse)
library(cowplot)
library(patchwork)
library(lubridate)
print(getwd())

set.seed(10)

# 0. Configure functions ----
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
  spectral_index_s30 <- rast(paste0(hls_dir, spectral_index_files_s30))
  spectral_index_s30 <- sapp(spectral_index_s30, assign_rast_time)
  
  spectral_index_l30 <- rast(paste0(hls_dir, spectral_index_files_l30))
  spectral_index_l30 <- sapp(spectral_index_l30, assign_rast_time)
  
  # Merge S30 and L30 stack
  spectral_index <- c(spectral_index_s30, spectral_index_l30)
  
  cat("done.")
  return(spectral_index)
}

# 1. Loading and preparing data: ----
UTM_TILE_ID <- "54WXE"
year <- 2020
severity_index <- "dNBR"

dnbr <- rast(
  sprintf("~/data/raster/hls/severity_rasters/%s_%s_%s.tif",
          UTM_TILE_ID, year,severity_index)) * 1000

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  "~/data/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg"
)
TEST_ID <- 14211 # fire ID for part of the large fire scar

selected_fire_perimeter <- fire_perimeters %>% 
  filter(fireid  == TEST_ID) %>%
  project(crs(dnbr))

# Buffer the selected fire perimeter
fire_perimeter_buffered <- buffer(selected_fire_perimeter, 1200)

plot(selected_fire_perimeter)
plot(fire_perimeter_buffered,add = TRUE)

binwidth <- 20

dnbr_in_perimeter <- dnbr %>% 
  mask(fire_perimeter_buffered, updatevalue = NA) %>% 
  crop(fire_perimeter_buffered)

p1 <- ggplot() +
  geom_spatraster(data = dnbr_in_perimeter) +
  scale_fill_viridis_c(option = "inferno",
                       na.value = "white",
                       name = "dNBR") +
  theme_cowplot()

p2 <- ggplot(data = dnbr_in_perimeter, aes(x = dNBR)) +
  geom_histogram(binwidth = binwidth, color = "black", fill = "indianred1") +
  xlim(c(-500,1000)) +
  labs(
    title = "Histogram of dNBR",
    x = "dNBR Values",
    y = "Frequency"
  ) +
  theme_cowplot()

pg <- p1 + p2

ggsave2(pg,filename = "figures/expl_dNBR_hist.png",width = 8,height = 6,bg = "white")

# 2. Execute spatial sampling ----
# Set dNBR bin size for sampling (dNBR * 1000)
binwidth <- 20

# get bins for the raster
bins <- seq(minmax(dnbr_in_perimeter)[1],minmax(dnbr_in_perimeter)[2],binwidth)

# Set proportion of sampled values per bin
nsample <- 1e4

# Create binned raster for sampling
rast_binned <- classify(dnbr_in_perimeter,
                        bins, 
                        include.lowest=TRUE,
                        brackets=TRUE)

# calculate number of pixels to sample in each bin (now 10% of available pixels)
frac_to_sample <- 0.01

# Create spatial point sample
sample_points <- sample_dnbr_points(rast_binned, sample_pct = frac_to_sample)

# Convert to sf
sample_points_sf <- st_as_sf(sample_points)

# Export points as gpkg
st_write(sample_points_sf, "~/data/feature_layers/sample_points.gpkg",
         layer = "sample_points", delete_layer = TRUE)

# 3. Run burn date assignment in python ----

# 4. Extract data ----

sample_points <- vect("~/data/feature_layers/sample_points_burn_date.gpkg") %>% 
  project(crs(dnbr)) %>% 
  mutate(ObservationID = 1:nrow(sample_points))

# Extract data at random points 
dnbr_sample <- terra::extract(dnbr, sample_points,
                              ID = FALSE, xy = TRUE) 

# Plot distribution of total raster vs. sample points
p3 <- ggplot(data = dnbr_sample, aes(x = dNBR)) +
  geom_histogram(binwidth = binwidth, color = "black", fill = "orange") +
  xlim(c(-500,1000)) +
  labs(
    title = "Histogram of dNBR (sampled)",
    x = "dNBR Values",
    y = "Frequency"
  ) +
  theme_cowplot()

pg <- p2 / p3 
ggsave2(pg, filename = "figures/expl_dNBR_hist_vs_sampled.png",
        width = 8, height = 6,bg = "white")

# Export dataframe to csv
# write.csv(dnbr_sample, file = "data/tables/dnbr_sample_test_roi.csv")

# 5. Check spatial autocorrelation of Random points ----
knn <- knn2nb(knearneigh(dnbr_sample[c("x","y")], k = 5))

# Convert neighbors to weights
weights <- nb2listw(knn, style = "W")

# Compute Moran's I
moran_i <- moran.test(dnbr_sample$dNBR, listw = weights)

ggplot() +
  geom_spatraster(data = mask(dnbr,fire_perimeter_buffered)) +
  scale_fill_viridis_c(option = "inferno",
                       na.value = "white",
                       name = "dNBR") +
  geom_spatvector(data = sample_points) +
  labs(title = sprintf("Randomly samlped points\nMoran's I: %.2f (p-value: %.4f)", 
                       moran_i$estimate[1],moran_i$p.value)) +                   
  theme_cowplot()

# ggsave2("figures/stratified_sample_points_moranI.png",bg = "white")

# 6. Stack spectral_index rasters ----
HLS_DIR <- "~/data/raster/hls/"
index_name = "NDMI"

image_stack <- stack_time_series(HLS_DIR,UTM_TILE_ID, index_name)

# create raster time series object
rt <- rts(image_stack, time(image_stack))

df_spectral_index <- terra::extract(rt, sample_points) %>%
  t() %>%
  as_tibble()

filename <- sprintf("%s_sampled_%s_%s.csv",index_name,UTM_TILE_ID,year)
write.csv2(df_spectral_index,paste0(HLS_DIR,filename))

# Load NDVI and NDMI time series
read_hls_data_frames <- function(index_name, UTM_TILE_ID,year,
                                 sample_points, dnbr_sample){
  filename <- sprintf("%s_sampled_%s_%s.csv",index_name,UTM_TILE_ID,year)
  
  df <- read.csv2(paste0(HLS_DIR,filename)) %>% 
    select(-1) %>% 
    mutate(dNBR = dnbr_sample$dNBR) %>% 
    as_tibble() %>% 
    mutate(burn_date = sample_points$burn_date)
  
  df_long <- df %>%
    mutate(ObservationID = 1:nrow(df)) %>%  # Add a column for observation IDs (row numbers)
    gather(key = "Time", value = index_name, 
           -c(ObservationID,dNBR,burn_date)) %>% 
    mutate(Time = as.POSIXct(Time, format = "X%Y.%m.%d.%H.%M.%S")) 
  
  return(df_long)
}

df_ndmi <- read_hls_data_frames("NDMI", UTM_TILE_ID,year,sample_points, dnbr_sample)
df_ndvi <- read_hls_data_frames("NDVI",UTM_TILE_ID,year,sample_points, dnbr_sample)

# Concatenate the data frames
df_hls <- df_ndmi %>% 
  mutate(NDVI = df_ndvi$index_name) %>% 
  rename(NDMI = index_name)



df_nobs <- df_hls %>% 
  select(-dNBR) %>% 
  group_by(ObservationID) %>% 
  summarize(valid_count = sum(!is.na(index_name)))

# Compute daily means
df_daily_spectral_index <- df_hls %>% 
  group_by(Time, ObservationID) %>% 
  summarise(DailyMean = mean(index_name, na.rm = TRUE),
            dNBR = first(dNBR),
            burn_date = first(burn_date))

pct_cutoff <- 0.5
thr_nobs <- quantile(df_nobs$valid_count,pct_cutoff , na.rm = TRUE)

ggplot(df_nobs) +
  geom_histogram(aes(x = valid_count)) +
  geom_vline(aes(xintercept = thr_nobs),
             color = "red", linetype = "dashed", size = 1) + 
  geom_text(aes(x = thr_nobs-10, y = 300), 
                label = sprintf("%sth percentile = %s",
                                pct_cutoff*100, round(thr_nobs, 2)),
            color = "red", vjust = -1) +
  labs(title = "Number of non-NA observations",
       subtitle = sprintf("Thresholds of valid data (%sth percentile)",
                          pct_cutoff*100)) +
  theme_cowplot()

ggsave2(sprintf("figures/Histogram_%s_observations_%sth_pctile.png",
                index_name,pct_cutoff*100),
        bg = "white",width = 10, height = 8)

# get nr. of valid observations
valid_counts_by_id <- df_daily_spectral_index %>%
  group_by(ObservationID) %>%
  summarise(valid_count = sum(!is.na(DailyMean)))

# Plot where number of observations per points
sample_points_filtered <- sample_points %>% 
  merge(valid_counts_by_id, by = "ObservationID")

(p <- ggplot() +
    # geom_spatraster(data = dnbr_in_perimeter, aes(fill = dNBR)) +
    geom_spatvector(data = sample_points_filtered, aes(color = valid_count)) +
    scale_color_viridis_c(option = "magma",
                         na.value = "transparent") +
    labs(fill = "Number of non-masked observations") +
    theme_cowplot()
)

ggsave2(p, 
        filename = "figures/Nobs_randompoints_map.png",
        bg = "white",width = 10, height = 8)

# Load filtered and formatted dataframe
df_filtered <- read.csv2(sprintf("data/tables/%s_filtered_%sth_pctile.csv",
                                index_name,pct_cutoff*100)) %>% 
  mutate(Time = as.Date(Time),
         burn_date = as.Date(burn_date))

# Filter out observations with fewer than X observations
df_filtered <- df_daily_spectral_index %>%
  inner_join(valid_counts_by_id, by = "ObservationID") %>%
  filter(valid_count >= thr_nobs,
         format(Time, "%Y") == "2020") %>% 
  select(-valid_count) %>% 
  # Flag pre- & post-fire observations
  mutate(
    Time = ymd_hms(Time),
    burn_date = ymd_hms(burn_date),
    BeforeBurnDate = Time < burn_date
  )

# Plot where these filtered points sit
sample_points_filtered <- sample_points %>% 
  filter(ObservationID %in% df_filtered$ObservationID)

(p <- ggplot() +
  geom_spatraster(data = dnbr_in_perimeter, aes(fill = dNBR)) +
  geom_spatvector(data = sample_points_filtered) +
  scale_fill_viridis_c(option = "magma",
                       na.value = "transparent") +
  labs(fill = "dNBR") +
  theme_cowplot()
)
ggsave2(p, 
        filename = sprintf("figures/Location_randompoints_%sth_pctile.png",pct_cutoff*100),
        bg = "white",width = 10, height = 8)

# Time series for all points
ggplot(df_filtered) +
  geom_point(aes(x = Time,y = DailyMean),alpha = .2) +
  labs(y = sprintf("Daily mean %s\n (per point location)",index_name)) +
  theme_cowplot()

ggsave2(sprintf("figures/Timeseries_%s_randompoints.png",index_name),
        bg = "white",width = 10, height = 8)

# Time series for some random points
rand_id <- sample(unique(df_filtered$ObservationID), 16)

df_filtered %>% 
  filter(ObservationID %in% rand_id) %>% 
  ggplot() +
    geom_point(aes(x = Time,y = DailyMean, color = BeforeBurnDate),alpha = .2) +
    labs(y = sprintf("Daily mean %s\n (per point location)",index_name)) +
    geom_vline(aes(xintercept = as.numeric(burn_date)), color = "red", linetype = "dashed") +
    facet_wrap(~ObservationID) +
    theme_cowplot()

ggsave2(sprintf("figures/%s_perpoint.png",index_name),
            bg = "white",width = 9, height = 8)
    
# Histogram of index values
ggplot(df_filtered) +
  geom_histogram(aes(x = DailyMean)) +
  geom_vline(xintercept = mean(df_filtered$DailyMean,na.rm = T),
             color = "red", linetype = "dashed", size = 1) + 
  geom_text(x = 0.2, y = 1e3, 
                label = paste("Mean = ", 
                              round(mean(df_filtered$DailyMean,na.rm = T), 2)), 
            color = "red", vjust = -1) +
  labs(title = sprintf("Histogram of %s values",index_name)) +
  theme_cowplot()

ggsave2(sprintf("figures/Histogram_%s.png",index_name),
        bg = "white",width = 10, height = 8)

# Scatterplot with dNBR color
ggplot(df_filtered) +
  geom_point(aes(x = DailyMean,y = dNBR),alpha = .2) +
  labs(x = sprintf("Daily mean %s\n (per point location)",index_name),
       y = "dNBR") +
 theme_cowplot()

# Write filtered data frame to CSV
write.csv2(df_filtered,
           sprintf("data/tables/%s_filtered_%sth_pctile.csv",
                  index_name,pct_cutoff*100))

# 7. Run pre-fire linear models ----
df_mod <- df_filtered %>%
  group_by(ObservationID) %>% 
  filter(BeforeBurnDate) %>% 
  do(mod = lm(DailyMean ~ Time, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2],
         Intercept = summary(mod)$coeff[1]) %>%
  select(-mod)

df_mod <- df_filtered %>%
  group_by(ObservationID) %>%
  select(dNBR) %>% 
  summarize(dNBR = first(dNBR)) %>% 
  select(-ObservationID) %>% 
  cbind(df_mod)
    
# Plot slope of index vs. dNBR
ggplot(df_mod) +
  geom_point(aes(x = Slope,y = dNBR),alpha = .2) +
  labs(x = sprintf("Slope of pre-fire %s",index_name),
       y = "dNBR") +
  theme_cowplot()

ggsave2(sprintf("figures/%s_vs_dNBR.png",index_name),
        bg = "white",width = 10, height = 8)

sample_points_merged <- sample_points %>% 
  filter(ObservationID %in% df_mod$ObservationID) %>% 
  merge(., df_mod, by = "ObservationID") 

(p <- ggplot() +
    geom_spatraster(data = dnbr_in_perimeter,aes(fill = dNBR)) +
    geom_spatvector(data = sample_points_merged, aes(color = Slope)) +
    scale_fill_viridis_c(option = "magma",
                         na.value = "transparent") +
    scale_color_viridis_c(option = "viridis",
                         na.value = "transparent") +
    labs(fill = "dNBR",
         color = sprintf("Slope %s vs. dNBR",index_name)) +
    theme_cowplot()
)
ggsave2(p, 
        filename = sprintf("figures/Map_Slope_per_Point_%s_vs_dNBR_%sth.png",
                           pct_cutoff*100,index_name),
        bg = "white",width = 10, height = 8)

# Report summary stats for pre-fire observations
df_filtered %>%
  group_by(ObservationID) %>%
  filter(BeforeBurnDate) %>%
  summarise(
    n_observations = n()
  ) %>%
  summarise(
    mean_obs = mean(n_observations),
    median_obs = median(n_observations),
    min_obs = min(n_observations),
    max_obs = max(n_observations)
  )
