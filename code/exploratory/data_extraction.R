library(terra)
library(sf)
library(tidyterra)
library(rts)

library(gt)
library(tidyverse)
library(cowplot)
library(patchwork)

print(getwd())

set.seed(10)

## Loading and preparing data:
dnbr <- rast(
  "~/data/raster/hls/dNBR/2020_dNBR.tif"
) * 1000

# Load features (fire perimeters and ROIs)
fire_perimeters_2020 <- vect(
  "~/data/feature_layers/fire_atlas/final_viirs2020.gpkg"
)
roi <- vect("~/data/feature_layers/roi.geojson")

fire_within_roi <- terra::intersect(fire_perimeters_2020, roi)


## Dissolve fire perimeters that lie in the ROI to one
fire_union_roi <- terra::aggregate(fire_within_roi) %>%
  project(crs(dnbr))

fire_perimeter_buffered <- buffer(fire_union_roi, 1200)

plot(fire_union_roi)
plot(fire_perimeter_buffered,add = TRUE)

binwidth <- 20

dnbr_in_perimeter <- mask(dnbr, fire_perimeter_buffered, updatevalue = NA)

p1 <- ggplot() +
  geom_spatraster(data = dnbr_in_perimeter) +
  scale_fill_viridis_c(option = "inferno",
                       na.value = "white",
                       name = "dNBR") +
  theme_cowplot()

p2 <- ggplot(data = dnbr_in_perimeter, aes(x = NBR)) +
  geom_histogram(binwidth = binwidth, color = "black", fill = "indianred1") +
  xlim(c(-500,1000)) +
  labs(
    title = "Histogram of dNBR",
    x = "dNBR Values",
    y = "Frequency"
  ) +
  theme_cowplot()

p1 + p2

# ggsave2("figures/expl_dNBR_hist.png",bg = "white")

## Extract pixels within fire perimeter
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


# Sample n pixels and mask with burn perimeter
random_points <- spatSample(rast_binned, nsample, 
                            # method = "stratified",
                            method = "random",
                            values = FALSE,
                            as.points = TRUE) %>% 
  mask(fire_perimeter_buffered) 

# Export points as gpkg
writeVector(random_points, filename = "~/data/feature_layers/random_points.shp")
random_points <- vect("~/data/feature_layers/random_points.shp")

# Extract data at random points 
dnbr_sample <- terra::extract(dnbr, random_points,
                              ID = FALSE, xy = FALSE) 

p3 <- ggplot(data = dnbr_sample, aes(x = NBR)) +
  geom_histogram(binwidth = binwidth, color = "black", fill = "orange") +
  xlim(c(-500,1000)) +
  labs(
    title = "Histogram of dNBR (sampled)",
    x = "dNBR Values",
    y = "Frequency"
  ) +
  theme_cowplot()

p2 / p3 
# ggsave2("figures/expl_dNBR_hist_vs_sampled.png",bg = "white")

# Export dataframe to csv
# write.csv(dnbr_sample, file = "data/tables/dnbr_sample_test_roi.csv")

## Check spatial autocorrelation of Random points:
library(spdep)

knn <- knn2nb(knearneigh(dnbr_sample[c("x","y")], k = 5))

# Convert neighbors to weights
weights <- nb2listw(knn, style = "W")

# Compute Moran's I
moran_i <- moran.test(dnbr_sample$NBR, listw = weights)

ggplot() +
  geom_spatraster(data = dnbr_in_perimeter) +
  scale_fill_viridis_c(option = "inferno",
                       na.value = "white",
                       name = "dNBR") +
  geom_spatvector(data = random_points) +
  labs(title = sprintf("Randomly samlped points\nMoran's I: %.2f (p-value: %.4f)", 
                       moran_i$estimate[1],moran_i$p.value)) +                   
  theme_cowplot()

# ggsave2("figures/stratified_random_points_moranI.png",bg = "white")

## Stack NDMI rasters
HLS_DIR <- "~/data/raster/hls/"

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

# List NDMI COGs
ndmi_files_s30 <- list.files(HLS_DIR,
                             pattern = "S30.T54WXE.*NDVI\\.tif$"
)
ndmi_files_l30 <- list.files(HLS_DIR,
                             pattern = "L30.T54WXE.*NDVI\\.tif$"
)

ndmi_s30_list <- lapply(paste0(HLS_DIR, ndmi_files_s30),rast)

ndmi_s30 <- rast(paste0(HLS_DIR, ndmi_files_s30[1:(length(ndmi_files_s30)-2)]))
ndmi_s30 <- sapp(ndmi_s30, assign_rast_time)

ndmi_l30_list <- lapply(paste0(HLS_DIR, ndmi_files_l30),rast)

ndmi_l30 <- rast(paste0(HLS_DIR, ndmi_files_l30))
ndmi_l30 <- sapp(ndmi_l30, assign_rast_time)
ndmi_l30 <- project(ndmi_l30, crs(ndmi_s30))

ndmi <- c(ndmi_s30, ndmi_l30)
ndmi

# create raster time series object
rt <- rts(ndmi, time(ndmi))

df_ndmi <- terra::extract(rt, random_points) %>%
  t() %>%
  as_tibble()

summary(df_ndmi)

write.csv2(df_ndmi,paste0(HLS_DIR,"ndvi_sampled.csv"))

# Load data sample
index_name <- "NDMI"

df <- read.csv2(sprintf("%s%s_sampled.csv",HLS_DIR,tolower(index_name))) %>% 
  select(-1) %>% 
  mutate(dNBR = dnbr_sample$NBR) %>% 
  as_tibble()

df_long <- df %>%
  mutate(ObservationID = 1:nrow(df)) %>%  # Add a column for observation IDs (row numbers)
  gather(key = "Time", value = index_name, -c(ObservationID,dNBR)) %>% 
  mutate(Time = as.POSIXct(Time, format = "X%Y.%m.%d.%H.%M.%S")) 

df_nobs <- df_long %>% 
  select(-dNBR) %>% 
  group_by(ObservationID) %>% 
  summarize(valid_count = sum(!is.na(index_name)))

thr_nobs <- quantile(df_nobs$valid_count, 0.75, na.rm = TRUE)

ggplot(df_nobs) +
  geom_histogram(aes(x = valid_count)) +
  geom_vline(aes(xintercept = thr_nobs),
             color = "red", linetype = "dashed", size = 1) + 
  geom_text(aes(x = thr_nobs-10, y = 300, 
                label = paste("75%-ile = ", round(thr_nobs, 2))), 
            color = "red", vjust = -1) +
  geom_vline(xintercept = 170,
             color = "blue", linetype = "dashed", size = 1) +  
  geom_text(aes(x = 170, y = 200, 
                label = "170"), 
            color = "blue", vjust = -1) +
  labs(title = "Number of non-NA observations",
       subtitle = "Thresholds of valid data (blue = visually, red = 75% percentile)") +
  theme_cowplot()

ggsave2(sprintf("figures/Histogram_%s_observations.png",index_name),
        bg = "white",width = 10, height = 8)

# Compute daily means
df_daily_ndmi <- df_long %>% 
  group_by(Time, ObservationID) %>% 
  summarise(DailyMean = mean(index_name, na.rm = TRUE),
            dNBR = first(dNBR))

# get nr. of valid observations
valid_counts_by_id <- df_daily_ndmi %>%
  group_by(ObservationID) %>%
  summarise(valid_count = sum(!is.na(DailyMean)))

# Filter out observations with fewer than X observations
df_filtered <- df_daily_ndmi %>%
  inner_join(valid_counts_by_id, by = "ObservationID") %>%
  filter(valid_count >= thr_nobs,
         format(Time, "%Y") == "2020") %>% 
  select(-valid_count)  

# Time series for all points
ggplot(df_filtered) +
  geom_point(aes(x = Time,y = DailyMean),alpha = .2) +
  labs(y = sprintf("Daily mean %s\n (per point location)",index_name)) +
  theme_cowplot()

ggsave2(sprintf("figures/Timeseries_%s_randompoints.png",index_name),
        bg = "white",width = 10, height = 8)

# Time series for 64 random points
rand_id <- sample(unique(df_filtered$ObservationID), 64)

df_filtered %>% 
  filter(ObservationID %in% rand_id) %>% 
  ggplot() +
    geom_point(aes(x = Time,y = DailyMean),alpha = .2) +
    labs(y = sprintf("Daily mean %s\n (per point location)",index_name)) +
    facet_wrap(~ObservationID) +
    theme_cowplot()

ggsave2(sprintf("figures/%s_perpoint.png",index_name),
            bg = "white",width = 16, height = 16)
    
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
 ¨¨ theme_cowplot()
ä