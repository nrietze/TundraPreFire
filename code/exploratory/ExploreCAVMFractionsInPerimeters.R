##################################################
## Project: Chapter 3
## Script purpose: retrieve spatiotemporal distribution of CAVM types in burn perimeters
## Date: September 1, 2024
## Author: Nils Rietze (nils.rietze@uzh.ch)
##################################################

library(terra)
library(tidyterra)
library(tidyverse)
library(cowplot)
# 
## Section: Load and prepare data ----
#.....................................

fire_perimeters <- vect("C:/data/3_fire_data/burned_area/siberia_talucci/data/SiberiaFires/SiberiaFires2001-2020_wgs.shp") 

cavm_perimeter <- vect("C:/data/6_vegetation/cavm/vegetation_zones_wgs84.shp")
cavm <- rast("C:/data/6_vegetation/cavm/raster_cavm_v1.tif") %>% 
  as.factor() %>% 
  project("EPSG:4326")

cavm_legend <- read.csv2("C:/data/6_vegetation/cavm/Raster CAVM legend.csv") %>% 
  select(c("Raster.code","Vegetation.Unit")) %>% 
  mutate(Raster.code = as.factor(Raster.code))

# intersect Talucci-perimeters with CAVM extent
fp_cavm <- terra::intersect(fire_perimeters,cavm_perimeter) %>% 
  mutate(ID = 1:nrow(.))

# intersect CAVM raster with perimeters and bind spatvect
df_cavm <- terra::extract(cavm,fp_cavm,method = "simple",
                          weights = TRUE) %>% 
  rename("Raster.code" = "Band_1")

# get area per class per perimeter
class_summary <- df_cavm %>%
  group_by(ID, Raster.code) %>%
  summarise(class_count = sum(weight), .groups = 'drop') 

# get total area per perimeter
total_summary  <- df_cavm %>%
  group_by(ID) %>%
  summarise(total_count = sum(weight), .groups = 'drop') 

# merge summaries
summary_with_total <- merge(class_summary, total_summary, by = "ID")

# get percentages
summary_with_total <- summary_with_total %>%
  mutate(percentage = (class_count / total_count) * 100) %>% 
  left_join(cavm_legend, by = "Raster.code")

# reshape dataframe
summary_wide <- summary_with_total %>%
  select(-c("Raster.code","class_count","total_count")) %>% 
  pivot_wider(names_from = Vegetation.Unit, values_from = percentage) 

# match with fire perimeters
fire_perimeters_with_percentage <- merge(fp_cavm, 
                                         summary_wide, by = "ID")

df_fire_perimeters_with_percentage <- as.data.frame(fire_perimeters_with_percentage)

df_fire_perimeters_with_percentage %>% 
  select(ID, FireYr, SizeHa, 
         any_of(cavm_legend$Vegetation.Unit[!is.na(cavm_legend$Vegetation.Unit)])) %>% 
  pivot_longer(
    cols = -c(ID, FireYr, SizeHa),  # Keep ID, FireYr, SizeHa as static columns
    names_to = "VegetationUnit",    # New column for Vegetation Unit names
    values_to = "Percentage"        # New column for the percentages
  ) %>% 
  ggplot(aes(x = FireYr,y = Percentage,group = FireYr)) +
  geom_boxplot() + 
  facet_wrap(~VegetationUnit) +
  theme_cowplot()

ggplot(data = df_fire_perimeters_with_percentage, aes(x = SizeHa,
                                                      y = G3,
                                                      color = FireYr)) +
  geom_point() +
  theme_cowplot()
  
  