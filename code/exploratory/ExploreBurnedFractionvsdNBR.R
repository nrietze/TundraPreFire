library(terra)
library(landscapemetrics)
library(tidyterra)
library(tidyverse)
library(cowplot)
library(scico)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# 1. Load data ----

# PlanetScope Burned area & fractions
DIR_BF <- "C:/Users/nrietze/Documents/1_PhD/8_CHAPTER2/SiberiaFires/data/geodata/raster"

FNAME_BA <- paste(DIR_BF, "burned_area/planet/LargeScarCenter_burned_area_top5TD.tif",sep = "/")
burned_area <- rast(FNAME_BA)

# PlanetScope top 5 predictors for burned area
FNAME_BA_PR <- paste(DIR_BF, "burned_area/planet/LargeScarCenter_predictors_top5TD.tif",sep = "/")
r_p <- rast(FNAME_BA_PR)
GEMI <- r_p$gemi

# Sentinel-2 dNBR
dNBR_s2 <- rast("data/raster/sentinel2/melanie_hodel/2020_dNBR_BA_Pixelremoved.tif") %>% 
  crop(ext(burned_area))
names(dNBR_s2) <- "dNBR"

burn_severity <- rast("data/raster/sentinel2/melanie_hodel/Classified_2020_BA.tif") %>% 
  crop(ext(burned_area))
names(burn_severity) <- "burn_severity"

# compute 20 m burned fractions
window_side_length <- res(dNBR_s2)[1]
burned_fraction <- ifel(burned_area == "burned", 1, 0) %>% 
  resample(.,dNBR_s2,"sum") / (window_side_length / res(burned_area)[1])**2
names(burned_fraction) <- "burned_fraction"

# 20m GEMI
r_p_20m <- r_p %>% resample(.,dNBR_s2,'bilinear')
GEMI_20m <- GEMI %>% resample(.,dNBR_s2,'bilinear')

# plot burned fractions
plot(burned_fraction, 
     col = scico::scico(20,palette = "lipari",direction = -1))
# plot dNBR
plot(dNBR_s2, 
     col = scico::scico(20,palette = "managua",direction = -1))
# plot burn severity classes
plot(burn_severity, 
     col = scico::scico(3,palette = "hawaii"))

# Plot dNBR vs. burned fraction
n <- 2e4
r <- c(dNBR_s2,burned_fraction,burn_severity,r_p_20m)
df <- spatSample(r,n) %>% 
  mutate(burn_severity = factor(burn_severity,
                                labels = c("unburned",
                                           "high severity",
                                           "low/moderate severity"))) %>% 
  as_tibble()

# scatterplot dNBR vs. burned fraction
ggplot(data = df, aes(x = dNBR, y = burned_fraction)) +
  geom_point(size = .4, alpha = .5) +
  theme_cowplot()

# scatterplot dNBR vs. GEMI
ggplot(data = df, aes(x = dNBR, y = gemi)) +
  geom_point(size = .4, alpha = .5) +
  theme_cowplot()

# boxplot burn class vs. burned fraction
ggplot(data = df,
       aes(x = burn_severity, y = burned_fraction, fill = burn_severity)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                   scale = "width",
                   size = 0.2,alpha = 0.8) +
  geom_point(aes(y = burned_fraction, color = burn_severity),
             position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
  geom_boxplot(lwd = 0.3, width = .2,outlier.shape = NA, alpha = 0.6) +
  labs(y = NULL, x = NULL) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = c(0,0)) +
  scale_fill_manual(values = scico(3, palette = "hawaii")) +
  scale_color_manual(values = scico(3, palette = "hawaii")) +
  theme_cowplot() + 
  theme(legend.position = 'none',
        aspect.ratio=1) 


## Section:  get high severity patch distribution ----
#.................................................
# burn_severity_2019 <- rast("data/raster/sentinel2/melanie_hodel/Cla")
burn_severity_2020 <- rast("data/raster/sentinel2/melanie_hodel/Classified_2020_BA.tif")
names(burn_severity_2020) <- "burn_severity"

check_landscape(burn_severity_2020)

# Compute patch areas
patch_area <- lsm_p_area(burn_severity_2020==1)

# Convert to square metres
patch_area$area_m2 <- patch_area$value * 1e4
ggplot(data = patch_area) +
  geom_histogram(aes(x = area_m2)) +
  scale_x_log10() +
  theme_cowplot()
