library(cowplot)
library(viridis)
library(colorspace)
library(tidyverse)
library(scales)
library(ggh4x)
library(latex2exp)
library(terra)
library(tidyterra)
library(spatstat.utils)
library(grid)
library(patchwork)
library(readr)

set.seed(10)

# 1. Load data and simulate NDMI time series ----
severity_index <- "dNBR"
FIRE_ID <- 14211

OS <- Sys.info()[['sysname']]
if(OS == "Linux"){
  
  HLS_DIR <- "~/scratch/raster/hls/"
  TABLE_DIR <- "~/data/tables/"
  DATA_DIR <- "~/data/"
} else {
  
  HLS_DIR <- "data/raster/hls/"
  TABLE_DIR <- "data/tables/"
  DATA_DIR <- "data/"
}


# Output directory for sample tables
OUT_DIR <- paste0(TABLE_DIR,"sampled_data/")

# Load lookup tables
final_lut <- read.csv(paste0(DATA_DIR,"tables/processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017) 

fire_attributes <- final_lut %>% 
  filter(fireid %in% FIRE_ID)

UTM_TILE_ID <- fire_attributes$opt_UTM_tile
year <- fire_attributes$tst_year

# Load real data from fire scar (fireid = 14211)
data_all <- read_csv2(
  paste0(DATA_DIR,
         sprintf("tables/model_dataframes/1pct/%s_model_dataframe.csv",FIRE_ID,year))
  )

# Load burn severity raster
optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

fname_optimal_severity_raster <- optimality_lut %>% 
  filter(fireid == FIRE_ID,
         severity_index == "dNBR") %>% 
  pull(fname_severity_raster)

rast_burn_severity <- rast(fname_optimal_severity_raster) * SCALE_FACTOR

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  paste0(DATA_DIR,"feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")
)

selected_fire_perimeter <- fire_perimeters %>% 
  filter(fireid  == FIRE_ID) %>%
  project(crs(rast_burn_severity))

dnbr_in_perimeter <- rast_burn_severity %>% 
  mask(selected_fire_perimeter, updatevalue = NA) %>% 
  crop(selected_fire_perimeter)

# Load sample points (with associated burn dates)
fname_sample_points <- paste0(DATA_DIR,
                              sprintf("feature_layers/1pct/%s_sample_points_1pct_burn_date.gpkg",
                               FIRE_ID))

sample_points <- vect(fname_sample_points) %>% 
  project(crs(rast_burn_severity)) %>% 
  mutate(ObservationID = 1:nrow(.))

# Define x-axis: Days since fire (from -50 to +30)
days_since_burn <- seq(-50, 30, by = 1)  # Pre-fire (-50 to 0), Post-fire (0 to 30)

# Function to generate smooth polynomial NDMI curves
simulate_ndmi <- function(days_since_burn, type) {
  ndmi <- numeric(length(days_since_burn))
  
  if (type == "a") {
    # Upward curve before fire, then drop & stabilize
    ndmi[days_since_burn < 0] <- 0.15 + 0.0012 * days_since_burn[days_since_burn < 0]^2 + 0.03 * days_since_burn[days_since_burn < 0]
    ndmi[days_since_burn >= 0] <- 0.12 + 0.05 * exp(-0.05 * days_since_burn[days_since_burn >= 0])
  }
  
  if (type == "b") {
    # Subtle curved stability before fire, then drop & stabilize
    ndmi[days_since_burn < 0] <- 0.18 - 0.0005 * days_since_burn[days_since_burn < 0]^2
    ndmi[days_since_burn >= 0] <- 0.14 + 0.03 * exp(-0.05 * days_since_burn[days_since_burn >= 0])
  }
  
  if (type == "c") {
    # Downward curve before fire, then drop & stabilize
    ndmi[days_since_burn < 0] <- 0.25 - 0.0015 * days_since_burn[days_since_burn < 0]^2 - 0.02 * days_since_burn[days_since_burn < 0]
    ndmi[days_since_burn >= 0] <- 0.06 + 0.03 * exp(-0.05 * days_since_burn[days_since_burn >= 0])
  }
  
  # Constrain NDMI values within -0.2 to 0.4
  ndmi <- pmax(-0.1, pmin(ndmi, 0.3)) +1
  
  return(data.frame(DaysSinceBurn = days_since_burn, NDMI = ndmi, Curve = type))
}

# Generate NDMI curves for each fire response type
ndmi_data <- bind_rows(
  simulate_ndmi(days_since_burn, "a"),
  simulate_ndmi(days_since_burn, "b"),
  simulate_ndmi(days_since_burn, "c")
)

ggplot(ndmi_data,aes(x = DaysSinceBurn,y = NDMI)) +
  geom_point() + 
  theme_cowplot()

# Create bars at days -4 to -1
bars <- ndmi_data %>% 
  group_by(Curve) %>% 
  filter(DaysSinceBurn >= -10 & DaysSinceBurn < 0) %>% 
  select(DaysSinceBurn, NDMI, Curve) %>%
  ungroup() %>% 
  mutate(ymin = 0)

# 1. Plot Figure 1a) - Plot map of fire scar ----
FONT_SIZE <- 18

point_index <- c(2025,278,190)

plot_points <- sample_points %>% 
  filter(ObservationID %in% point_index) %>% 
  mutate(Curve = c("a","b","c"))

fig_1a <- ggplot() +
  geom_spatraster(data = dnbr_in_perimeter) +
  geom_spatvector(data = plot_points,
                  aes(color = Curve), size = 7,show.legend = FALSE) +
  scale_fill_continuous_diverging(palette = "Broc", 
                                   na.value = "white",rev = FALSE,
                                  name = "Burn Severity (dNBR)") +
  scale_color_viridis_d(end = 0.8, option = "inferno") +
  labs(title = "a) Sampling data in fire perimeter.") +
  theme_map(FONT_SIZE)

# 2. Plot Figure 1b) - NDMI curves ----
xlims <- c(-14,4)
curve_names <- c(
  "a" = "Low \nmoisture levels",
  "b" = "Moderate \nmoisture levels",
  "c" = "High \nmoisture levels"
)

fig_1b <- ggplot(ndmi_data, aes(x = DaysSinceBurn, y = NDMI, color = Curve)) +
  geom_line(linewidth = 1.2) +  # Smooth polynomial curves
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",linewidth = 1) +  # Fire event at day 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black",linewidth = 1) +  # NDMI baseline
  geom_vline(xintercept = -10, linetype = "dashed",
             color = "grey40",alpha = 0.7) +  
  # Bars for -4 to -1 days
  geom_rect(data = bars, aes(xmin = DaysSinceBurn - 0.5, xmax = DaysSinceBurn + 0.5, 
                             ymin = ymin, ymax = NDMI, fill = Curve), 
            color = NA, alpha = 0.7) +
  # geom_point(data = bars, aes(x = DaysSinceBurn,y = NDMI, color = Curve),
  #            size = 4 , shape = 18) +
  labs(title = "b) Spline fitting and summing up...",
       x = "Days since burn", y = "Metric for vegetation condition (here NDMI)") +
  scale_color_viridis_d(end = 0.8, option = "inferno") +
  scale_fill_viridis_d(end = 0.8, option = "inferno") +
  lims(x = xlims, 
       #y = c(-0.2, 0.4)
       ) +
  facet_grid(. ~ Curve,labeller = as_labeller(curve_names)) +
  theme_cowplot(FONT_SIZE) +
  theme(legend.position = "none")

annotation <- annotate("label", x = -5, y = 0.25, 
                       label = TeX("$S_{NDMI,10} = \\sum_{t=-10}^{10} {NDMI_t}$",
                                   output = "character"),
                       parse = TRUE,
                       size = 5, color = "#932667FF")

(fig_1b <- fig_1b + at_panel(annotation, PANEL == 2))

# 4. Plot Figure 1c) - cumulative sum of NDMI curves ----
point_vals <- bars %>% 
  group_by(Curve) %>%  
  summarize(NDMI_cumsum = sum(NDMI)) %>% 
  mutate(dNBR = c(450,150,10))

fig_1c <- ndmi_data %>% 
  filter(DaysSinceBurn < 0) %>% 
  group_by(Curve) %>% 
  mutate(NDMI_cumsum = revcumsum(NDMI)) %>% 
  ungroup() %>% 
  ggplot(aes(x = DaysSinceBurn, y = NDMI_cumsum, color = Curve)) +
  geom_point(data = point_vals,aes(x = -10, y = NDMI_cumsum),
             shape = 19,size = 6, alpha = 0.7) +
  # geom_line(linewidth = 1.2) +  
  geom_point(size = 2) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",
             linewidth = 1) +  # Fire event at day 0
  geom_vline(xintercept = -10, linetype = "dashed",
             color = "grey40",alpha = 0.7) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
  labs(title = "c) ...pre-fire vegetation conditions.",
       x = "Days since burn",
       y = TeX("Cumulative NDMI ($S_{NDMI}$)")) +
  scale_color_viridis_d(end = 0.8, option = "inferno") +
  scale_fill_viridis_d(end = 0.8, option = "inferno") +
  lims(x = c(xlims[1],0),
       y = c(-1,20)
       ) +
  theme_cowplot(FONT_SIZE) +
  theme(legend.position = "none")

# 5. Plot Figure 1d) - scatterplot for 10d pre fire ----
fig_1d <- ggplot(data_all,aes(x = NDMI.d_prefire_10 ,y = dnbr*1000)) +
  geom_point(alpha = .2, color = "gray30") +
  geom_point(data = point_vals, aes(x = NDMI_cumsum, y = dNBR, color = Curve),
             shape = 19,size = 6, alpha = 0.7)+
  geom_smooth(method='lm',se = FALSE, color = "black") +
  scale_color_viridis_d(end = 0.8, option = "inferno") +
  labs(title = "d) Testing relationship between vegetation condition and burn severity.",
       x = TeX("Cumulative NDMI 10 days before burn ($S_{NDMI,10}$)"),
       y = "Burn Severity (here dNBR)") +
  theme_cowplot(FONT_SIZE) +
  theme(legend.position = "none")

# 6. Align subplots ----
layout <- "
#AAA#
BBBCC
DDDDD
"

pg <- fig_1a + fig_1b + fig_1c + fig_1d +
  plot_layout(design = layout)

ggsave2(pg, "figures/Figure_1.png",bg = "white",width = 16, height = 14)

