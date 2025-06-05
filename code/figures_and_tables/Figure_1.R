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
burn_severity_index <- "dNBR_corr"
FIRE_ID <- 14205

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

# Load real data from fire scar
data_all <- read_csv2(
  paste0(DATA_DIR,
         sprintf("tables/model_dataframes/1pct/%s_model_dataframe.csv",FIRE_ID))
  )

# Load burn severity raster
optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

fname_optimal_severity_raster <- optimality_lut %>% 
  filter(fireid == FIRE_ID,
         severity_index == "dNBR") %>% 
  pull(fname_severity_raster)

rast_burn_severity <- rast(gsub("dNBR",burn_severity_index,fname_optimal_severity_raster))

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  paste0(DATA_DIR,"feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")
)

selected_fire_perimeter <- fire_perimeters %>% 
  filter(fireid  == FIRE_ID) %>%
  project(crs(rast_burn_severity))

dnbr_in_perimeter <- rast_burn_severity %>% 
  mask(ext(selected_fire_perimeter), updatevalue = NA) %>% 
  crop(ext(selected_fire_perimeter))

# Load sample points (with associated burn dates)
fname_sample_points <- paste0(DATA_DIR,
                              sprintf("feature_layers/1pct/%s_sample_points_1pct_burn_date.gpkg",
                               FIRE_ID))

sample_points <- vect(fname_sample_points) %>% 
  project(crs(rast_burn_severity)) %>% 
  mutate(ObservationID = 1:nrow(.))

# Define x-axis: Days since fire (from -50 to +30)
days_since_burn <- seq(-50, 30, by = 1)  # Pre-fire (-50 to 0), Post-fire (0 to 30)

# 2. Construct Figure 1 ----
FONT_SIZE <- 18

## a. Figure 1a) - Overview map ----
fig_1a <- ggdraw() + 
  draw_image('figures/Figure_1a.png', scale = 1) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

## b. Figure 1b) - burned area over time ----
(fig_1b <- ggplot(topN_fires,aes(x = tst_year, y = farea *100)) +
    geom_bar(stat = 'identity', fill = "#730000") +
    labs(x = "Fire Year",
         y = "Total annual burned area of studied fires (ha)") +
    theme_cowplot(FONT_SIZE) +
    theme(legend.position = "bottom"))

if (SAVE_FIGURES){
  ggsave2(fig_1b,filename = "figures/Figure_1b.png", bg ="white",
          width = 6, height = 7)
}  

## c. Plot Figure 1c) - Plot map of fire scar ----
# point_index <- c(78,74,127)
point_index <- df_allfires %>%
  filter(fireid == FIRE_ID) %>%
  distinct(ObservationID) %>% sample_n(3) %>% 
  pull()

plot_points <- sample_points %>% 
  filter(ObservationID %in% point_index)

dx <- (ext(dnbr_in_perimeter)[2] - ext(dnbr_in_perimeter)[1])
dy <- (ext(dnbr_in_perimeter)[4] - ext(dnbr_in_perimeter)[3])
xtext <- ext(dnbr_in_perimeter)[1] + dx *.2
ytext <- ext(dnbr_in_perimeter)[4] - dy *.3

(fig_1c <- ggplot() +
    geom_spatraster(data = dnbr_in_perimeter) +
    geom_spatvector(data = plot_points,
                    aes(color = as.factor(point_index)),
                    size = 7,show.legend = FALSE) +
    geom_spatvector(data = selected_fire_perimeter, fill = NA,
                    color = "black",linewidth = 1,show.legend = FALSE) +
    # annotate perimeter
    annotate("segment",
             x = xtext - 500, xend = 630000, y = ytext - 1300, yend = 7920000,
             colour = "black",linewidth = 1) +
    geom_text(aes(x = xtext - 800, y = ytext,
                  label = 'Fire \nperimeter',fontface = 'bold'),
              size = 5,
              colour = 'black') +
    scale_fill_continuous_diverging(palette = "Broc", na.value = "white",
                                    rev = FALSE,name = "Burn Severity (dNBR)",
                                    guide = guide_colourbar(
                                      title.position = "top",
                                      title.hjust = 0.5,
                                      barwidth = unit(20, "lines"))) +
    scale_color_viridis_d(end = 0.8, option = "inferno") +
    # labs(title = "a) Sampling data in fire perimeter.") +
    theme_map(FONT_SIZE) +
    theme(legend.position = "bottom",
          legend.justification = "center")
)

if (SAVE_FIGURES){
  ggsave2(fig_1c,filename = "figures/Figure_1c.png",
          width = 8, height = 6,bg = "white")
}

## d. Plot Figure 1d) - NDMI curves ----
point_index <- c(335 ,520 ,301)

df_transformed <- df_allfires %>% 
  filter(fireid == FIRE_ID, ObservationID %in% point_index) %>% 
  select(-matches("NDVI|LST|cum")) %>% 
  pivot_longer(
    cols = starts_with("NDMI.d_prefire_"),
    names_to = "d_prefire_name",
    values_to = "NDMI_fit"
  ) %>%
  mutate(d_prefire = as.integer(str_extract(d_prefire_name, "\\d+"))) %>%
  filter(d_prefire > 0 & d_prefire <= 30)

df_filtered_obs <- df_allfires %>% 
  filter(fireid == FIRE_ID, ObservationID %in% point_index) %>%
  mutate(d_prefire = burn_doy - doy,
         ObservationID = factor(ObservationID, levels = point_index),) %>%
  filter(d_prefire > 0 & d_prefire <= 30)

facet_labels <- setNames(paste(c("Low","Moderate","High"),"\nmoisture levels"),point_index)

(fig_1d <- ggplot() +
  geom_point(data = df_filtered_obs, aes(x = d_prefire, y = DailyMeanNDMI,
                                         color = as.factor(ObservationID)),
             show.legend = FALSE) +
  geom_line(data = df_transformed, aes(x = d_prefire, y = NDMI_fit), 
            color = "blue") +
  geom_vline(xintercept = 0,lty = "dashed") +
  facet_wrap(~ ObservationID, labeller = as_labeller(facet_labels)) +
  scale_color_viridis_d(end = 0.8, option = "inferno") +
  labs(
    x = "Days Before Fire", 
    y = "NDMI" ) +
  scale_x_reverse() +
  theme_cowplot(FONT_SIZE))

if (SAVE_FIGURES){
  ggsave2(fig_1d,filename = "figures/Figure_1d.png",
          width = 8, height = 4,bg = "white")
}

## e. Plot Figure 1e) - scatterplot for 10d pre fire ----
(fig_1e <- ggplot() +
   geom_point(data = df_allfires %>% filter(fireid == FIRE_ID),
              aes(x = NDMI.d_prefire_12, y = dnbr_corr),
              color = "gray80", alpha = .2) +
   geom_point(aes(x = c(-.15, 0, .15), y = c(.3, .15, -.15),
                  color = as.factor(sort(point_index,decreasing = T))), size = 10) +
   scale_color_viridis_d(end = 0.8, option = "inferno",direction = -1) +
   labs(
     x = "NDMI 12 days before fire",
     y = "Burn Severity (dNBR)"
   ) +
   theme_cowplot(FONT_SIZE) +
   theme(legend.position = "none"))

if (SAVE_FIGURES){
  ggsave2(fig_1e,filename = "figures/Figure_1e.png",
          width = 14, height = 6,bg = "white")
}

# 5. Align subplots ----
layout <- "
AABB
CCDD
EE##
"

pg <- fig_1a + fig_1b + fig_1c + fig_1d + fig_1e +
  plot_layout(design = layout)

ggsave2(pg, "figures/Figure_1.png",bg = "white",width = 14, height = 18)


# old plot code ----
fig_1d <- ggplot(ndmi_data, aes(x = DaysSinceBurn, y = NDMI, color = Curve)) +
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

# 4. Plot Figure 1c) - cumulative sum of NDMI curves
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
