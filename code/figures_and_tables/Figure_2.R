library(cowplot)
library(viridis)
library(colorspace)
library(tidyverse)
library(scales)
library(latex2exp)
library(terra)
library(geostats)
library(tidyterra)
library(grid)
library(patchwork)
library(readr)
library(pbapply)
library(data.table)
set.seed(10)

# 1. Load data and configure stuff ----
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

FONT_SIZE <- 20

# Output directory for sample tables
OUT_DIR <- paste0(TABLE_DIR,"sampled_data/")

# Load lookup tables
processing_lut <- read.csv(paste0(TABLE_DIR,"processing_LUT.csv")) %>%  # overall LUT
  filter(tst_year >= 2017) 

# Load burn severity raster
optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)

# Load features (fire perimeters and ROIs)
fire_perimeters <- vect(
  paste0(DATA_DIR,"feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg")
)

# Build global model data table
topN_fires <- fire_perimeters %>%
  filter(tst_year >= 2017) %>% 
  arrange(desc(farea)) %>% 
  slice_head(n = 25) 

TEST_ID <- topN_fires$fireid

if (length(TEST_ID)>0){
  subset_lut <- filter(processing_lut, fireid %in% TEST_ID)
}

SAMPLE_ID <- c(14211,
               14553,
               20922)

PlotSeverityMapHist <- function(cropped_severity_raster,fire_perimeter,
                                ONLY_MAP = FALSE, FONT_SIZE = 18){
  burn_severity_index <- names(cropped_severity_raster)
  
  binwidth <- .02
  
  # Plot dNBR map
  p1 <- ggplot() +
    geom_spatraster(data = cropped_severity_raster) +
    geom_spatvector(data = fire_perimeter, color = "white",fill = NA,
                    size = 2) +
    scale_fill_viridis_c(option = "inferno",
                         na.value = "white",
                         name = burn_severity_index) +
    theme_cowplot(FONT_SIZE)
  
  if (ONLY_MAP){
    return(p1)
  }
  
  # Plot histogram for this map
  xvar <- as.symbol(burn_severity_index)
  p2 <- ggplot(data = cropped_severity_raster, aes(x = !!xvar)) +
    geom_histogram(binwidth = binwidth, color = "black", fill = "indianred1") +
    xlim(c(-.5,1)) +
    labs(
      title = sprintf("Histogram of %s",burn_severity_index),
      x = sprintf("%s Values",burn_severity_index),
      y = "Frequency"
    ) +
    theme_cowplot(FONT_SIZE)
  
  pg <- p1 + p2
  return(pg)
}

# 2. Extract burn severity densities in fire perimeters ----
extract_burn_severity <- function(fire_perimeter, optimality_lut,
                                  FIRE_ID, burn_severity_index,
                                  PLOT_MAP_HISTOGRAM = FALSE,ONLY_MAP = FALSE,
                                  SAVE_DF = FALSE){
  
    fname_optimal_severity_raster <- optimality_lut %>% 
      filter(fireid == FIRE_ID,
             severity_index == "dNBR") %>% 
      pull(fname_severity_raster)
    
    if (OS == "Windows"){
      fname_optimal_severity_raster <- gsub("/home/nrietz/scratch","data",
                                            fname_optimal_severity_raster)
    }
      
    rast_burn_severity <- rast(gsub("dNBR",burn_severity_index,
                                    fname_optimal_severity_raster))
    
    selected_fire_perimeter <- fire_perimeters %>% 
      filter(fireid  == FIRE_ID) %>%
      project(crs(rast_burn_severity))
    
    # Crop raster
    dnbr_in_perimeter <- rast_burn_severity %>% 
      mask(selected_fire_perimeter, updatevalue = NA) %>%
      crop(ext(selected_fire_perimeter))
    
    if (PLOT_MAP_HISTOGRAM){
      p <- PlotSeverityMapHist(dnbr_in_perimeter, selected_fire_perimeter,
                               ONLY_MAP = ONLY_MAP)
      ggsave2(p,filename = sprintf("figures/burn_severity_maps/%s_%s_map_histogram.png",
                         FIRE_ID, burn_severity_index),
              bg = "white", width = 10,height = 8)
      
      return(p)
    }
    
    df_out <- dnbr_in_perimeter %>% 
      as.data.frame() %>% 
      mutate(FIRE_ID = FIRE_ID)
    
    if (SAVE_DF){
      write_csv2(df_out,paste0(OUT_DIR,
                               sprintf("%s_%s_in_perimeter.csv",
                                       FIRE_ID, burn_severity_index)) )
    }
    
    return(df_out)
  
}

burn_severity_index <- "dNBR_corr"

burn_severity_list <- pblapply(subset_lut$fireid, extract_burn_severity,
         fire_perimeter = fire_perimeter, optimality_lut = optimality_lut,
         burn_severity_index = burn_severity_index, 
         PLOT_MAP_HISTOGRAM = FALSE,
         SAVE_DF = FALSE)

df_burn_severity <- rbindlist(burn_severity_list)

df_burn_severity <- df_burn_severity %>% 
  mutate(VIZ_FIRES = as.factor(ifelse(FIRE_ID %in% SAMPLE_ID, FIRE_ID, NA)),
         FIRE_ID = as.factor(FIRE_ID))

# 3. Plot density curves ----
xvar <- as.symbol(burn_severity_index)
sample_colors <- c("#E16A86", "#50A315", "#009ADE")

(p_dens <- ggplot(data = df_burn_severity, aes(x = !!xvar, group = FIRE_ID)) +
    geom_line(stat = "density", linewidth = 1, alpha = 0.2, color = "gray70") +
    geom_line(
      data = filter(df_burn_severity, !is.na(VIZ_FIRES)),
      aes(x = !!xvar, group = VIZ_FIRES, color = VIZ_FIRES),
      stat = "density", linewidth = 2, alpha = 0.8
    ) +
    labs(
      x = "Burn severity (dNBR)",
      y = "Density",
      color = "",
      subtitle = "d)"
    ) +
    scale_color_manual(values = sample_colors, 
                       labels = c("lowland tundra (2020)",
                                  "upland tundra (2020)",
                                  "lowland tundra (2019)"),
                       guide = guide_legend(ncol = 1)) +
    theme_cowplot(FONT_SIZE) +
    theme(
      legend.position = c(0.75, 0.95),
      legend.justification = c(0.5, 1),
      legend.background = element_blank(), 
      legend.key = element_blank(),
      legend.title = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.subtitle=element_text(size=FONT_SIZE,face="bold")
    )
)

# 4. Plot maps ----
plot_burn_severity_map <- function(fire_perimeter, optimality_lut,
                                  FIRE_ID, burn_severity_index,
                                  FONT_SIZE = 18){
  
  fname_optimal_severity_raster <- optimality_lut %>% 
    filter(fireid == FIRE_ID,
           severity_index == "dNBR") %>% 
    pull(fname_severity_raster)
  
  if (OS == "Windows"){
    fname_optimal_severity_raster <- gsub("/home/nrietz/scratch","data",
                                          fname_optimal_severity_raster)
  }
  
  rast_burn_severity <- rast(gsub("dNBR",burn_severity_index,
                                  fname_optimal_severity_raster))
  
  selected_fire_perimeter <- fire_perimeters %>% 
    filter(fireid  == FIRE_ID) %>%
    project(crs(rast_burn_severity))
  
  # Crop raster
  cropped_severity_raster <- rast_burn_severity %>% 
    crop(ext(selected_fire_perimeter))
  
  # create inverse mask for plotting
  extent_rect <- as.polygons(ext(cropped_severity_raster))
  crs(extent_rect) <- crs(cropped_severity_raster)
  plot_mask <- erase(extent_rect, selected_fire_perimeter)
  
  p <- ggplot() +
    geom_spatraster(data = cropped_severity_raster, show.legend = TRUE) +
    geom_spatvector(data = selected_fire_perimeter, color = "white",fill = NA,
                    size = 2) +
    geom_spatvector(data = plot_mask, color = NA,fill = "white",alpha = 0.8) +
    scale_fill_viridis_c(option = "inferno",
                         na.value = "white",
                         name = "Burn severity (dNBR)",
                         limits = c(-0.5, 1.1)) +
    theme_map(FONT_SIZE)
  
  # create scale bar
  plot_limits <- ggplot_build(p)$layout$panel_params[[1]]
  xlim <- plot_limits$x_range
  ylim <- plot_limits$y_range
  
  # Calculate position for upper right corner
  xtext <- xlim[1] + 0.2 * diff(xlim)  # 15% from left edge
  ytext <- ylim[1] + 0.1 * diff(ylim)   # 10% from bottom edge
  
  # Length of scale bar (in map units; 2000 = 2 km if 1 unit = 1 meter)
  scale_length <- 2000
  
  p <- p +
    # Scale bar (rectangle at calculated position)
    geom_rect(
      aes(xmin = xtext - scale_length, xmax = xtext,
          ymin = ytext, ymax = ytext + diff(ylim)*0.01), # height of bar
      fill = 'gray20'
    ) +
    # Label ('2 km') above the bar
    geom_text(
      aes(x = xtext - scale_length/2, y = ytext + diff(ylim)*0.05, label = '2 km'),
      size = 5, fontface = 'bold', color = 'gray20'
    ) +
    # North arrow (arrow shape using text or geom_segment)
    geom_text(
      aes(x = xtext + diff(xlim)*0.1, y = ytext + diff(ylim)*0.02, label = 'N'),
      size = 7, fontface = 'bold', color = 'gray20'
    ) +
    geom_text(
      aes(x = xtext + diff(xlim)*0.1, y = ytext + diff(ylim)*0.08, label = '↑'),
      size = 10, color = 'gray20'
    )
  
  return(p)

}

p1 <- plot_burn_severity_map(fire_perimeter, optimality_lut,
                              SAMPLE_ID[1], burn_severity_index) +
    theme(panel.border = element_rect(colour = sample_colors[1], 
                                      fill=NA, linewidth=2),
          plot.subtitle=element_text(size=FONT_SIZE,face="bold")) +
  labs(subtitle = "a)")

p2 <- plot_burn_severity_map(fire_perimeter, optimality_lut,
                              SAMPLE_ID[2], burn_severity_index) +
    theme(panel.border = element_rect(colour = sample_colors[2], 
                                      fill=NA, linewidth=2),
          plot.subtitle=element_text(size=FONT_SIZE,face="bold")) +
  labs(subtitle = "b)")

p3 <- plot_burn_severity_map(fire_perimeter, optimality_lut,
                              SAMPLE_ID[3], burn_severity_index) +
    theme(panel.border = element_rect(colour = sample_colors[3], 
                                      fill=NA, linewidth=2),
          plot.subtitle=element_text(size=FONT_SIZE, face="bold")) +
  labs(subtitle = "c)")

# 5. Create Figure 2 ----
pg <- p1 + p2 + p3 + p_dens + 
  # plot_annotation(tag_levels = 'a',tag_suffix = ")") +
  plot_layout(guides = 'collect') &
  theme(plot.tag.position = c(-0.03, 1),
        plot.tag = element_text(size = FONT_SIZE, hjust = 0, vjust = 0),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 20))

ggsave2(pg,filename = "figures/manuscript/Figure_2.png",
        bg = "white", width = 13,height = 10)


# 6. Check semivariogram & dNBR spatial heterogeneity ----
FIRE_ID <- SAMPLE_ID[1]

fname_optimal_severity_raster <- optimality_lut %>% 
  filter(fireid == FIRE_ID,
         severity_index == "dNBR") %>% 
  pull(fname_severity_raster)

if (OS == "Windows"){
  fname_optimal_severity_raster <- gsub("/home/nrietz/scratch","data",
                                        fname_optimal_severity_raster)
}

rast_burn_severity <- rast(gsub("dNBR",burn_severity_index,
                                fname_optimal_severity_raster))

selected_fire_perimeter <- fire_perimeters %>% 
  filter(fireid  == FIRE_ID) %>%
  project(crs(rast_burn_severity))

# Crop raster
cropped_severity_raster <- rast_burn_severity %>% 
  crop(ext(selected_fire_perimeter))

## a. fit semivariogram ----
points_df <- as.data.frame(mask(cropped_severity_raster, selected_fire_perimeter,
                                updatevalue = NA),
                           xy = TRUE, na.rm = TRUE)
colnames(points_df)[3] <- "value"

sampled_df <- sample_frac(points_df,.1)
sv <- semivariogram(x = sampled_df$x,y = sampled_df$y,
                    z = sampled_df$value,
                    bw = 30, nb = 100)

# plot semivariogram
df_sv <- data.frame(distance = sv$h, gamma = sv$sv)

ggplot(df_sv, aes(x = distance,y = gamma)) +
  geom_point() + 
  geom_vline(aes(xintercept = sv$snr[3]), lty = "dashed", color = "gray40") +
  annotate("text", x = sv$snr[3], y = max(df_sv$gamma),
           label = sprintf("Effective range = %.0f m",sv$snr[3]),
           size = 5, color = "grey40") +
  labs(x = "Distance (m)", y = "Semivariance (spherical)") +
  theme_cowplot()

ggsave2(filename = sprintf("figures/dNBR_corr_semivariogram_%s.png",FIRE_ID),
        width = 10,height = 10, bg = "white")


## b. Calculate focal statistics ----
mu_dnbr_in_perimeter <- focal(cropped_severity_raster,w = 33,fun = "mean",na.policy = "all",na.rm = T)
names(mu_dnbr_in_perimeter) <- "dNBR_mean"
sd_dnbr_in_perimeter <- focal(cropped_severity_raster,w = 33,fun = "sd",na.policy = "all",na.rm = T)
names(sd_dnbr_in_perimeter) <- "dNBR_sd"
cv_dnbr_in_perimeter <- sd_dnbr_in_perimeter  / mu_dnbr_in_perimeter
names(cv_dnbr_in_perimeter) <- "dNBR_cv"

df_cv <- c(sd_dnbr_in_perimeter,
           mu_dnbr_in_perimeter,
           cv_dnbr_in_perimeter) %>% 
  mask(selected_fire_perimeter, updatevalue = NA) %>% 
  as.data.frame()

# Plot std.dev. vs. mean dNBR
p <- df_cv %>% drop_na() %>% sample_frac(.2) %>% 
  ggplot(aes(x = dNBR_mean,y = dNBR_sd)) + 
  geom_point(alpha = 0.05) + 
  labs(x = "Focal mean dNBR", y = "std. dev. dNBR") + 
  theme_cowplot(18)

p2 <- ggExtra::ggMarginal(p, type="density")

ggsave2(p2, filename = sprintf("figures/dNBR_corr_sd_vs_mean_focal_33_%s.png",FIRE_ID),
        width = 10,height = 10, bg = "white")

# Plot CV vs.dNBR
df_cv %>% sample_frac(.2) %>% 
  ggplot(aes(x = dNBR_mean,y = dNBR_cv)) + 
  geom_point(alpha = 0.1) + 
  labs(x = "Focal mean dNBR", y = "CV dNBR") + 
  lims(y = c(-5,5)) +
  theme_cowplot(18)

# 7. Burn severity vs. burned area ----
df1 <- df_burn_severity %>% left_join(fire_perimeters %>% 
                                        as.data.frame() %>% 
                                        select(fireid, farea) %>% 
                                        mutate(FIRE_ID = as.factor(fireid)),
                                      by = "FIRE_ID")
df2 <- df1 %>% summarize(avg_dnbr = mean(!!sym(burn_severity_index)),
                         sd_dnbr = sd(!!sym(burn_severity_index)),
                         cv_dnbr = sd(!!sym(burn_severity_index)) /mean(!!sym(burn_severity_index)),
                         farea = first(farea),
                         .by = FIRE_ID)

ggplot(df2, aes(x = farea,y = avg_dnbr)) +
  geom_smooth(method = "lm",color = "black") +
  geom_point(size = 5, alpha = .7) + 
  scale_x_log10() +
  labs(x = "Burned area (ha)", y = "Average burn severity (dNBR)") +
  theme_cowplot(18)

ggsave(sprintf("figures/burned_area_vs_mean_%s_relationship.png",burn_severity_index),
       width = 8, height = 8, bg = "white")

mod1 <- lm(avg_dnbr ~ log10(farea), data = df2)
summary(mod1)


# 8. Check variation of dNBR distributions of one UTM tile ----
flist <- list.files(path = paste0(HLS_DIR,"/severity_rasters/"),
                    pattern = "^dNBR_corr_54WXE_2019.*tif$", full.names = T)

dnbr_rasters <- rast(flist) 
names(dnbr_rasters) <- lapply(sources(dnbr_rasters), 
                              function(x) str_split(basename(tools::file_path_sans_ext(x)),"_",simplify = T)[4])

perimeter_in_54wxe <- fire_perimeters %>% 
  filter(fireid == 21286) %>% 
  project(crs(dnbr_rasters))
  
cropped_dnbr_rasters <- dnbr_rasters %>% 
  crop(perimeter_in_54wxe) %>% 
  mask(perimeter_in_54wxe,updatevalue = NA)

df_dnbr_rasters <- as.data.frame(cropped_dnbr_rasters) %>% 
  drop_na() %>% 
  as_tibble() %>% 
  pivot_longer(cols = tidyr::everything(), values_to = "dNBR", names_to = "date") %>% 
  mutate(date = ymd(date),
         date_fct = as.factor(date)) %>% 
  filter(date >= ymd(paste(c(perimeter_in_54wxe$ted_year,
                             perimeter_in_54wxe$ted_month,
                             perimeter_in_54wxe$ted_day),collapse = "-")))

ggplot(df_dnbr_rasters, aes(x = dNBR, group = date_fct, color = date_fct)) +
  geom_density() + 
  scale_color_scico_d(palette = "lajolla") +
  theme_cowplot(FONT_SIZE)
