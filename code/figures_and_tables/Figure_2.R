library(cowplot)
library(viridis)
library(colorspace)
library(tidyverse)
library(scales)
library(latex2exp)
library(terra)
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
  
    rast_burn_severity <- rast(gsub("dNBR",burn_severity_index,
                                    fname_optimal_severity_raster))
    
    selected_fire_perimeter <- fire_perimeters %>% 
      filter(fireid  == FIRE_ID) %>%
      project(crs(rast_burn_severity))
    
    # Crop raster
    dnbr_in_perimeter <- rast_burn_severity %>% 
      # mask(selected_fire_perimeter, updatevalue = NA) %>% 
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
      color = ""
    ) +
    scale_color_manual(values = sample_colors, 
                       labels = c("lowland tundra (2020)",
                                  "upland tundra (2020)",
                                  "lowland tundra (2019)"),
                       guide = guide_legend(ncol = 1)) +
    theme_cowplot() +
    theme(
      legend.position = c(0.75, 0.95),
      legend.justification = c(0.5, 1),
      legend.background = element_blank(), 
      legend.key = element_blank(),
      legend.title = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank()
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
    geom_spatraster(data = cropped_severity_raster, show.legend = FALSE) +
    geom_spatvector(data = selected_fire_perimeter, color = "white",fill = NA,
                    size = 2) +
    geom_spatvector(data = plot_mask, color = NA,fill = "white",alpha = 0.8) +
    scale_fill_viridis_c(option = "inferno",
                         na.value = "white",
                         name = burn_severity_index,
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

(p1 <- plot_burn_severity_map(fire_perimeter, optimality_lut,
                              SAMPLE_ID[1], burn_severity_index) +
    theme(panel.border = element_rect(colour = sample_colors[1], 
                                      fill=NA, linewidth=2)) )

(p2 <- plot_burn_severity_map(fire_perimeter, optimality_lut,
                              SAMPLE_ID[2], burn_severity_index) +
    theme(panel.border = element_rect(colour = sample_colors[2], 
                                      fill=NA, linewidth=2)) )

(p3 <- plot_burn_severity_map(fire_perimeter, optimality_lut,
                              SAMPLE_ID[3], burn_severity_index) +
    theme(panel.border = element_rect(colour = sample_colors[3], 
                                      fill=NA, linewidth=2)) )
  

# 5. Create Figure 2 ----
pg <- p1 + p2 + p3 + p_dens

ggsave2(pg,filename = "figures/Figure_2.png",
        bg = "white", width = 10,height = 10)


