library(tidyverse)
library(cowplot)
library(patchwork)
library(colorspace)
library(gt)
set.seed(10)

# 1. Load & format lmer results ----
burn_severity_index <- "dnbr_corr"

fname_results_rds <- sprintf("data/models/results_lmer_%s_nocumsum_fireid_rastdoy.rds",burn_severity_index)

results <- readRDS(fname_results_rds)

results_clean <- results %>%
  mutate(
    # Extract base predictor: NDVI, NDMI, LST (or keep as-is for Intercept, etc.)
    base_predictor = case_when(
      str_detect(predictor, "^NDVI\\.d_prefire_") ~ "NDVI",
      str_detect(predictor, "^NDMI\\.d_prefire_") ~ "NDMI",
      str_detect(predictor, "^lst_pred_d_prefire_") ~ "LST",
      TRUE ~ predictor  # Keep Intercept or others as-is
    ),
    
    # Extract numeric day from "d_prefire_" part
    day_from_name = as.integer(str_extract(predictor, "(?<=d_prefire_)[0-9]+")),
    
    # Replace with parsed day if available
    days_before_fire = coalesce(day_from_name, days_before_fire)
  ) %>%
  select(-day_from_name) 

FONT_SIZE <- 18

# 2. Plot model results ----

## a. R2 curves ----
max_r2_day <- results_clean %>%
  filter(r_squared_marg == max(r_squared_marg, na.rm = TRUE)) %>%
  select(days_before_fire, r_squared_marg) %>% 
  distinct()

max_day <- max_r2_day$days_before_fire
max_r2 <- max_r2_day$r_squared_marg

s <- sprintf("Maximum marginal R² is \n %s days before the fire", max_day)

ylims <- c(0, 0.25)

(p_r2_marg <- ggplot(results) + 
    geom_point(aes(x = days_before_fire, y = r_squared_marg)) +
    geom_vline(aes(xintercept = max_day), lty = "dashed") +
    annotate("text", x = 20, y = .05,
             label = s,
             size = 5, color = "grey40") +
    geom_curve(aes(x = 17, y = .05, 
                   xend = max_day, yend = max_r2 -0.005),
               arrow = arrow(length = unit(0.08, "inch")), linewidth = 1,
               color = "grey40", curvature = 0.3) +
    scale_x_reverse() +
    scale_y_continuous(breaks = seq(ylims[1],ylims[2],.05),
                       labels = seq(ylims[1],ylims[2],.05),
                       limits = ylims
    ) +
    labs(x = "Days before fire", y = expression("marginal"~R^2)) + 
    theme_cowplot(FONT_SIZE) )

max_r2_cond_day <- results_clean %>%
  filter(r_squared_cond == max(r_squared_cond, na.rm = TRUE)) %>%
  select(days_before_fire, r_squared_cond) %>% 
  distinct()

max_cond_day <- max_r2_cond_day$days_before_fire
max_cond_r2 <- max_r2_cond_day$r_squared_cond

s_cond <- sprintf("Maximum conditional R² is \n %s days before the fire", max_cond_day)

ylims <- c(0.6, 0.8)

(p_r2_cond <- ggplot(results) + 
    geom_point(aes(x = days_before_fire, y = r_squared_cond)) +
    geom_vline(aes(xintercept = max_cond_day), lty = "dashed") +
    annotate("text", x = max_cond_day + 8, y = max_cond_r2 + .05,
             label = s_cond,
             size = 5, color = "grey40") +
    geom_curve(aes(x = max_cond_day + 8, y = max_cond_r2 + .035, 
                   xend = max_cond_day, yend = max_cond_r2 + 0.005),
               arrow = arrow(length = unit(0.08, "inch")), linewidth = 1,
               color = "grey40", curvature = 0.3) +
    scale_x_reverse() +
    labs(x = "Days before fire", y = "conditional R²") + 
    scale_y_continuous(breaks = seq(ylims[1],ylims[2],.05),
                       labels = seq(ylims[1],ylims[2],.05),
                       limits = ylims
    ) +
    theme_cowplot(FONT_SIZE) )

if (SAVE_FIGURES){
  ggsave2(p_r2_cond,filename = sprintf("figures/Figure_S_R2cond_%s.png",burn_severity_index),
          width = 8,height = 8,bg = "white")
}

## b.  Plot effect sizes over time ----
mod_labs <- c(
  burn_doy_scaled = "Date of burn \n(Day of year)",
  NDVI = "Vegetation greenness \n(NDVI HLS)", 
  NDMI = "Vegetation moisture \n(NDMI HLS)",
  LST = "Land surface temperature \n(LST Landsat-8)",
  raster_doy_scaled = "Season of burn severity raster \n(Day of year)",
  elevation = "Elevation",
  slope = "Terrain slope",
  northness = "Northness",
  eastness = "Eastness"
)

p_eff <- (results_clean %>%
            mutate(
              base_predictor = factor(base_predictor, levels = names(mod_labs)),
              significant = p_val < 0.05,
              alpha_val = ifelse(significant, 1, 0.5)
            ) %>%
            filter(base_predictor %in% c("NDVI", "NDMI", "LST",
                                         "elevation", "slope", "northness", "eastness", 
                                         "burn_doy_scaled")) %>%
            ggplot(aes(
              x = days_before_fire,
              y = coef,
              color = base_predictor,
              alpha = alpha_val,
              group = base_predictor
            )) +
            geom_point() +
            geom_errorbar(aes(
              ymin = coef - std_err,
              ymax = coef + std_err
            ), width = 0.5) +
            geom_hline(yintercept = 0, lty = "dashed") +
            scale_x_reverse() +
            facet_wrap(~ base_predictor,
                       scales = "free_y", 
                       labeller = as_labeller(mod_labs),
                       nrow = 2) +
            labs(x = "Days before fire", y = "Effect size (coefficient)") +
            scale_alpha_identity() +
            theme_cowplot(FONT_SIZE) +
            theme(legend.position = "none"))

p_eff / p_r2_marg +
  plot_annotation(tag_levels = 'a',tag_suffix = ")") +
  plot_layout(heights = c(3, 2))

ggsave2(filename = sprintf("figures/manuscript/Figure_3_%s.png",burn_severity_index),
        width = 14,height = 8,bg = "white")


# 3. Table SYY ----
burn_severity_indices <- c("dnbr_corr","rdnbr_corr","rbr","dgemi")

extract_r2 <- function(burn_severity_index){
  fname_results_rds <- sprintf("data/models/results_lmer_%s_nocumsum_fireid_rastdoy.rds",burn_severity_index)

  results <- readRDS(fname_results_rds)
  
  results_clean <- results %>%
    mutate(
      # Extract base predictor: NDVI, NDMI, LST (or keep as-is for Intercept, etc.)
      base_predictor = case_when(
        str_detect(predictor, "^NDVI\\.d_prefire_") ~ "NDVI",
        str_detect(predictor, "^NDMI\\.d_prefire_") ~ "NDMI",
        str_detect(predictor, "^lst_pred_d_prefire_") ~ "LST",
        TRUE ~ predictor  # Keep Intercept or others as-is
      ),
      
      # Extract numeric day from "d_prefire_" part
      day_from_name = as.integer(str_extract(predictor, "(?<=d_prefire_)[0-9]+")),
      
      # Replace with parsed day if available
      days_before_fire = coalesce(day_from_name, days_before_fire)
    ) %>%
    select(-day_from_name) 

  max_r2_day <- results_clean %>%
    filter(r_squared_marg == max(r_squared_marg, na.rm = TRUE)) %>%
    select(days_before_fire, r_squared_marg) %>% 
    distinct() %>% 
    rename(days_before_fire_marg = days_before_fire)
  
  max_r2_cond_day <- results_clean %>%
    filter(r_squared_cond == max(r_squared_cond, na.rm = TRUE)) %>%
    select(days_before_fire, r_squared_cond) %>% 
    distinct() %>% 
    rename(days_before_fire_cond = days_before_fire)
  
  combined <- bind_cols(max_r2_day, max_r2_cond_day) %>%
    add_column(burn_severity_index = burn_severity_index,.before = 1)
  
  return(combined)
  
}

df <- map_dfr(burn_severity_indices, extract_r2)

(tab_SYY <- df %>% 
    mutate(burn_severity_index = recode(
      burn_severity_index,
      "dnbr_corr" = "dNBR",
      "rdnbr_corr" = "RdNBR",
      "rbr" = "RBR",
      "dgemi" = "dGEMI"
    )) %>% 
    gt() %>% 
    fmt_number(
      n_sigfig = 2
    ) %>% 
    cols_label(
      burn_severity_index = md("**Burn Severity Index**"),
      days_before_fire_marg = md("**Day of maximum R^2^<sub>marginal</sub>**"),
      days_before_fire_cond = md("**Day of maximum R^2^<sub>conditional</sub>**"),
      r_squared_marg = md("**Maximum R^2^<sub>marginal</sub>**"),
      r_squared_cond = md("**Maximum R^2^<sub>conditional</sub>**")
    ))

gtsave(tab_SYY,filename = "data/tables/Table_SYY.html")
