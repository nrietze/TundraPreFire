library(terra)
library(tidyterra)
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)
library(broom)
library(tictoc)
library(future)
library(pbapply)
library(gt)
set.seed(10)

# 1. Load models ----
model <- readr::read_rds("data/models/brm_benchmarking_20pct_day1.rds")

mod_files <- list.files("data/models/",pattern = "*day1.rds",full.names = T)

all_models <- lapply(mod_files,read_rds)

# 3. Model assessment ----
# Combine dataframes from different sampling runs and add identifier column "run"
all_data_frames <- do.call(rbind, lapply(all_models, function(x) x$data %>% mutate(run = nrow(x$data))))

all_data_frames <- rbind(all_data_frames, df_subset %>% 
                          select(colnames(all_data_frames)[1:10]) %>% 
                          mutate(run = nrow(df_subset)))

ggplot(all_data_frames) +
  geom_histogram(aes(x=dnbr), color = "tomato",fill = "salmon") +
  labs(x = "Burn severity (dNBR)",
     title = "Distribution of dNBR for benchmark subset") +
  facet_wrap(~run)+
  theme_cowplot()

ggsave("figures/hist_dNBR_percentage_benchmark.png",bg = "white",
       width = 8, height = 8)

  
## a. Report posterior checks ----
model_summaries <- lapply(all_models, function(model) {
  list(
    ess = bayestestR::effective_sample(model),
    r2  = performance::r2_bayes(model)
  )
})

# View effective sample sizes from model 1
model_summaries[[1]]$ess %>% gt()

# View R² for model 1
r2s <- unlist(lapply(model_summaries, function(x) {x$r2$R2_Bayes_marginal}))

ggplot() +
  geom_point(aes(x = factor(c("1","20","30","40","50")),y = r2s),
             size = 5) +
  labs(x = "Percentage sampled from total data frame", 
       y = "marginal pseudo R-squared",
       title = "R-squared for Day 1 before fire") +
  theme_cowplot(18)

ggsave("figures/R2_day1_benchmark.png",bg = "white",
       width = 8, height = 8)


model <- all_models[[1]]

# Report effective sample size
bayestestR::effective_sample(model) %>% gt()

# Get pseudo R2
performance::r2_bayes(model)

# model sumamry
summary(model)

# plots
plot(model)

pp_check(model)

pp_check(model, type='error_scatter_avg')
pp_check(model, type='error_hist')

mod_labs <-c(elevation = 'Elevation',
             slope = 'Slope',
             northness = 'Northness',
             eastness = 'Eastness',
             doy = "Burn timing (Day of Year)",
             cumsum_lst_d_prefire_10 = "Land surface temperature \n(LST Landsat)",
             NDVI.d_prefire_10 = "Greenness (NDVI HLS)", 
             NDMI.d_prefire_10 = "Vegetation moisture content \n(NDMI HLS)")

## b. Plot fixed effects ----
(p <- model %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   mutate(`.variable` = gsub("b_", "", `.variable`) ) %>%
   # mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>%
   # build plot 
   ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.95),size = 4,fill = "#E8CEB6") +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   # scale_y_discrete(labels = mod_labs) +
   labs(x = "Coefficient estimate \n(scaled)",
        y = "") +
   theme_cowplot(18)
)


# Get random effects
ranef(model)$fireid %>% 
  as_tibble(rownames = "fireid")

## c. Plot random effects ----
(p_re <- model %>%
   gather_draws(`r_fireid\\[.*`, regex = TRUE) %>%          # filter estimated effects
   mutate(`.variable` = factor(`.variable`)) %>% 
   ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.05,0.95),size = 0.5,fill = '#b5ccb9') +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   theme_minimal_hgrid() + 
   theme(panel.grid = element_blank(),
         legend.position = "none") +
   labs(x = "Posterior estimates (standardized)", y = ""))


## d. Plot all mode effects over time ----

# Function to extract posteriors andfor each day and plot per predictor
plot_predictor_effects <- function(models, days, predictors, palette = NULL) {
  
  predictor_regex <- paste0("b_", predictors, collapse = "|")
  
  # Extract draws for all predictors
  all_effects <- map2_dfr(models, days, ~{
    gather_draws(.x, !!sym(predictor_regex), regex = TRUE) %>%
      mutate(day_before_fire = .y)
  })
  
  # Clean predictor names
  all_effects <- all_effects %>%
    mutate(predictor = gsub("^b_", "", .variable))
  
  # Add significance flag (does 95% CI exclude 0?)
  significance_flags <- all_effects %>%
    group_by(day_before_fire, predictor) %>%
    summarise(
      lower = quantile(.value, 0.025),
      upper = quantile(.value, 0.975),
      sig = !(lower < 0 & upper > 0),
      .groups = "drop"
    )
  
  # Join with full draws
  all_effects <- left_join(all_effects, significance_flags, by = c("day_before_fire", "predictor"))
  
  # Optional color palette
  if (is.null(palette)) {
    palette <- scales::hue_pal()(length(predictors))
    names(palette) <- predictors
  }
  
  # Plot with colored predictors
  ggplot(all_effects, aes(x = day_before_fire, y = .value,
                          fill = predictor, alpha = sig)) +
    stat_halfeye(.width = 0.95, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3)) +
    scale_fill_manual(values = palette) +
    labs(
      x = "Days before fire",
      y = "Posterior estimate",
      fill = "Predictor",
      alpha = "Significant"
    ) +
    scale_x_reverse() +
    facet_wrap(~ predictor, scales = "free_y") +
    theme_cowplot(14)
}

plot_predictor_effects(models = all_models, days = 1:30, predictors = predictors)

ggsave("figures/dnbr_effect_sizes_brms_benchmark_10pct.png",bg = "white",
       height = 10, width = 14)
