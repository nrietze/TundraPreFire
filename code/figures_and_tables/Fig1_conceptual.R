library(ggplot2)
library(cowplot)
library(viridis)

# 1. Simulate NDMI time series ----
# Define parameters
doy <- 165:200  # Day of Year range
burn_day <- 180 # Fire occurs at DOY 180

# Function to simulate NDMI curves
simulate_ndmi <- function(doy, burn_day, type) {
  pre_fire <- doy < burn_day
  post_fire <- doy >= burn_day
  
  ndmi <- numeric(length(doy))
  
  if (type == "a") {
    # Increase before fire, drop and level out
    ndmi[pre_fire] <- seq(0.1, 0.3, length.out = sum(pre_fire)) + rnorm(sum(pre_fire), 0, 0.05)
    ndmi[post_fire] <- seq(-0.1, -0.05, length.out = sum(post_fire)) + rnorm(sum(post_fire), 0, 0.05)
  }
  
  if (type == "b") {
    # Stable before fire, slight drop after fire
    ndmi[pre_fire] <- rep(0.1, sum(pre_fire)) + rnorm(sum(pre_fire), 0, 0.03)
    ndmi[post_fire] <- seq(0.0, 0.05, length.out = sum(post_fire)) + rnorm(sum(post_fire), 0, 0.03)
  }
  
  if (type == "c") {
    # Decrease before fire, slight drop after fire, then level out
    ndmi[pre_fire] <- seq(0.2, 0.0, length.out = sum(pre_fire)) + rnorm(sum(pre_fire), 0, 0.05)
    ndmi[post_fire] <- seq(-0.1, 0, length.out = sum(post_fire)) + rnorm(sum(post_fire), 0, 0.05)
  }
  
  # Constrain NDMI values between -0.2 and 0.4
  ndmi <- pmax(-0.2, pmin(ndmi, 0.4))
  
  return(data.frame(DaysSinceBurn = doy - burn_day, NDMI = ndmi, Curve = type))
}

# Simulate NDMI for each curve
ndmi_data <- rbind(
  simulate_ndmi(doy, burn_day, "a"),
  simulate_ndmi(doy, burn_day, "b"),
  simulate_ndmi(doy, burn_day, "c")
)

# Create bars at days -4 to -1
bars <- expand.grid(DaysSinceBurn = -4:-1, Curve = c("a", "b", "c")) %>%
  mutate(ymin = -0.2, ymax = 0.4)  # Bar height

# Figure 1 a)
ggplot() +
  geom_smooth(data = ndmi_data, aes(x = DaysSinceBurn, y = NDMI, color = Curve),
              method = "loess", span = 0.3, se = FALSE, size = 1) + 
  geom_rect(data = bars, aes(xmin = DaysSinceBurn - 0.5, xmax = DaysSinceBurn + 0.5, 
                             ymin = ymin, ymax = ymax, fill = Curve), 
            color = NA, alpha = 0.4) +  # Bars for -4 to -1 days
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  labs(title = "Conceptual Methodology",
       x = "Days since burn", y = "NDMI",
       color = "Curve Type") +
  scale_color_viridis_d(end = 0.8, option = "inferno") +
  lims(x = c(-55,30),y = c(-0.2,0.4)) +
  facet_wrap(~Curve) +
  theme_cowplot()
