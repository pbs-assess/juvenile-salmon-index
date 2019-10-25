### Practice fitting survey catch data with sdmTMB
## Oct. 23 2019

library(tidyverse)
library(sdmTMB)
library(ggplot2)

browseVignettes("sdmTMB")

jchin <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  #remove stations that aren't present in both dataset
  filter(stableStation == "Y",
         #focus on summer for now
         # month %in% c(6, 7),
         #remove largest values so it converges
         !ck_juv > 100
         ) %>%
  #scale UTM coords
  mutate(xUTM_start = xUTM_start / 10000,
         yUTM_start = yUTM_start / 10000,
         jdayZ = scale(jday),
         bottomZ = scale(avg_bottom_depth)) %>% 
  #remove extra vars and stock ppn data
  select(-c(date, stableStation, samp_catch:SEAK))

jchin_spde <- make_spde(jchin$xUTM_start, jchin$yUTM_start, n_knots = 150)
plot_spde(jchin_spde)


qcs_grid %>% 
  nrow()

jchin %>% 
  filter(is.na(jdayZ))

ggplot(jchin, aes(x = jdayZ, y = ck_juv)) +
  geom_point()

## Develop index 
# Fit GLMM without covariates
m1 <- sdmTMB(ck_juv ~ 0 + as.factor(year), 
  data = jchin,
  time = "year", 
  spde = jchin_spde, 
  silent = FALSE,
  anisotropy = TRUE, 
  include_spatial = TRUE,
  ar1_fields = FALSE,
  family = nbinom2(link = "log"))
m2 <- sdmTMB(ck_juv ~ 0 + as.factor(year) + jdayZ, 
             data = jchin,
             time = "year", 
             spde = jchin_spde, 
             silent = FALSE,
             anisotropy = TRUE, 
             include_spatial = TRUE,
             ar1_fields = FALSE,
             family = nbinom2(link = "log"))

dum <- jchin %>% 
  mutate(resid = residuals(m1))

# Generate prediction grid 
surv_grid <- expand.grid(
  year = unique(jchin$year),
  jdayZ = seq(-2, 1.75, length.out = 100),
  X = seq(from = min(jchin$xUTM_start),
          to = max(jchin$xUTM_start),
          by = 2),
  Y = seq(from = min(jchin$yUTM_start),
          to = max(jchin$yUTM_start),
          by = 2)
)

predictions <- predict(m2, newdata = surv_grid, return_tmb_object = TRUE)

# Stolen from vignette...
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

# Predictions incorporating all fixed and random effects
plot_map(predictions$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

# Just fixed effects (not very meaningful here...)
plot_map(predictions$data, "exp(est_non_rf)") +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

# Spatial random effects (temporally stable factors driving changes in abundance)
plot_map(predictions$data, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

# Spatiotemporal random effects (dynamic drivers)
plot_map(predictions$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

ind <- get_index(predictions, bias_correct = FALSE)

scale <- 2 * 2  # 2 x 2 km grid 
ggplot(ind, aes(year, est*scale)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate')
