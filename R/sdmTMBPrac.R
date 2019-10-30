### Practice fitting survey catch data with sdmTMB
## Oct. 23 2019

library(tidyverse)
library(sdmTMB)
library(ggplot2)

# browseVignettes("sdmTMB")

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
         jdayZ = as.vector(scale(jday)[,1]),
         jdayZ2 = jdayZ^2,
         bottomZ = as.vector(scale(avg_bottom_depth)[,1])) %>% 
  #remove extra vars and stock ppn data
  select(-c(date, stableStation, samp_catch:SEAK))

jchin_spde <- make_spde(jchin$xUTM_start, jchin$yUTM_start, n_knots = 150)
plot_spde(jchin_spde)

ggplot(jchin, aes(x = jdayZ2, y = ck_juv)) +
  geom_point()

## Develop index 
# Fit GLMM without covariates
# m1 <- sdmTMB(ck_juv ~ 0 + as.factor(year), 
#   data = jchin,
#   time = "year", 
#   spde = jchin_spde, 
#   silent = FALSE,
#   anisotropy = TRUE, 
#   include_spatial = TRUE,
#   ar1_fields = FALSE,
#   family = nbinom2(link = "log"))

## Daily model
mDay <- readRDS(here::here("data", "modelFits", "dayModel.rds"))
mDay <- sdmTMB(ck_juv ~ 0 + as.factor(year) + jdayZ + jdayZ2,
                 data = jchin,
                 time = "year",
                 spde = jchin_spde,
                 silent = FALSE,
                 anisotropy = TRUE,
                 include_spatial = TRUE,
                 ar1_fields = FALSE,
                 family = nbinom2(link = "log"))
saveRDS(mDay, here::here("data", "modelFits", "dayModel.rds"))

# Prediction grid (removing subannual daily effect)
# Inappropriate because generating predictions for land masses but fine for now
# surv_grid_days <- expand.grid(
#   X = seq(from = round(min(jchin$xUTM_start)),
#           to = round(max(jchin$xUTM_start)),
#           by = 1),
#   Y = seq(from = round(min(jchin$yUTM_start)),
#           to = round(max(jchin$yUTM_start)),
#           by = 1),
#   year = unique(jchin$year)
# ) %>%
#   mutate(jdayZ = 0,
#          jdayZ2 = 0)

surv_grid_days <- readRDS(here::here("data", "spatialData", 
                                     "trimmedSurveyGrid.rds")) %>% 
  expand(nesting(X, Y), year = unique(jchin$year)) %>%
  mutate(jdayZ = 0,
         jdayZ2 = 0)

# Too many gaps in space
# surv_grid_days <- jchin %>%
#   expand(nesting(xUTM_start, yUTM_start), year) %>%
#   mutate(jdayZ = 0,
#          jdayZ2 = 0) %>%
#   rename(X = xUTM_start, Y = yUTM_start)

pred_m <- predict(mDay, newdata = surv_grid_days, return_tmb_object = TRUE)
glimpse(pred_m$data)

## Ignore depth for now
# mDepth <- sdmTMB(ck_juv ~ 0 + as.factor(year) + bottomZ, 
#              data = jchin,
#              time = "year", 
#              spde = jchin_spde, 
#              silent = FALSE,
#              anisotropy = TRUE, 
#              include_spatial = TRUE,
#              ar1_fields = FALSE,
#              family = nbinom2(link = "log"))
# saveRDS(mDepth, here::here("data", "modelFits", "depthModel.rds"))
# mDepth <- readRDS(here::here("data", "modelFits", "depthModel.rds"))
# surv_grid <- jchin %>%
#   expand(nesting(xUTM_start, yUTM_start, bottomZ), year) %>% 
#   rename(X = xUTM_start, Y = yUTM_start)

# waiting on patch for neg binomial residuals
dum <- jchin %>% 
  mutate(resid = residuals(m1))

# Stolen from vignette...
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile(width = 1, height = 1) +
    facet_wrap(~year) +
    coord_fixed()
}

# Predictions incorporating all fixed and random effects
plot_map(pred_m$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

# Just fixed effects (not very meaningful here...)
plot_map(pred_m$data, "exp(est_non_rf)") +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

# Spatial random effects (temporally stable factors driving changes in abundance)
plot_map(pred_m$data, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

# Spatiotemporal random effects (dynamic drivers)
plot_map(pred_m$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

ind <- get_index(pred_m, bias_correct = FALSE)

scale <- 2 * 2  # 2 x 2 km grid 
ggplot(ind, aes(year, est*scale)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate')
 
