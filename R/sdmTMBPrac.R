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
         month %in% c(6, 7),
         #remove largest values so it converges
         !ck_juv > 100,
         !is.na(avg_bottom_depth)
         ) %>%
  #scale UTM coords
  mutate(xUTM_start = xUTM_start / 10000,
         yUTM_start = yUTM_start / 10000,
         jdayZ = as.vector(scale(jday)[,1]),
         bottomZ = as.vector(scale(avg_bottom_depth)[,1])) %>% 
  #remove extra vars and stock ppn data
  select(-c(date, stableStation, samp_catch:SEAK))

jchin_spde <- make_spde(jchin$xUTM_start, jchin$yUTM_start, n_knots = 150)
plot_spde(jchin_spde)

ggplot(jchin, aes(x = bottomZ, y = ck_juv)) +
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


mDepth <- sdmTMB(ck_juv ~ 0 + as.factor(year) + bottomZ, 
             data = jchin,
             time = "year", 
             spde = jchin_spde, 
             silent = FALSE,
             anisotropy = TRUE, 
             include_spatial = TRUE,
             ar1_fields = FALSE,
             family = nbinom2(link = "log"))
# saveRDS(mDepth, here::here("data", "modelFits", "depthModel.rds"))
mDepth <- readRDS(here::here("data", "modelFits", "depthModel.rds"))

dum <- jchin %>% 
  mutate(resid = residuals(m1))
# waiting on patch for neg binomial residuals

# Generate prediction grid 
surv_grid <- jchin %>%
  expand(nesting(xUTM_start, yUTM_start, bottomZ), year) %>% 
  rename(X = xUTM_start, Y = yUTM_start)

pred_m <- predict(mDepth, newdata = surv_grid, return_tmb_object = TRUE)

# Stolen from vignette...
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
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

ind <- get_index(pred_m2, bias_correct = FALSE)

scale <- 2 * 2  # 2 x 2 km grid 
ggplot(ind, aes(year, est*scale)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate')


## Structure of daily model incorrect - likely need to fit separate models
# by survey
# mDay <- sdmTMB(ck_juv ~ 0 + as.factor(year) + jdayZ, 
#                  data = jchin,
#                  time = "year", 
#                  spde = jchin_spde, 
#                  silent = FALSE,
#                  anisotropy = TRUE, 
#                  include_spatial = TRUE,
#                  ar1_fields = FALSE,
#                  family = nbinom2(link = "log"))
# saveRDS(mDay, here::here("data", "modelFits", "dayModel.rds"))

# daySeq <- c(-1.5, -0.2, 1.1)
# surv_grid_days <- jchin %>% 
#   expand(nesting(xUTM_start, yUTM_start), year, daySeq) %>% 
#   rename(X = xUTM_start, Y = yUTM_start, jdayZ = daySeq)
# pred_m <- predict(mDay, newdata = surv_grid_days, return_tmb_object = TRUE)

