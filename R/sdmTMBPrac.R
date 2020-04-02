### Practice fitting survey catch data with sdmTMB
## Oct. 23 2019

library(tidyverse)
library(sdmTMB)
library(ggplot2)

# browseVignettes("sdmTMB")

bridge <- readRDS(here::here("data", "ipes_hs_merged_bridge.rds")) 
# jchin <- readRDS(here::here("data", "juvCatchGSI_reg4.rds"))

# use bridge data alone for now since ignoring stock composition
jchin <- bridge %>%
  #scale UTM coords
  mutate(xUTM_start = xUTM_start / 10000,
         yUTM_start = yUTM_start / 10000,
         yday_z = as.vector(scale(yday)[,1]),
         yday_z2 = yday_z^2,
         season = as.factor(
           case_when(
            month %in% c("2", "3") ~ "winter",
            month %in% c("5", "6", "7", "8") ~ "summer",
            month %in% c("9", "10", "11" , "12") ~ "fall")
           ),
         time_f = as.factor(time_f)
         ) %>% 
  #remove extra vars and stock ppn data
  dplyr::select(station_id, stable_station:date, time_f:dur, head_depth, ck_juv)

#generate two datasets:
# 1) limited to spatially overlapping tows and <20 m headrope depth, but includes
# inlets/nocturnal tows and all seasons
# 2) limited to above, plus daytime tows only in June/July, in overlapping strata
jchin1 <- jchin %>% 
  filter(stable_station == "1",
         head_depth < 21)
jchin2 <- jchin %>% 
  filter(synoptic == "1",
         head_depth < 21,
         month %in% c(6, 7),
         time_f == "day")

jchin_spde1 <- make_spde(jchin1$xUTM_start, jchin1$yUTM_start, n_knots = 150)
jchin_spde2 <- make_spde(jchin2$xUTM_start, jchin2$yUTM_start, n_knots = 150)
plot_spde(jchin_spde1)
plot_spde(jchin_spde2)

ggplot(jchin1, aes(x = yday, y = ck_juv, colour = as.factor(month))) +
  geom_point()
ggplot(jchin2, aes(x = yday, y = ck_juv, colour = as.factor(month))) +
  geom_point()


## Develop index 
## Daily model
dir.create(file.path(here::here("data", "modelFits")))
m1_nb <- sdmTMB(ck_juv ~ 0 + as.factor(year) + jdayZ + jdayZ2,
                 data = jchin,
                 time = "year",
                 spde = jchin_spde,
                 silent = FALSE,
                 anisotropy = TRUE,
                 include_spatial = TRUE,
                 ar1_fields = FALSE,
                 family = nbinom2(link = "log"))
saveRDS(m1_nb, here::here("data", "modelFits", "day_nb.rds"))

## Daily model incorporating seasonal effects and quadratics
m2_nb <- sdmTMB(ck_juv ~ 0 + as.factor(year) + season:jdayZ + season:jdayZ2,
               data = jchin,
               time = "year",
               spde = jchin_spde,
               silent = FALSE,
               anisotropy = TRUE,
               include_spatial = TRUE,
               ar1_fields = FALSE,
               family = nbinom2(link = "log"))
saveRDS(m2_nb, here::here("data", "modelFits", "day_season_nb.rds"))

## as above but with tweedie
m2_tw <- sdmTMB(ck_juv ~ 0 + as.factor(year) + season:jdayZ + season:jdayZ2,
             data = jchin,
             time = "year",
             spde = jchin_spde,
             silent = FALSE,
             anisotropy = TRUE,
             include_spatial = TRUE,
             ar1_fields = FALSE,
             family = tweedie(link = "log"))
saveRDS(m2_tw, here::here("data", "modelFits", "day_season_nb.rds"))

# examine residuals
jchin$resid_tw <- residuals(m2_tw)
hist(jchin$resid_tw)
ggplot(jchin, aes(xUTM_start, yUTM_start, col = resid_tw)) + 
  scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()

qqnorm(jchin$resid_tw)
abline(a = 0, b = 1, col = "red")


# Prediction grid (removing subannual daily effect)
surv_grid <- readRDS(here::here("data", "trimmedSurveyGrid.rds")) %>% 
  #scale down
  mutate(X = X / 10000,
         Y = Y / 10000) %>% 
  expand(nesting(X, Y), year = unique(jchin$year), 
         season = unique(jchin$season)) %>%
  mutate(jdayZ = 0,
         jdayZ2 = 0)

pred_m <- predict(m2_tw, newdata = surv_grid, return_tmb_object = TRUE)
glimpse(pred_m$data)


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

ind <- get_index(pred_m, bias_correct = FALSE)

scale <- 2 * 2  # 2 x 2 km grid 
ggplot(ind, aes(year, est*scale)) + 
  geom_line() +
  geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate')
 

## Map of set locations as reference to prediction grid
ggplot(jchin) +
  geom_point(aes(x = xUTM_start, y = yUTM_start)) +
  coord_fixed(xlim = c(37, 82), ylim = c(530, 588), ratio = 1.3)
