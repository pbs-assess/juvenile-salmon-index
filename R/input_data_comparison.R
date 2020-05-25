### Fit subsets of survey catch data to determine how changes in sampling
# effort will impact annual abundance estimates
# Model fits and some figures from day/night and spatial strat sections are 
# incorporated into subset_comparison.Rmd
## April 29, 2020

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
         season = as.factor(
           case_when(
            month %in% c("2", "3") ~ "winter",
            month %in% c("5", "6", "7", "8") ~ "summer",
            month %in% c("9", "10", "11" , "12") ~ "fall")
           ),
         yday_z = as.vector(scale(yday)[,1]),
         yday_z2 = yday_z^2,
         time_f = as.factor(time_f),
         juv_cpue = ck_juv / dur,
         offset = dur
         ) %>% 
  #remove extra vars and stock ppn data
  dplyr::select(station_id, stable_station:date, time_f, yday, yday_z, yday_z2,
                offset, 
                month:dur, head_depth,
                ck_juv, juv_cpue)


#-------------------------------------------------------------------------------
## Test for day/night impacts
ipes_only <- jchin %>% 
  filter(synoptic == "1",
         head_depth < 21,
         month %in% c(6, 7),
         !year == "2014",
         dataset == "IPES")

library(glmmTMB)

#fit various models with distributions appropriate for zero-inflated integer 
#data
fit_zibin <- glmmTMB(ck_juv ~ time_f + (1|year), data = ipes_only, 
                     ziformula = ~1 + (1|year), family = nbinom2)
fit_pois <- glmmTMB(ck_juv ~ time_f + (1|year), data = ipes_only, 
                     ziformula = ~1, family = poisson)
summary(fit_zibin)
summary(fit_pois)
AIC(fit_zibin, fit_pois)

# make predictions
preds <- predict(fit_zibin,
        newdata = data.frame(time_f = unique(ipes_only$time_f),
                             year = NA),
        se.fit = TRUE,
        type = "response") 
data.frame(time_f = unique(ipes_only$time_f),
                       fit = preds$fit,
                       se = preds$se.fit) %>% 
  mutate(low = fit + qnorm(0.025) * se,
         up = fit + qnorm(0.975) * se) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = time_f, y = fit, ymin = low, ymax = up)) +
  labs(x = "Sampling Time", y = "Predicted Abundance Per Tow") +
  ggsidekick::theme_sleek()


## alternative format using brms
library(brms)

fit_zibinB <- brm(ck_juv ~ time_f + (1|year), data = ipes_only, 
                  family = zero_inflated_negbinomial(),
                  chains = 4, cores = 4, control = list(adapt_delta = 0.97))
plot(marginal_effects(fit_zibinB))
summary(fit_zibinB)

predsB <- predict(fit_zibinB,
                 newdata = data.frame(time_f = unique(ipes_only$time_f),
                                      year = NA),
                 se.fit = TRUE,
                 type = "response")

data.frame(time_f = unique(ipes_only$time_f),
           fit = predsB[ , "Estimate"],
           se = predsB[ , "Est.Error"]) %>% 
  mutate(low = fit + qnorm(0.025) * se,
         up = fit + qnorm(0.975) * se) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = time_f, y = fit, ymin = low, ymax = up)) +
  labs(x = "Sampling Time", y = "Predicted Abundance Per Tow") +
  ggsidekick::theme_sleek()



# ------------------------------------------------------------------------------
## Test for spatial stratification impacts

#generate two datasets:
# 1) limited to spatially overlapping tows and <20 m headrope depth, but includes
# inlets/nocturnal tows and all seasons
# 2) limited to above, plus daytime tows only in June/July, in overlapping strata
jchin1 <- jchin %>% 
  filter(stable_station == "1",
         month %in% c(6, 7),
         head_depth < 21,
         !juv_cpue > 500 #remove large values that may be skewing results
         ) %>% 
  mutate(trim_size = "large") 
jchin2 <- jchin %>% 
  filter(synoptic == "1",
         head_depth < 21,
         month %in% c(6, 7),
         !year == "2014"
         ) %>% 
  mutate(trim_size = "small") 


jchin %>% 
  filter(stable_station == "1",
         month %in% c(6, 7),
         head_depth < 21,
         juv_cpue > 500 #remove large values that may be skewing results
  ) 
# compare spatial distributions
library(ggmap)
nAm <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))
rbind(jchin1, jchin2) %>% 
  ggplot(.) +
  geom_point(aes(x = start_long, y = start_lat, color = dataset)) +
  geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region),
           color = "black", fill = "gray80") +
  coord_fixed(xlim = c(-129.5, -123), ylim = c(48, 52), ratio = 1.3) +
  facet_wrap(~trim_size) +
  ggsidekick::theme_sleek()
  

jchin1_spde <- make_spde(jchin1$xUTM_start, jchin1$yUTM_start, n_knots = 150)
jchin2_spde <- make_spde(jchin2$xUTM_start, jchin2$yUTM_start, n_knots = 150)
plot_spde(jchin1_spde)
plot_spde(jchin2_spde)

ggplot(jchin1, aes(x = yday, y = juv_cpue, colour = as.factor(dataset))) +
  geom_point()
ggplot(jchin2, aes(x = yday, y = juv_cpue, colour = as.factor(dataset))) +
  geom_point()

jchin2 %>% 
  group_by(year) %>% 
  summarize(n = length(unique(station_id)), 
            mean_catch = mean(ck_juv)) %>% 
  print(n = Inf)

## Develop index 
# Yearly model with tweedie distribution 
dir.create(file.path(here::here("data", "modelFits")))
# m1_tw <- sdmTMB(juv_cpue ~ 0 + as.factor(year) + yday_z + yday_z2,
#                 data = jchin1,
#                 time = "year",
#                 spde = jchin1_spde,
#                 silent = FALSE,
#                 anisotropy = TRUE,
#                 include_spatial = TRUE,
#                 ar1_fields = FALSE,
#                 family = tweedie(link = "log"))
# saveRDS(m1_tw, here::here("data", "modelFits", "ck1_yrInt_tw.rds"))
# negative binomial model with effort offset favored by SA
m1_nb <- sdmTMB(ck_juv ~ 0 + as.factor(year) + yday_z + yday_z2 + offset,
             data = jchin1,
             time = "year",
             spde = jchin1_spde,
             silent = FALSE,
             anisotropy = TRUE,
             include_spatial = TRUE,
             ar1_fields = FALSE,
             family = nbinom2(link = "log"))
saveRDS(m1_nb, here::here("data", "modelFits", "ck1_yrInt_nb.rds"))

# m2_tw <- sdmTMB(juv_cpue ~ 0 + as.factor(year) + yday_z + yday_z2,
#                 data = jchin2,
#                 time = "year",
#                 spde = jchin2_spde,
#                 silent = FALSE,
#                 anisotropy = TRUE,
#                 include_spatial = TRUE,
#                 ar1_fields = FALSE,
#                 family = tweedie(link = "log"))
# saveRDS(m2_tw, here::here("data", "modelFits", "ck2_yrInt_tw.rds"))
m2_nb <- sdmTMB(ck_juv ~ 0 + as.factor(year) + yday_z + yday_z2 + offset,
                data = jchin2,
                time = "year",
                spde = jchin2_spde,
                silent = FALSE,
                anisotropy = TRUE,
                include_spatial = TRUE,
                ar1_fields = FALSE,
                family = nbinom2(link = "log"))
saveRDS(m2_nb, here::here("data", "modelFits", "ck2_yrInt_nb.rds"))


# examine residuals
check_res <- function(mod, dat) {
  dat$resid <- residuals(mod)
  #histogram
  hist(dat$resid)
  # qqplot
  qqnorm(dat$resid)
  abline(a = 0, b = 1, col = "red")
  # map of resids
  ggplot(dat, aes(xUTM_start, yUTM_start, col = resid)) + 
    scale_colour_gradient2() +
    geom_point() + facet_wrap(~year) + coord_fixed()
}

# check_res(m1_tw, jchin1)
# check_res(m2_tw, jchin2)
check_res(m1_nb, jchin1)
check_res(m2_nb, jchin2)


# Make predictions 
# generate a spatial grid for each dataset (necessary since they have different
# extents)
source(here::here("R", "prepPredictionGrid.R"))
jchin_list <- list(jchin1, jchin2) 
pred_grid_list <- map(jchin_list, function(x) {
  x %>% 
    mutate(xUTM_start = xUTM_start * 10000, #scale up to match reality
           yUTM_start = yUTM_start * 10000) %>% 
    make_pred_grid() %>% #make the grid
    mutate(X = X / 10000, #scale back up
           Y = Y / 10000) %>% 
    expand(nesting(X, Y), #make a predictive dataframe
           year = unique(x$year)) %>%
    mutate(yday_z = 0,
           yday_z2 = 0,
           offset = mean(x$dur))
})
saveRDS(pred_grid_list[[1]], here::here("data", "pred_grids",
                                        "large_ck_grid.rds"))
saveRDS(pred_grid_list[[2]], here::here("data", "pred_grids",
                                        "small_ck_grid.rds"))

pred_m1 <- predict(m1_nb, newdata = pred_grid_list[[1]], return_tmb_object = TRUE)
pred_m2 <- predict(m2_nb, newdata = pred_grid_list[[2]], return_tmb_object = TRUE)
glimpse(pred_m1$data)

# Compare predictions from the two datasets
ind1 <- get_index(pred_m1, bias_correct = FALSE) %>% 
  mutate(dataset = "large")
ind2 <- get_index(pred_m2, bias_correct = FALSE) %>% 
  #add gap for 2014
  rbind(., data.frame(year = 2014, est = NA, lwr = NA, upr = NA, log_est = NA,
                      se = NA, max_gradient = NA, bad_eig = FALSE)) %>% 
  mutate(dataset = "small")  
inds <- rbind(ind1, ind2)

scale <- 2 * 2  # 2 x 2 km grid 
ggplot(inds, aes(year, est*scale)) + 
  geom_line(aes(group = 1)) +
  geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4) +
  xlab('Year') + ylab('Abundance Estimate') +
  facet_wrap(~dataset)

## Ignore following for now

# Stolen from vignette...
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

# Predictions incorporating all fixed and random effects
plot_map(pred_m2$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

# Just fixed effects (not very meaningful here...)
plot_map(pred_m2$data, "exp(est_non_rf)") +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

# Spatial random effects (temporally stable factors driving changes in abundance)
plot_map(pred_m2$data, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

# Spatiotemporal random effects (dynamic drivers)
plot_map(pred_m2$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

