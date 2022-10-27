## Spline Knots
# Quick evaluation of effects of spline type and knot specification on marginal
# effects predictions
# Oct. 27, 2022

library(tidyverse)

# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  filter(species == "CHUM")

week_dat <- data.frame(
  week = seq(0, 52, length.out = 100),#seq(5, 48, length.out = 100),
  dist_to_coast_km = median(dat$dist_to_coast_km),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year = 2011L
)

fit1 <- mgcv::gam(
  n_juv ~ 1 +
    s(week, bs = "tp", k = 5) +
    survey_f,
  data = dat,
  family = mgcv::nb())
fit2 <- mgcv::gam(
  n_juv ~ 1 +
    s(week, bs = "cc", k = 5) +
    survey_f,
  data = dat,
  family = mgcv::nb())
fit3 <- mgcv::gam(
  n_juv ~ 1 +
    s(week, bs = "cc", k = 5) +
    survey_f,
  data = dat,
  family = mgcv::nb(),
  knots = list(week = c(0, 52)))

p1 <- predict(fit1, newdata = week_dat, se_fit = T, re_form = NA)
p2 <- predict(fit2, newdata = week_dat, se_fit = T, re_form = NA)
p3 <- predict(fit3, newdata = week_dat, se_fit = T, re_form = NA)

week_pred <- week_dat %>% 
  mutate(
    p1 = as.numeric(p1), p2 = as.numeric(p2), p3 = as.numeric(p3)
  ) %>% 
  pivot_longer(
    cols = c(p1:p3), names_to = "model", values_to = "est"
  )

ggplot(week_pred, 
       aes(week, exp(est), colour = model)) +
  geom_line() 
