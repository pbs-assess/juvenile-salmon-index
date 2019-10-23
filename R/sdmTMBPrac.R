### Practice fitting survey catch data with sdmTMB
## Oct. 23 2019

library(tidyverse)
library(sdmTMB)
library(ggplot2)

jchin <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  #remove stations that aren't present in both dataset
  filter(stableStation == "Y",
         #focus on summer for now
         # month %in% c(6, 7),
         #remove largest values
         !ck_juv > 100) %>%
  #remove extra vars and stock ppn data
  select(-c(date, stableStation, samp_catch:SEAK))

jchin_spde <- make_spde(jchin$xUTM_start, jchin$yUTM_start, n_knots = 100)
plot_spde(jchin_spde)


## Develop index - FAILS TO CONVERGE 
# Fit GLMM without spatiotemporal random fields or covariates
m1 <- sdmTMB(ck_juv ~ 0 + as.factor(year), data = jchin,
             time = "year", spde = jchin_spde,
             family = tweedie(link = "log"))

plot(ck_juv ~ jday, jchin)


## Fit spatial trend models
m1 <- sdmTMB(ck_juv ~ 1, data = jchin,
             spde = jchin_spde, family = tweedie(link = "log"),
             spatial_trend = TRUE, time = "year",
             spatial_only = FALSE)

hist(pcod$density)
hist(jchin$ck_juv)
