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

range(jchin$xUTM_start)
range(jchin$yUTM_start)
range(jchin$ck_juv)
hist(jchin$ck_juv)
nrow(jchin)

jchin$xUTM_start <- jchin$xUTM_start/10000
jchin$yUTM_start <- jchin$yUTM_start/10000
jchin_spde <- make_spde(jchin$xUTM_start, jchin$yUTM_start, n_knots = 150)
plot_spde(jchin_spde)

## Develop index - FAILS TO CONVERGE 
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
m1


# plot(ck_juv ~ jday, jchin)
