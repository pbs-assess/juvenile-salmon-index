## Fit GSI data to multinomial model
# Nov. 28, 2019
# Note: Currently assumes certain assignment

library(tidyverse)

# identify focal stations from different dataset (FIX EVENTUALLY)
juv <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  filter(stableStation == "Y")
gsiLongAgg <- readRDS(here::here("data", "longGSI_reg4.rds"))
gsiLongAgg %>% 
  filter(station_id %in% unique(juv$station_id)) %>% 
  group_by(Region4Name) %>% 
  tally()

## FILTER AND ADJUST 
gsiWide <- gsiLongAgg %>% 
  select(-ship_fl, -c(aggProb:maxProb)) %>%
  mutate(present = 1) %>% 
  pivot_wider(., names_from = Region4Name, values_from = present) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  glimpse()

gsi <- gsiWide %>% 
  filter(station_id %in% unique(juv$station_id),
         age == "J") %>%
  #scale UTM coords
  mutate(
    season = as.factor(
      case_when(
        month %in% c("2", "3") ~ "winter",
        month %in% c("5", "6", "7", "8") ~ "summer",
        month %in% c("9", "10", "11" , "12") ~ "fall")
    )
  ) %>% 
  #remove extra vars and stock ppn data
  select(fish_number, jday, season, year, SalSea:SEAK) %>% 
  arrange(year)

glimpse(gsi)

## Obs
y_obs <- gsi %>% 
  filter(year %in% ("2018")) %>% 
  select(SalSea:SEAK) %>% 
  as.matrix()
  
library(TMB)
compile("R/multinomialPractice/multinomial_generic.cpp")
dyn.load(dynlib("R/multinomialPractice/multinomial_generic"))

## Data and parameters
# .X <- cbind(1, X) #predictor with intercept
.X <- cbind(1, X) #predictor with intercept
data <- list(cov=.X, y_obs = y_obs)
parameters <- list(betas = matrix(data = 0, nrow = ncol(.X), 
                                  ncol = ncol(y_obs) - 1))