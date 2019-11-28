## Fit GSI data to multinomial model
# Nov. 28, 2019
# Note: Currently assumes certain assignment

library(tidyverse)

juv <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  filter(stableStation == "Y") %>%
  #scale UTM coords
  mutate(#xUTM = xUTM_start / 10000,
         #yUTM = yUTM_start / 10000,
         jdayZ = as.vector(scale(jday)[,1]),
         season = as.factor(
           case_when(
             month %in% c("2", "3") ~ "winter",
             month %in% c("5", "6", "7", "8") ~ "summer",
             month %in% c("9", "10", "11" , "12") ~ "fall")
         )
  ) %>% 
  #remove extra vars and stock ppn data
  select(station_id, dataset, jday:year, jdayZ, season, start_lat, start_long, 
         ck_juv:SEAK) %>% 
  arrange(jdayZ)

compMat <- 
