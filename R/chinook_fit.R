### Juvenile Chinook model fit
## Fit sdmTMB to 
## 1) Develop standardized index of abundance
## 2) Estimate vessel effects (requires developing interannual smooth or 
## AR-term, otherwise confounded with year effects)
## May 26, 2022


library(tidyverse)
library(sdmTMB)
library(ggplot2)


dat_trim <- readRDS(here::here("data", "chin_catch_sbc.rds"))
