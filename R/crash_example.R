library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)


# make subdirectories for storage
dir.create("data/fits", recursive = TRUE, showWarnings = FALSE)

# downscale data and predictive grid
dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds")) 

# downscale data and predictive grid
dat <- dat_in %>% 
  mutate(
    year_f = as.factor(year),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3),
    scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    scale_yday = scale(as.numeric(yday))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night),
    bias = ifelse(year > 2015, "yes", "no")) %>% 
  filter(species == "coho") %>%
  droplevels()

spde <- make_mesh(
  dat,
  c("utm_x_1000", "utm_y_1000"),
  cutoff = 50
) 

fit <- sdmTMB(
  n_juv ~ 0 + season_f,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  # spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = FALSE,
  groups = "season_f",
  silent = FALSE
)

pp <- predict(fit, newdata = dat)






# Build a mesh to implement the SPDE approach:
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)

rr <- sample.int(3, size = nrow(pcod_2011), replace = T)
pcod_2011$bins <- as.factor(rr)

grid2 <- purrr::map(
  rr,
  ~ qcs_grid %>% 
    mutate(
      bins = .x
    )) %>% 
  bind_rows() %>% 
  mutate(bins = as.factor(bins))

# Fit a Tweedie spatial random field GLMM with a smoother for depth:
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod_2011, mesh = mesh,
  family = tweedie(link = "log"),
  groups = "bins",
  spatial = "on"
)
# Predict on new data:
p <- predict(fit, newdata = grid2)
