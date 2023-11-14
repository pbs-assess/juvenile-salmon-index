## Survey Design Impacts
# Evaluate changes in IPES survey design (number of stations and day/night
# format)
# 1. Randomly sample 2000 "stations" within the predictive grid domain
# 2. Use species specific models to simulate catch estimates (day and night) at 
# each station for 3 years (n = 10-20 iterations initially)
# 3. Sample the simulated catch with different survey designs
#   a. 50 stations each sampled day/night
#   b. 35 stations each sampled day/night
#   c. 100 stations sampled randomly (2/3 day, 1/3 night)
#   d. 70 stations sampled randomly (2/3 day, 1/3 night)
# 4. Fit the model to each dataset and estimate index
# 5. Quantify changes in CV of index estimates across survey scenarios
# ASSUMPTIONS:
# - no seasonality
# - day/night samples separate
# - spatial clustering of stations not accounted for


library(sdmTMB)
library(sdmTMBextra)
library(tidyverse)

# ensure branch with future sims is installed
devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "future-sims")

dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds")) 

yr_key <- data.frame(
  year = c(unique(dat_in$year), 2023, 2024, 2025)
) %>% 
  arrange(year) %>% 
  mutate(
    sox_cycle = rep(1:4, length.out = 28L) %>% 
      as.factor(),
    pink_cycle = rep(1:2, length.out = 28L) %>% 
      as.factor()
  )

# downscale data and predictive grid
dat <- dat_in %>% 
  mutate(
    year_f = as.factor(year),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3 * 500),
    scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    day_night = as.factor(day_night)) %>% 
  filter(!species == "steelhead") %>% 
  droplevels() %>% 
  left_join(., yr_key, by = "year") %>% 
  filter(species == "sockeye", 
         season_f == "su")

# mesh
dat_coords <- dat %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(2, 10) * 500,
  cutoff = 30,
  offset = c(10, 50)
)  
spde <- make_mesh(
  dat,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 


fit <- sdmTMB(
  n_juv ~ day_night + survey_f + scale_dist + scale_depth + sox_cycle,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  # spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  extra_time = c(2016, 2020),
  # groups = "season_f",
  # control = sdmTMBcontrol(
  #   map = list(
  #     ln_tau_Z = factor(
  #       rep(1, times = length(unique(dat$season_f)))
  #     )
  #   )
  # ),
  silent = FALSE
)


## SIMULATE DATA ---------------------------------------------------------------

extra_yrs <- c(2023, 2024, 2025)

# append extra_yrs to epsilon
p <- get_pars(fit)
eps <- p$epsilon_st
eps2 <- array(dim = dim(eps) + c(0, length(extra_yrs), 0))
eps2[,,1] <- cbind(
  eps[,,1], 
  matrix(NA, nrow = dim(eps)[1], ncol = length(extra_yrs))
)

# predictive grid
grid_in <- readRDS(
  here::here("data", "spatial", "pred_ipes_grid.RDS")
  )$ipes_grid %>% 
  mutate(
    utm_x_1000 = X / 1000,
    utm_y_1000 = Y / 1000,
    dist_to_coast_km = shore_dist / 1000,
    cell_id = row_number(),
    target_depth = 0,
    survey_f = as.factor("ipes"),
    scale_depth = (target_depth - mean(dat$target_depth)) /
      sd(dat$target_depth),
    scale_dist = (dist_to_coast_km - mean(dat$dist_to_coast_km)) /
      sd(dat$dist_to_coast_km)
  ) 

full_grid <- replicate_df(
  grid_in, "day_night", time_values = unique(dat$day_night)
  ) %>% 
  replicate_df(
    ., "year", time_values = extra_yrs
  ) %>% 
  mutate(
    full_grid = TRUE,
    effort = mean(dat$effort)
  ) 

set.seed(123)
samp_foo <- function(draws, type = c("paired", "unpaired")) {
  out_dat <- vector(mode = "list", length = length(extra_yrs))
  if (type == "paired") {
    for (i in seq_along(extra_yrs)) {
      out_dat[[i]] <- data.frame(
        cell_id = sample.int(nrow(grid_in), draws),
        year = extra_yrs[i]
        )
    }
  }
  if (type == "unpaired") {
    for (i in seq_along(extra_yrs)) {
      out_dat[[i]] <- data.frame(
        cell_id = sample.int(nrow(grid_in), draws),
        day_night = sample(c(rep("DAY", .6 * draws),
                             rep("NIGHT", .4 * draws))),
        year = extra_yrs[i]
      )
    }
  }
  tt <- bind_rows(out_dat)
} 



future_obs <- full_grid %>% 
  filter(cell_id %in% scen_1) %>% 
  select(intersect(colnames(.), colnames(dat))) 
all_obs <- rbind(
  dat %>% 
    select(colnames(future_obs)),
  future_obs
) %>% 
  mutate(
    full_grid = FALSE
  )

# join observations to full grid when simulating
new_dat <- rbind(
  full_grid %>% 
    select(colnames(all_obs)),
  all_obs
) %>% 
  left_join(
    ., yr_key, by = "year"
  ) 
new_mesh <- make_mesh(new_dat, c("utm_x_1000", "utm_y_1000"),
                      mesh = spde$mesh)

b <- tidy(fit, "ran_pars")
s <- sdmTMB_simulate(
  ~ day_night + survey_f + scale_dist + scale_depth + sox_cycle,
  data = new_dat,
  mesh = new_mesh,
  range = b$estimate[b$term == "range"],
  sigma_E = b$estimate[b$term == "sigma_E"],
  sigma_O = b$estimate[b$term == "sigma_O"],
  phi = b$estimate[b$term == "phi"],
  B = unname(coef(fit)),
  offset = new_dat$effort,
  time = "year",
  spatiotemporal = "rw",
  family = sdmTMB::nbinom2(),
  seed = 2927819,
  fixed_re = list(omega_s = p$omega_s, epsilon_st = eps2) 
)
s2 <- sB %>% 
  mutate(
  scale_depth = new_dat$scale_depth,
  scale_dist = new_dat$scale_dist,
  survey_f = new_dat$survey_f,
  full_grid = new_dat$full_grid,
  day_night = new_dat$day_night,
  sox_cycle = new_dat$sox_cycle,
  effort = new_dat$effort
)

sim_grid <- s2 %>% 
  filter(full_grid == TRUE,
         day_night == "DAY")
sim <- s2 %>% 
  filter(full_grid == FALSE,
         year > max(dat$year))

# ggplot(sim, aes(utm_x_1000, utm_y_1000, colour = mu)) + geom_point() +
#   scale_color_viridis_c(trans = "log10") +
#   facet_wrap(~year)+
#   coord_fixed()
# ggplot(sim_grid %>% filter(day_night=="DAY"), 
#        aes(utm_x_1000, utm_y_1000, fill = mu)) + geom_raster() +
#   scale_fill_viridis_c(trans = "sqrt") +
#   facet_wrap(~year) +
#   coord_fixed()

sim_future <- sim  %>% 
  select(year, utm_x_1000, utm_y_1000, day_night, survey_f, n_juv = observed,
         effort,
         scale_dist, scale_depth, sox_cycle) %>%  
  mutate(simulated_data = TRUE)
combined_dat <- bind_rows(
  dat %>% 
    mutate(simulated_data = FALSE) %>% 
    select(colnames(sim_future)), 
  sim_future
)

new_mesh2 <- make_mesh(combined_dat, c("utm_x_1000", "utm_y_1000"), mesh = spde$mesh)
fit2 <- sdmTMB(
  n_juv ~ day_night + survey_f + scale_dist + scale_depth + sox_cycle,
  offset = combined_dat$effort,
  data = combined_dat,
  mesh = new_mesh2,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  extra_time = c(2016, 2020),
  silent = FALSE
)


sp_scalar <- (1 * (13 / 1000)) * 500

g <- replicate_df(
  grid_in, "year", time_values = unique(fit2$data$year)
) %>% 
  mutate(
    day_night = "DAY"
  )  %>% 
  left_join(., yr_key, by = "year")
pred2 <- predict(fit2, newdata = g, return_tmb_object = TRUE)
ind <- get_index(pred2, bias_correct = TRUE, area = sp_scalar) 

true <- group_by(sim_grid, year) %>%  
  summarise(true_abund = sum(mu) * sp_scalar) 

ggplot(ind , 
       aes(year, est, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_line(aes(x = year, y = abundance), 
            data = true %>% 
              filter(year > 2021),
            inherit.aes = FALSE,
            colour = "red") +
  ylim(0, NA) +
  ylab("biomass")

ind$cv <- sqrt(exp(ind$se^2) - 1)
  
ind_out <- ind %>% 
  left_join(., true, by = "year") %>% 
  mutate(
    mean_est = mean(est),
    geo_mean_est = exp(mean(log_est)),
    cv = sqrt(exp(se^2) - 1),
    re = (est - true_abund) / true_abund,
    coverage = ifelse(
      lwr < true_abund & upr > true_abund, 1, 0
    )
  ) 

