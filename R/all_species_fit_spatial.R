### Juvenile all species fit 
## Fit spatial model for each species to ensure reasonable convergence and
## to test for effects of survey domain on fixed effect estimates
## Oct 13, 2022

library(tidyverse)
library(sdmTMB)
library(ggplot2)

# downscale data 
dat_full <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    scaled_vol_swept = volume_km3 * 1000,
    effort = log(scaled_vol_swept),
    week = lubridate::week(date)
  ) %>% 
  filter(!species == "steelhead",
         !year_f %in% c("1995",
                      "1996",
                      "1997"),
         lat < 53) %>% 
  droplevels()



# subset data that excludes northern regions to examine potential seasonal 
# effects (once that is excluded also remove 2020 due to small samples)
dat_sub <- dat_full %>% 
  filter(lat < 51.75,
         !year_f == "2020") %>% 
  droplevels() 


## coordinates should be common to all species
dat_coords <- dat_full %>% 
  filter(species == "pink") %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()
dat_coords_sub <- dat_sub %>% 
  filter(species == "pink") %>% 
  select(utm_x_1000, utm_y_1000) %>%
  as.matrix()

## use INLA mesh based on SA recommendations and model selection (see notes)
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(1, 5) * 500,
  cutoff = 20,
  offset = c(20, 200)
)
inla_mesh_raw_sub <- INLA::inla.mesh.2d(
  loc = dat_coords_sub,
  max.edge = c(1, 5) * 500,
  cutoff = 20,
  offset = c(20, 200)
)

spde <- make_mesh(
  dat_full %>% filter(species == "pink"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 
spde_sub <- make_mesh(
  dat_sub %>% filter(species == "pink"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw_sub
) 


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


# switch to tbls
dat_tbl <- dat_full %>% 
  group_by(species) %>% 
  group_nest()
dat_tbl_sub <- dat_sub %>% 
  group_by(species) %>% 
  group_nest()


sp_mod_list1 <- furrr::future_map(
  dat_tbl %>% pull(data),
  ~ {
    sdmTMB(
      n_juv ~ 0 +
        year_f +
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night +
        survey_f,
      offset = .x$effort,
      data = .x,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      anisotropy = TRUE,
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)
summary(sp_mod_list1[[1]])
purrr::map(sp_mod_list1, sanity)


sp_mod_list2 <- furrr::future_map(
  dat_tbl_sub %>% pull(data),
  ~ {
    sdmTMB(
      n_juv ~ 0 +
        year_f +
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night +
        survey_f,
      offset = .x$effort,
      data = .x,
      mesh = spde_sub,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      anisotropy = TRUE,
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        # nlminb_loops = 2
        newton_loops = 1
      ),
      silent = FALSE
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)
summary(sp_mod_list2[[1]])
purrr::map(sp_mod_list2, sanity)


## PREDICTIONS -----------------------------------------------------------------

## fixed effects predicted to be most sensitive, check those first

week_dat <- data.frame(
  week = seq(0, 52, length.out = 100),
  dist_to_coast_km = median(dat_full$dist_to_coast_km),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year_f = "2012"
)
dist_dat <- data.frame(
  week = median(dat_full$week),
  dist_to_coast_km = seq(min(dat_full$dist_to_coast_km), 
                         max(dat_full$dist_to_coast_km), length.out = 100),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year_f = "2012"
)
day_dat <- data.frame(
  week = median(dat_full$week),
  dist_to_coast_km = median(dat_full$dist_to_coast_km),
  day_night = c("DAY", "NIGHT"),
  target_depth = 0,
  survey_f = "hss",
  year_f = "2012"
)
target_dat <- data.frame(
  week = median(dat_full$week),
  dist_to_coast_km = median(dat_full$dist_to_coast_km),
  day_night = "DAY",
  target_depth = seq(min(dat_full$target_depth), 
                     max(dat_full$target_depth), length.out = 100),
  survey_f = "hss",
  year_f = "2012"
)
survey_dat <- data.frame(
  week = median(dat_full$week),
  dist_to_coast_km = median(dat_full$dist_to_coast_km),
  day_night = "DAY",
  target_depth = 0,
  survey_f = c("hss", "ipes"),
  year_f = "2012"
)

pred_tbl <- tibble(
  var = c("week",
          "dist_to_coast_km",
          "day_night",
          "target_depth",
          "survey_f"),
  data = list(week_dat, dist_dat, day_dat, target_dat, survey_dat),
  plot = c("line", "line", "dot", "line", "dot")
)


# predict with each model set
pred_list <- vector(length = nrow(pred_tbl), mode = "list")
for (i in seq_along(pred_tbl$var)) {
  dum <- furrr::future_map2(
    sp_mod_list1, dat_tbl$species, function(x , sp) {
      predict(x, newdata = pred_tbl$data[[i]], se_fit = T, re_form = NA) %>% 
        mutate(species = sp,
               dataset = "full")
    }
  ) %>%
    bind_rows()
  dum_sub <- furrr::future_map2(
    sp_mod_list2, dat_tbl$species, function(x , sp) {
      predict(x, newdata = pred_tbl$data[[i]], se_fit = T, re_form = NA) %>% 
        mutate(species = sp,
               dataset = "sub")
    }
  ) %>%
    bind_rows()
  pred_list[[i]] <- rbind(dum, dum_sub)
}


pred_tbl$pred_dat <- purrr::map(
  pred_list,
  ~ {
    .x %>% 
      group_by(species, dataset) %>% 
      mutate(
        exp_est = exp(est),
        max_est = max(exp_est),
        scale_est = exp_est / max_est,
        up = (est + 1.96 * est_se),
        lo = (est - 1.96 * est_se),
        exp_up = exp(up),
        exp_lo = exp(lo),
        scale_up = exp(up) / max_est,
        scale_lo = exp(lo) / max_est
      )
  }
)

p_dat <- pred_tbl %>% 
  select(var, pred_dat) %>% 
  unnest(cols = pred_dat) 


pp_foo <- function(dat, ...) {
  ggplot(
    dat,
    mapping = aes(!!!ensyms(...))
  ) +
    ggsidekick::theme_sleek() +
    ylab("Abundance Index")
}


pp_foo(
  dat = p_dat %>% filter(var == "week"),
  x = "week", y = "scale_est", ymin = "scale_lo", 
  ymax = "scale_up", fill = "dataset"
) +
  ggsidekick::theme_sleek() +
  ylab("Abundance Index") + 
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, ncol = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Week")

pp_foo(
  dat = p_dat %>% filter(var == "target_depth"),
  x = "target_depth", y = "scale_est", ymin = "scale_lo", 
  ymax = "scale_up", fill = "dataset"
) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, ncol = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Target Headrope Depth (m)") 


pp_foo(
  dat = p_dat %>% filter(var == "dist_to_coast_km"),
  x = "dist_to_coast_km", y = "scale_est", ymin = "scale_lo", 
  ymax = "scale_up", fill = "dataset"
) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, ncol = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Distance to Nearest Coastline (km)") 


pp_foo(
  dat = p_dat %>% filter(var == "day_night"),
  x = "day_night", y = "est", ymin = "lo", 
  ymax = "up", fill = "dataset"
) +
  geom_pointrange(shape = 21,
                  position = position_dodge(0.3)) +
  facet_wrap(~species, nrow = 1) +
  xlab("Diel Effects") +
  theme(
    legend.position = "none"
  )

pp_foo(
  dat = p_dat %>% filter(var == "survey_f"),
  x = "survey_f", y = "est", ymin = "lo", 
  ymax = "up", fill = "dataset"
) +
  geom_pointrange(shape = 21,
                  position = position_dodge(0.3)) +
  facet_wrap(~species, nrow = 1) +
  xlab("Survey Effects") +
  theme(
    legend.position = "none"
  )
