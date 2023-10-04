## Focus on recovering parameters for coho only


library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)

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
    day_night = as.factor(day_night)) %>% 
  filter(species == "coho") %>% 
  droplevels()


## mesh shared among species
dat_coords <- dat %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(2, 10) * 500,
  cutoff = 30,
  offset = c(10, 50)
)  
plot(inla_mesh_raw)
spde <- make_mesh(
  dat,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 
plot(spde)

# mvrfrw only
dat_tbl_mvrfrw <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_only.rds")
)
fit_nb2 <- dat_tbl_mvrfrw$fit[[3]]

fit_nb2_poly <- sdmTMB(
  n_juv ~ 0 + survey_f + day_night + poly(scale_yday, 2) + scale_dist +
    scale_depth,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat$season_f)))
      )
    )
  ),
  silent = FALSE
)

fit_nb2_poly <- sdmTMB(
  n_juv ~ 0 + survey_f + day_night + poly(scale_yday, 2) + scale_dist +
    scale_depth,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  share_range = FALSE,
  # control = sdmTMBcontrol(
  #   map = list(
  #     ln_tau_Z = factor(
  #       rep(1, times = length(unique(dat$season_f)))
  #     )
  #   )
  # ),
  silent = FALSE
)

dat2 <- dat %>% 
  group_by(season_f) %>% 
  mutate(
    scale_yday = scale(as.numeric(yday))[ , 1],
  ) %>% 
  ungroup()
fit_nb2_yday <- sdmTMB(
  n_juv ~ 0 + survey_f + day_night + season_f + scale_yday:season_f + 
    scale_dist + scale_depth,
  offset = dat2$effort,
  data = dat2,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat2$season_f)))
      )
    )
  ),
  silent = FALSE
)
saveRDS(fit_nb2_yday, here::here("data", "fits", "coho_mvrfrw_yday_int.rds"))

fit_nb2_yday2 <- sdmTMB(
  n_juv ~ 0 + survey_f + day_night + #season_f + 
    scale_yday:season_f + 
    scale_dist + scale_depth,
  offset = dat2$effort,
  data = dat2,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat2$season_f)))
      )
    )
  ),
  silent = FALSE
)
saveRDS(fit_nb2_yday2, here::here("data", "fits", "coho_mvrfrw_yday_int2.rds"))

fit_nb2_yday3 <- sdmTMB(
  n_juv ~ survey_f + day_night + season_f + scale_yday:season_f + 
    scale_dist + scale_depth,
  offset = dat2$effort,
  data = dat2,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  silent = FALSE
)
saveRDS(fit_nb2_yday3, here::here("data", "fits", "coho_mvrfrw_yday_int3.rds"))

dat3 <- dat2 %>% filter(!season_f == "sp")
spde2 <- make_mesh(
  dat3,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 
fit_nb2_yday4 <- sdmTMB(
  n_juv ~ survey_f + day_night + season_f + scale_yday:season_f + 
    scale_dist + scale_depth,
  offset = dat3$effort,
  data = dat3,
  mesh = spde2,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  silent = FALSE
)
saveRDS(fit_nb2_yday4, here::here("data", "fits", "coho_mvrfrw_yday_int4.rds"))



fit_nb2_smooth <- sdmTMB(
  n_juv ~ 0 + survey_f + day_night + s(scale_yday, k = 4) + scale_dist +
    scale_depth,
  # knots = list(yday = c(0, 365)),
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
  # control = sdmTMBcontrol(
  #   map = list(
  #     ln_tau_Z = factor(
  #       rep(1, times = length(unique(dat$season_f)))
  #     )
  #   )
  # ),
  silent = FALSE
)


# compare to original fits
dat_tbl_mvrfrw <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_only.rds")
)
fit_nb2 <- dat_tbl_mvrfrw$fit[[3]]
dat_tbl_fe <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw.rds")
) 
fit_nb2_fe <- dat_tbl_fe$fit[[3]]


# plot residuals 
dat$resid1 <- resid(fit_nb2)
dat$resid2 <- resid(fit_nb2_yday)

ggplot(dat) +
  geom_boxplot(aes(x = as.factor(year), y = resid2)) +
  facet_wrap(~season_f)
plot(dat$resid1 ~ dat$resid2)


### Compare fits to original and simulate data 
sims_list <- readRDS(here::here("data", "fits", "nb_mcmc_draws_nb2_mvrfrw.rds"))
coho_sims <- sims_list[[3]]
dum_coho <- dat %>% 
  mutate(sim_catch = coho_sims[ , 1])

fit_sim <-  sdmTMB(
  sim_catch ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth,
  offset = dum_coho$effort,
  data = dum_coho,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat$season_f)))
      )
    )
  ),
  silent = FALSE
)



### Spatial predictions for original and simulated fit

# predictive grid
grid_list <- readRDS(here::here("data", "spatial", "pred_ipes_grid.RDS")) %>% 
  purrr::map(
    .,
    ~ {.x %>% 
        mutate(utm_x_1000 = X / 1000,
               utm_y_1000 = Y / 1000,
               dist_to_coast_km = shore_dist / 1000)}
  )

year_season_key <- expand.grid(
  year = unique(dat$year) %>% sort(),
  season_f = unique(dat$season_f) %>% sort()
) %>% 
  arrange(year, season_f) 

summer_grid <- grid_list$ipes_grid
fall_grid <- grid_list$wcvi_grid

pred_grid_list <- rbind(
  year_season_key %>% mutate(survey_f = "hss"),
  year_season_key %>% mutate(survey_f = "ipes")
) %>% 
  mutate(
    year_f = as.factor(year),
    id = row_number(),
    survey_f = survey_f %>% as.factor()
  ) %>% 
  split(., .$id)
spatial_grid <- pred_grid_list %>%
  purrr::map(., function (x) {
    fall_grid %>%
      mutate(
        year = x$year,
        year_f = x$year_f,
        survey_f = x$survey_f,
        season_f = x$season_f,
        ys_index = x$ys_index,
        target_depth = 0,
        day_night = "DAY",
        scale_depth = (target_depth - mean(dat$target_depth)) /
          sd(dat$target_depth),
        scale_dist = (dist_to_coast_km  - mean(dat$dist_to_coast_km )) /
          sd(dat$dist_to_coast_km ),
        scale_yday = 0
      )
  }) %>%
  bind_rows() %>%
  filter(
    !survey_f == "ipes"
  ) %>%
  select(-c(depth, slope))


fit_list <- list(fit_nb2_yday3, fit_sims_list[[1]])
nms <- c("original", "sim")


spatial_preds_list <- furrr::future_map2(
  fit_list,
  nms,
  ~ {
    # separate predictions to generate upsilon estimates
    pp <- predict(.x,
                  newdata = spatial_grid,
                  return_tmb_report = TRUE)
    p_out <- predict(.x,
                     newdata = spatial_grid,
                     se_fit = FALSE, re_form = NULL)

    p_out %>%
      mutate(dataset = .y,
             upsilon_stc = as.numeric(pp$proj_upsilon_st_A_vec))
  }
)

spatial_preds <- spatial_preds_list %>%
  bind_rows() %>%
  filter(!season_f == "sp")


year_seq <- sample(seq(2012, 2020, by = 1), 6)
sp_scalar <- 1 * (13 / 1000)
sub_spatial <- spatial_preds %>% 
  filter(year %in% year_seq) %>%
  group_by(dataset) %>% 
  mutate(
    grid_est = sp_scalar * exp(est),
    scale_est = grid_est / max(grid_est)
  ) %>% 
  ungroup() 

# check summed abundance (similar patterns to index)
spatial_preds %>% 
  group_by(year_f, season_f, dataset) %>% 
  summarize(
    sum_est = sum(est)
  ) %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_point(aes(x = year_f, y = sum_est, colour = dataset)) +
  facet_wrap(~season_f)


scale_max <- quantile(sub_spatial$scale_est, 0.99)
sub_spatial$scale_est2 <- ifelse(
  sub_spatial$scale_est > scale_max, 
  scale_max, 
  sub_spatial$scale_est
)

plot_map_raster <- function(dat, column = est) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c()  +
    ggsidekick::theme_sleek() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "top") +
    scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
    scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0))
}

summer_p <- plot_map_raster(dat = sub_spatial %>% 
                  filter(season_f == "su",
                         year %in% year_seq),
                column = est) +
  facet_grid(dataset~year) +
  labs(title = "Summer Spatial Preds (log space)")

fall_p <- plot_map_raster(dat = sub_spatial %>% 
                              filter(season_f == "wi",
                                     year %in% year_seq),
                            column = est) +
  facet_grid(dataset~year) +
  labs(title = "Fall Spatial Preds (log space)")

# fall difference in estimate
sub_spatial %>% 
  filter(season_f == "wi",
         year %in% year_seq) %>% 
  select(X, Y, year_f, est, dataset) %>% 
  pivot_wider(names_from = dataset, values_from = est) %>% 
  mutate(diff_est = sim - original) %>% 
  plot_map_raster(dat = .,
                  column = diff_est) +
  facet_wrap(~year_f) +
  scale_fill_gradient2(name = "Spatial\nField")


summer_epsilon <- plot_map_raster(dat = sub_spatial %>% 
                                  filter(season_f == "su",
                                         year %in% year_seq),
                                column = upsilon_stc) +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(dataset~year) +
  labs(title = "Summer Epsilon")
fall_epsilon <- plot_map_raster(dat = sub_spatial %>% 
                                          filter(season_f == "wi",
                                                 year %in% year_seq),
                                        column = upsilon_stc) +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(dataset~year) +
  labs(title = "Fall Epsilon")

omega <- plot_map_raster(dat = sub_spatial %>% 
                           filter(season_f == "wi",
                                  year %in% year_seq),
                         column = omega_s) +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(~dataset) +
  labs(title = "Fall Epsilon")

# fall difference in epsilon
sub_spatial %>% 
  filter(season_f == "wi",
         year %in% year_seq) %>% 
  select(X, Y, year_f, upsilon_stc, dataset) %>% 
  pivot_wider(names_from = dataset, values_from = upsilon_stc) %>% 
  mutate(diff_epsilon = sim - original) %>% 
  plot_map_raster(dat = .,
                  column = diff_epsilon) +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_wrap(~year_f)



# spatial RF
omega_season <- sub_spatial %>% 
  # values will be duplicated across years, seasons and surveys; select one
  filter(year_f == year_seq[1], season_f == "wi") %>% 
  select(X, Y, dataset, starts_with("zeta_s_season_")) %>%
  pivot_longer(cols = starts_with("zeta"), values_to = "omega_est", 
               names_to = "season", names_prefix = "zeta_s_season_f") %>% 
  filter(!season == "sp") %>% 
  mutate(season = fct_recode(season, "summer" = "su", "fall" = "wi"))
plot_map_raster(dat = omega_season,
                column = omega_est) +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(dataset~season)  +
  labs(title = "Omega")


# fixed parameter estimates
fix_pars <- purrr::pmap(
  list(fit_list, nms),
  function(x, y) {
    tidy(x, effects = "fixed", conf.int = T) %>% 
      mutate(
        dataset = y
      ) 
  }
) %>% 
  bind_rows() %>% 
  filter(term %in% c("season_fsp", "season_fsu", "season_fwi", "survey_fipes", 
                     "day_nightNIGHT")) %>% 
  mutate(
    # term2 = ifelse(grepl("season", term), "Season", term),
    term = fct_recode(as.factor(term), "Spring" = "season_fsp", 
                      "Summer" = "season_fsu", "Fall" = "season_fwi",
                      "IPES Survey" = "survey_fipes", 
                      "Nocturnal Sampling" = "day_nightNIGHT")
  )
ran_pars <- purrr::pmap(
  list(fit_list, nms),
  function(x, y) {
    pars1 <- tidy(x, effects = "ran_par", conf.int = T) 
    
    # pull upsilon estimate separately (not currently generated by predict)
    est <- as.list(x$sd_report, "Estimate", report = TRUE)
    se <- as.list(x$sd_report, "Std. Error", report = TRUE)
    upsilon <- data.frame(
      term = "sigma_upsilon",
      log_est = est$log_sigma_U,
      log_se = se$log_sigma_U
    ) %>% 
      mutate(
        estimate = exp(log_est),
        std.error = exp(log_se),
        conf.low = exp(log_est + (stats::qnorm(0.025) * log_se)),
        conf.high = exp(log_est + (stats::qnorm(0.975) * log_se)) 
      ) %>% 
      select(term, estimate, std.error, conf.low, conf.high)
    
    rbind(pars1, upsilon) %>% 
      mutate(
        dataset = y
      )
  }
) %>% 
  bind_rows()
ran_pars$term <- fct_recode(ran_pars$term, "sigma_omega" = "sigma_Z")

rbind(fix_pars, ran_pars) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = dataset, y = estimate, ymin = conf.low, 
                      ymax = conf.high)) +
  facet_wrap(~term, scales = "free_y") +
  ggsidekick::theme_sleek()



# simulate from year FE
set.seed(456)
object <- fit_nb2_yday4
samp <- sample_mle_mcmc(object, mcmc_iter = 120, mcmc_warmup = 100)

obj <- object$tmb_obj
random <- unique(names(obj$env$par[obj$env$random]))
pl <- as.list(object$sd_report, "Estimate")
fixed <- !(names(pl) %in% random)
map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
obj_mle <- object
obj_mle$tmb_obj <- obj
obj_mle$tmb_map <- map
sims_fe <- simulate(obj_mle, mcmc_samples = sdmTMBextra::extract_mcmc(samp), 
                    nsim = 10)


# spread to list then fit 
sims_list <- apply(sims_fe, 2, as.list)
fit_sims_list <- furrr::future_map(
  sims_list,
  function(x) {
    dum_in <- dat3 %>% 
      mutate(
        sim_catch = as.numeric(x)
      )
    sdmTMB(
      sim_catch ~ survey_f + day_night + season_f + 
        scale_yday:season_f + 
        scale_dist + scale_depth,
      offset = dum_in$effort,
      data = dum_in,
      mesh = spde2,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      time = "year",
      spatiotemporal = "rw",
      anisotropy = TRUE,
      groups = "season_f",
      silent = FALSE
    )
  }
)


# check parameter estimates
pars <- purrr::map(
  fit_sims_list, function (x) {
    fix <- tidy(x, effects = "fixed") 
    ran <- tidy(x, effects = "ran_pars") 
    
    # pull upsilon estimate separately (not currently generated by predict)
    est <- as.list(x$sd_report, "Estimate", report = TRUE)
    se <- as.list(x$sd_report, "Std. Error", report = TRUE)
    upsilon <- data.frame(
      term = "sigma_U",
      log_est = est$log_sigma_U,
      log_se = se$log_sigma_U
    ) %>% 
      mutate(
        estimate = exp(log_est),
        std.error = exp(log_se) 
      ) %>% 
      select(term, estimate, std.error)
    
    rbind(fix, ran, upsilon) %>% 
      # add unique identifier for second range term
      group_by(term) %>% 
      mutate(
        par_id = row_number(),
        term = ifelse(par_id > 1, paste(term, par_id, sep = "_"), term)
      ) %>% 
      ungroup() 
  }
) %>% 
  bind_rows()


# as above but for fitted model
fix <- tidy(fit_nb2_yday4, effects = "fixed")
ran <- tidy(fit_nb2_yday4, effects = "ran_pars")

# pull upsilon estimate separately (not currently generated by predict)
est <- as.list(fit_nb2_yday4$sd_report, "Estimate", report = TRUE)
se <- as.list(fit_nb2_yday4$sd_report, "Std. Error", report = TRUE)
upsilon <- data.frame(
  term = "sigma_U",
  log_est = est$log_sigma_U,
  log_se = se$log_sigma_U
) %>% 
  mutate(
    estimate = exp(log_est),
    std.error = exp(log_se) 
  ) %>% 
  select(term, estimate, std.error)

orig_pars <- rbind(fix, ran, upsilon) %>% 
  # add unique identifier for second range term
  group_by(term) %>% 
  mutate(
    par_id = row_number(),
    term = ifelse(par_id > 1, paste(term, par_id, sep = "_"), term)
  ) %>% 
  ungroup()

ggplot() +
  geom_boxplot(data = pars %>%
                 filter(!grepl("year", term)),
               aes(x = term, y = estimate)) +
  geom_point(data = orig_pars %>% 
               filter(!grepl("year", term)),
             aes(x = term, y = estimate), colour = "red") +
  facet_wrap(~term, scales = "free", ncol = 3) +
  labs(y = "Parameter Estimate", x =  "Species") +
  ggsidekick::theme_sleek() 



# check fall index

index_grid_hss <- readRDS(here::here("data", "index_hss_grid.rds"))
index_grid_hss$scale_yday <- 0
sp_scalar <- 1 * (13 / 1000)

sim_ind_list_fall <- furrr::future_map(
  fit_sims_list,
  function (x) {
    pp <- predict(x,
                  newdata = index_grid_hss %>%
                    filter(season_f == "wi") %>% 
                    droplevels(),
                  se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
    get_index(pp, area = sp_scalar, bias_correct = TRUE) %>%
      mutate(season_f = "wi")
  }
)

true_pred <- predict(fit_nb2_yday4,
                     newdata = index_grid_hss %>%
                       filter(season_f == "wi") %>% 
                       droplevels(),
                     se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
true_ind_dat <- get_index(true_pred, area = sp_scalar, bias_correct = TRUE) %>%
  mutate(season_f = "wi")
true_ind_dat$year_f <- as.factor(true_ind_dat$year)



# sim_ind_list_summer <- readRDS(
#   here::here("data", "fits", "sim_fit",
#              paste("coho", "_summer_sim_index_final_mvrfrw.rds", sep = ""))
# )
# sim_ind_list_fall <- readRDS(
#     here::here("data", "fits", "sim_fit",
#                paste("coho", "_fall_sim_index_final_mvrfrw.rds", sep = ""))
#   )
# 
# 
# true_index_list <- readRDS(
#   here::here("data", "fits", "season_index_list_mvrfrw.rds")
# )
# 
# true_ind_dat <- true_index_list %>%
#   bind_rows() %>% 
#   mutate(
#     season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
#     # survey = case_when(
#     #   season == "fall" & year %in% fall_years ~ "sampled",
#     #   season == "summer" & year %in% summer_years ~ "sampled",
#     #   TRUE ~ "no survey"
#     # ),
#     log_lwr = log_est - (1.96 * se),
#     log_upr = log_est + (1.96 * se),
#     year_f = as.factor(year)
#   ) %>% 
#   filter(species == "coho")

n_iter <- length(sim_ind_list_fall) #how many sims 

sample_day <- dat %>% 
  group_by(year, season_f) %>% 
  summarize(mean_yday = mean(yday, na.rm = TRUE))

sim_ind_dat <- rbind(
  tibble(
    iter = seq(1, n_iter, by = 1),
    sim_index = sim_ind_list_fall
  )#,
  # tibble(
  #   iter = seq(1, n_iter, by = 1),
  #   sim_index = sim_ind_list_summer
  # )
) %>%
  unnest(cols = "sim_index") %>%
  left_join(
    ., 
    true_ind_dat %>% 
      select(year, year_f, season_f, true_log_est = log_est),
    by = c("year", "season_f")) %>% 
  # left_join(., sample_day, by = c("year", "season_f")) %>% 
  mutate(
    # season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
    resid_est = true_log_est - log_est
    # resid_est = est - exp(true_log_est)
  )

ggplot() +
  geom_boxplot(data = sim_ind_dat,
               aes(x = year_f, y = log_est)) +
  geom_point(data = true_ind_dat,
             aes(x = year_f, y = log_est), colour = "red") +#, colour = survey)) +
  # facet_grid(~season, scales = "free_y") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top") +
  scale_x_discrete(breaks = seq(2000, 2020, by = 5)) +
  labs(y = "Log Estimated Abundance") +
  theme(axis.title.x = element_blank())

ggplot() +
  geom_boxplot(data = sim_ind_dat,
               aes(x = year_f, y = resid_est)) +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  # facet_wrap(~season_f, scales = "free_y") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top") +
  scale_x_discrete(breaks = seq(2000, 2020, by = 5)) +
  labs(y = "Index Residuals") +
  theme(axis.title.x = element_blank())

ggplot(sample_day %>% filter(season_f == "wi")) +
  geom_point(aes(x = as.factor(year), y = mean_yday))

sim_ind_dat %>% 
  filter(!is.na(mean_yday)) %>% 
  group_by(year_f, mean_yday, season_f) %>% 
  summarize(mean_resid = median(resid_est)) %>% 
  ggplot(data = .) +
  geom_point(aes(x = mean_yday, y = mean_resid, fill = year_f), shape = 21) +
  ggsidekick::theme_sleek() +
  facet_wrap(~season_f, scales = "free_x")



## compare spatial distribution of 2012-2022 fall observed samples to sim


dum_in <- dat2 %>% 
  mutate(
    sim_catch = as.numeric(sims_list[[1]])
  ) %>% 
  pivot_longer(cols = c(sim_catch, n_juv)) 

ggplot(dum_in %>% filter(year > 2012, season_f == "wi")) +
  geom_point(aes(x = utm_x_1000, y = utm_y_1000, size = value)) +
  facet_grid(name ~ year_f) #+
#  scale_size_continuous(trans = "log")
