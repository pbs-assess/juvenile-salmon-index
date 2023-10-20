## Focus on recovering parameters for coho only

devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "mvrfrw")

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)


ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = 4L)


# make subdirectories for storage
dir.create("data/fits", recursive = TRUE, showWarnings = FALSE)

# downscale data and predictive grid
dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds")) 


# # year season key
# year_season_key <- expand.grid(
#   year = unique(dat_in$year) %>% sort(),
#   season_f = unique(dat_in$season_f) %>% sort()
# ) %>% 
#   arrange(year, season_f) %>% 
#   mutate(
#     # make consecutive index for model fitting
#     ys_index = seq(1, nrow(.), by = 1),
#     year_season_f = paste(year, season_f, sep = "_") %>% as.factor()
#   )


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
  # left_join(., year_season_key %>% select(ys_index, year_season_f)) %>% 
  filter(species == "coho") %>%
  droplevels()


# same specs as normal
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

# current model struc used for other sp
fit <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist + scale_depth,
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
  silent = FALSE,
  do_fit = FALSE
)
fit2 <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
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
fit3 <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist + scale_depth,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  # spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  silent = FALSE
)


saveRDS(fit, here::here("data", "fits", "coho_spatial_varying_nb2_mvrfrw_only.rds"))
saveRDS(fit3, here::here("data", "fits", "coho_spatial_varying_nb2_mvrfrw3.rds"))
fit <- readRDS(here::here("data", "fits", "coho_spatial_varying_nb2_mvrfrw_only.rds"))
fit3 <- readRDS(here::here("data", "fits", "coho_spatial_varying_nb2_mvrfrw3.rds"))


# original model fit with RW through time and extra years
# specify missing season-year levels
# missing_surveys <- year_season_key$ys_index[!(year_season_key$ys_index %in% 
#                                                 dat$ys_index)]

# fit4 <- sdmTMB(
#   n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
#     scale_depth,
#   offset = dat$effort,
#   data = dat,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   spatial = "off",
#   spatial_varying = ~ 0 + season_f + year_f,
#   time_varying = ~ 1,
#   time_varying_type = "rw0",
#   time = "ys_index",
#   spatiotemporal = "off",
#   anisotropy = TRUE,
#   extra_time = missing_surveys,
#   control = sdmTMBcontrol(
#     map = list(
#       # 1 per season, fix all years to same value
#       ln_tau_Z = factor(
#         c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
#       )
#     )
#   ),
#   silent = FALSE
# )


# simulate from year FE
set.seed(456)

# f <- here::here("data", "fits", "nb_mcmc_draws_nb2_mvrfrw_coho.rds")
# if (!file.exists(f)) {
#   # hard coded with large values for cluster
#   object <- fit
#   samp <- sample_mle_mcmc(
#     object, mcmc_iter = 220L, mcmc_warmup = 200L, mcmc_chains = 50L,
#     stan_args = list(thin = 5L, cores = 50L)
#   )
#   # samp <- sample_mle_mcmc(object, mcmc_iter = 105, mcmc_warmup = 100)
# 
#   obj <- object$tmb_obj
#   random <- unique(names(obj$env$par[obj$env$random]))
#   pl <- as.list(object$sd_report, "Estimate")
#   fixed <- !(names(pl) %in% random)
#   map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
#   obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
#   obj_mle <- object
#   obj_mle$tmb_obj <- obj
#   obj_mle$tmb_map <- map
#   sim_out <- simulate(obj_mle, mcmc_samples = sdmTMBextra::extract_mcmc(samp), 
#                       nsim = 200)
# 
#   saveRDS(sim_out, f)
# } else {
#   sim_out <- readRDS(f)
# }

# hard coded with large values for cluster
object <- fit3
samp <- sample_mle_mcmc(object, mcmc_iter = 110, mcmc_warmup = 100)

obj <- object$tmb_obj
random <- unique(names(obj$env$par[obj$env$random]))
pl <- as.list(object$sd_report, "Estimate")
fixed <- !(names(pl) %in% random)
map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
obj_mle <- object
obj_mle$tmb_obj <- obj
obj_mle$tmb_map <- map
sim_out <- simulate(obj_mle, mcmc_samples = sdmTMBextra::extract_mcmc(samp), 
                    nsim = 10)
saveRDS(sim_out, here::here("data", "nb_mcmc_draws_nb2_mvrfrw_coho.rds"))
sim_out <- readRDS(here::here("data", "nb_mcmc_draws_nb2_mvrfrw_coho.rds"))



# fit model to sims
gc()
dir.create(here::here("data", "fits", "sim_fit"), showWarnings = FALSE)


# spread to list then fit 
sims_list <- apply(sim_out, 2, as.list)

f <- here::here("data", "fits", "sim_fit", "coho_nb2_mvrfrw.rds")
if (!file.exists(f)) {
  fit_sims_list <- furrr::future_map(
    sims_list,
    function(x) {
      dum_in <- dat %>% 
        mutate(
          sim_catch = as.numeric(x)
        )
      sdmTMB(
        sim_catch ~ 0 + season_f + day_night + survey_f + scale_dist + scale_depth,
        offset = dum_in$effort,
        data = dum_in,
        mesh = spde,
        family = sdmTMB::nbinom2(),
        spatial = "on",
        # spatial_varying = ~ 0 + season_f,
        time = "year",
        spatiotemporal = "rw",
        anisotropy = TRUE,
        groups = "season_f",
        # control = sdmTMBcontrol(
        #   map = list(
        #     ln_tau_Z = factor(
        #       rep(1, times = length(unique(dum_in$season_f)))
        #     )
        #   )
        # ),
        silent = TRUE
      )
    }
  )
  saveRDS(fit_sims_list, f)
} else {
  fit_sims_list <- readRDS(f)
}


### NOTE!!! Below will only be for coho. Might be easiest for you to generate
## the dataframe of simulated parameter estimates for this model and the
## other species and I can generate the box plot for the supplement.

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

# export (see equivalent snippet in par_recovery_sim.R)
dir.create(here::here("data", "preds"), showWarnings = FALSE)
saveRDS(pars, here::here("data", "preds", "sim_pars_mvrfrw.rds"))


# as above but for fitted model
fix <- tidy(fit2, effects = "fixed")
ran <- tidy(fit2, effects = "ran_pars")

# pull upsilon estimate separately (not currently generated by predict)
est <- as.list(fit2$sd_report, "Estimate", report = TRUE)
se <- as.list(fit2$sd_report, "Std. Error", report = TRUE)
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


## generate indices

index_grid_hss <- readRDS(here::here("data", "index_hss_grid.rds")) %>% 
  mutate(day_night = as.factor(day_night),
         trim = ifelse(
           season_f == "wi" & utm_y_1000 < 5551, "yes", "no"
         )) %>%
  #subset to northern domain
  filter(trim == "no")


# trim_grid <- index_grid_hss %>% 
#   filter(year == "1998" & season_f == "wi" & utm_y_1000 > 5550)
# ggplot(trim_grid) +
#   geom_raster(aes(x = utm_x_1000, y = utm_y_1000, fill = dist_to_coast_km)) +
#   scale_fill_continuous()
# 

sp_scalar <- 1 * (13 / 1000)

index_grid_hss$day_night <- factor(index_grid_hss$day_night, 
                                   levels = levels(dat$day_night))
stopifnot(identical(levels(dat$season_f), levels(index_grid_hss$season_f)))
stopifnot(identical(levels(dat$year_f), levels(index_grid_hss$year_f)))
stopifnot(identical(levels(dat$day_night), levels(index_grid_hss$day_night)))
stopifnot(identical(levels(dat$survey_f), levels(index_grid_hss$survey_f)))


future::plan(future::multisession, workers = 4L)
sim_ind_list_summer <- furrr::future_map(
  fit_sims_list[1:10],
  function (x) {
    pp <- predict(x,
                  newdata = index_grid_hss %>%
                    filter(season_f == "su"),
                  se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
    get_index(pp, area = sp_scalar, bias_correct = TRUE) %>%
      mutate(season_f = "su",
             species = "coho")
  }
)
saveRDS(
  sim_ind_list_summer ,
  here::here("data", "fits", "sim_fit", "coho_summer_sim_index_final_mvrfrw.rds")
)
sim_ind_list_fall <- furrr::future_map(
  fit_sims_list[1:10],
  function (x) {
    pp <- predict(x,
                  newdata = index_grid_hss %>%
                    filter(season_f == "wi"),
                  se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
    get_index(pp, area = sp_scalar, bias_correct = TRUE) %>%
      mutate(season_f = "wi",
             species = "coho")
  }
)
saveRDS(
  sim_ind_list_fall ,
  here::here("data", "fits", "sim_fit", "coho_fall_sim_index_final_mvrfrw.rds")
)


sim_ind_list_summer <- readRDS(
  here::here("data", "fits", "sim_fit", "coho_summer_sim_index_final_mvrfrw.rds")
)
sim_ind_list_fall <- readRDS(
  here::here("data", "fits", "sim_fit", "coho_fall_sim_index_final_mvrfrw.rds")
)


# true indices 
# summer_pred <- predict(
#   fit3, 
#   newdata = index_grid_hss %>%
#     filter(season_f == "su") ,
#   se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE
# )
# true_index_summer <- get_index(
#   summer_pred, area = sp_scalar, bias_correct = TRUE
# ) %>%
#   mutate(season_f = "su",
#          species = "coho")
# fall_pred <- predict(
#   fit3, 
#   newdata = index_grid_hss %>%
#     filter(season_f == "wi") ,
#   se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE
# )
# true_index_fall <- get_index(
#   fall_pred, area = sp_scalar, bias_correct = TRUE
# ) %>%
#   mutate(season_f = "wi",
#          species = "coho")
# 
# true_ind_dat <- rbind(true_index_summer, true_index_fall) %>% 
#   mutate(
#     season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
#     year_f = as.factor(year)
#   )


# combine sims and join
# n_iter <- length(sim_ind_list_summer)
# sim_ind_dat <- tibble(
#   iter = rep(seq(1, n_iter, by = 1), times = 2),
#   sim_index = c(sim_ind_list_summer, sim_ind_list_fall)
# ) %>%
#   unnest(cols = "sim_index") %>%
#   left_join(
#     ., 
#     true_ind_dat %>% 
#       select(species, year, year_f, season_f, true_log_est = log_est),
#     by = c("species", "year", "season_f")) %>% 
#   mutate(
#     season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
#     resid_est = true_log_est - log_est
#   )

n_iter <- length(sim_ind_list_all)
sim_ind_dat <- tibble(
  iter = rep(seq(1, n_iter, by = 1)),
  sim_index = sim_ind_list_all
) %>%
  unnest(cols = "sim_index") %>%
  left_join(
    ., 
    true_ind_dat %>% 
      select(species, year, season_f, true_log_est = log_est),
    by = c("species", "year", "season_f")) %>% 
  mutate(
    season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
    resid_est = true_log_est - log_est
  )


ggplot() +
  geom_boxplot(data = sim_ind_dat,
               aes(x = as.factor(year), y = log_est)) +
  geom_point(data = true_ind_dat,
             aes(x = as.factor(year), y = log_est),
             colour = "red") +
  facet_wrap(~season, scales = "free_y") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top") +
  scale_x_discrete(breaks = seq(2000, 2020, by = 5)) +
  labs(y = "Log Estimated Abundance") +
  theme(axis.title.x = element_blank())


ggplot() +
  geom_boxplot(data = sim_ind_dat,
               aes(x = year_f, y = resid_est)) +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  facet_wrap(~season, scales = "free_y") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top") +
  scale_x_discrete(breaks = seq(2000, 2020, by = 5)) +
  labs(y = "Index Residuals") +
  theme(axis.title.x = element_blank())



pars <- readRDS(here::here("data", "preds", "sim_pars_mvrfrw.rds"))


p_foo <- function(x, dat_in) {
  pp <- predict(x,
                newdata = dat_in,
                return_tmb_report = TRUE)
  p_out <- predict(x,
                   newdata = dat_in,
                   se_fit = FALSE, re_form = NULL)
  
  p_out %>%
    mutate(upsilon_stc = as.numeric(pp$proj_upsilon_st_A_vec))
}   


## compare spatial predictions
fall_pred_sim <- p_foo(
  x = fit_sims_list[[1]], 
  dat_in = index_grid_hss %>%
    filter(season_f == "wi",
           utm_y_1000 > 5500)
)
fall_pred_real <- p_foo(
  x = fit, 
  dat_in = index_grid_hss %>%
    filter(season_f == "wi",
           utm_y_1000 > 5500)
)

fall_pred <- fall_pred_sim %>% 
  mutate(sim_minus_real = est - fall_pred_real$est,
         sim_minus_real_eps = upsilon_stc - fall_pred_real$upsilon_stc,
         sim_minus_real_fe = est_non_rf - fall_pred_real$est_non_rf)

ggplot(fall_pred %>% filter(year > 2010),
       aes(x = utm_x_1000, y = utm_y_1000, fill = sim_minus_real_eps)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~ year)


summer_pred_sim <- p_foo(
  x = fit_sims_list[[1]], 
  dat_in = index_grid_hss %>%
    filter(season_f == "su",
           utm_y_1000 > 5550)
)
summer_pred_real <- p_foo(
  x = fit, 
  dat_in = index_grid_hss %>%
    filter(season_f == "su",
           utm_y_1000 > 5500)
)

summer_pred <- summer_pred_sim %>% 
  mutate(sim_minus_real = est - summer_pred_real$est,
         sim_minus_real_eps = upsilon_stc - summer_pred_real$upsilon_stc)

ggplot(summer_pred %>% filter(year > 2010),
       aes(x = utm_x_1000, y = utm_y_1000, fill = sim_minus_real_eps)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~ year)
