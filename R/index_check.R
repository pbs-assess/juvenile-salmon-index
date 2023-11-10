## Check sockeye RW version


library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)


dat_tbl_rw  <- readRDS(
 here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_rw.rds")
) %>% 
  filter(species %in% c("sockeye", "pink"))

dat_tbl  <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_only.rds")
) %>% 
  filter(species %in% c("sockeye", "pink"))

# true index
index_grid_hss <- readRDS(here::here("data", "index_hss_grid.rds")) %>% 
  mutate(day_night = as.factor(day_night),
         trim = ifelse(
           season_f == "wi" & utm_y_1000 < 5551, "yes", "no"
         )) %>%
  #subset to northern domain
  filter(trim == "no") %>% 
  left_join(., yr_key, by = "year")

sp_scalar <- 1 * (13 / 1000)


## fit cyclic models
yr_key <- data.frame(
  year = unique(dat_tbl$data[[2]]$year)
) %>% 
  arrange(year) %>% 
  mutate(
    sox_cycle = rep(1:4, length.out = 25L) %>% 
      as.factor(),
    pink_cycle = rep(1:2, length.out = 25L) %>% 
      as.factor()
  )
sox_dat <- dat_tbl$data[[2]] %>% 
  left_join(., yr_key, by = "year")
sox_fit_cycle <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth + sox_cycle,
  offset = sox_dat$effort,
  data = sox_dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  # time_varying = ~ 1,
  time = "year",
  spatiotemporal = "rw",
  # time_varying_type = "rw0",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(sox_dat$season_f)))
      )
    )
  ),
  silent = FALSE
)

pink_dat <- dat_tbl$data[[1]] %>% 
  left_join(., yr_key, by = "year")
pink_fit_cycle <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth + pink_cycle,
  offset = pink_dat$effort,
  data = pink_dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  # time_varying = ~ 1,
  time = "year",
  spatiotemporal = "rw",
  # time_varying_type = "rw0",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(pink_dat$season_f)))
      )
    )
  ),
  silent = FALSE
)

AIC(pink_fit_cycle, dat_tbl$fit[[1]], dat_tbl_rw$fit[[1]])
AIC(sox_fit_cycle, dat_tbl$fit[[2]], dat_tbl_rw$fit[[2]])



# fit RW indices
cycle_ind_list <- furrr::future_map2(
  dat_tbl_rw$species[4:5], list(pink_fit_cycle, sox_fit_cycle), 
  function(x, y) {
    summer_pred <- predict(
      y,
      newdata = index_grid_hss %>%
        filter(season_f == "su") ,
      se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE
    )
    true_index_summer <- get_index(
      summer_pred, area = sp_scalar, bias_correct = TRUE
    ) %>%
      mutate(season_f = "su",
             species = x)
    fall_pred <- predict(
      y,
      newdata = index_grid_hss %>%
        filter(season_f == "wi") ,
      se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE
    )
    true_index_fall <- get_index(
      fall_pred, area = sp_scalar, bias_correct = TRUE
    ) %>%
      mutate(season_f = "wi",
             species = x)
    list(true_index_summer, true_index_fall)
  }
)

simp_ind_list <- furrr::future_map2(
  dat_tbl$species, dat_tbl$fit, 
  function(x, y) {
    summer_pred <- predict(
      y,
      newdata = index_grid_hss %>%
        filter(season_f == "su") ,
      se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE
    )
    true_index_summer <- get_index(
      summer_pred, area = sp_scalar, bias_correct = TRUE
    ) %>%
      mutate(season_f = "su",
             species = x)
    fall_pred <- predict(
      y,
      newdata = index_grid_hss %>%
        filter(season_f == "wi") ,
      se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE
    )
    true_index_fall <- get_index(
      fall_pred, area = sp_scalar, bias_correct = TRUE
    ) %>%
      mutate(season_f = "wi",
             species = x)
    list(true_index_summer, true_index_fall)
  }
)


rw_index_dat <- purrr::map(rw_ind_list, bind_rows) %>%
  bind_rows() %>% 
  mutate(model = "rw")

simp_index_dat <- purrr::map(simp_ind_list, bind_rows) %>%
  bind_rows() %>% 
  mutate(model = "vanilla")

cyc_index_dat <- purrr::map(cycle_ind_list, bind_rows) %>%
  bind_rows() %>% 
  mutate(model = "cycle")

orig_index_dat <- readRDS(
  here::here("data", "season_index_list_mvrfrw.rds")
) %>% 
  bind_rows() %>% 
  mutate(model = "year_fe")

rbind(rw_index_dat, orig_index_dat) %>% 
  rbind(., cyc_index_dat) %>% 
  filter(species %in% c("pink", "sockeye")) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = year, y = log_est, ymin = (log_est - 1.96*se), 
                      ymax = (log_est + 1.96*se), fill = model),
                  position = position_dodge(width = 1),
                   shape = 21) +
  facet_grid(species~season_f, scales = "free_y") +
  ggsidekick::theme_sleek()



## recover pars
cyc_fit_list <- list(pink_fit_cycle, sox_fit_cycle)
dd <- dat_tbl %>% 
  filter(species %in% c("pink", "sockeye")) %>% 
  mutate(
    cyc_fit = cyc_fit_list
  )

sims_list <- furrr::future_map(
  dd$cyc_fit, function (x) {
    object <- x
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
    simulate(obj_mle, mcmc_samples = sdmTMBextra::extract_mcmc(samp), nsim = 10)
  }
)
dd$sims <- sims_list


sim_tbl <- purrr::pmap(
  list(dd$species, dd$cyc_fit, dd$sims),
  function (sp, fits, x) {
    tibble(
      species = sp,
      iter = seq(1, ncol(x), by = 1),
      sim_dat = lapply(seq_len(ncol(x)), function(i) {
        # use fits instead of data in tibble because mismatched from extra_time
        fits$data %>% 
          mutate(sim_catch = x[ , i]) 
      })
    )
  }
) %>% 
  bind_rows() %>% 
  left_join(., 
            dd %>% select(species, cyc_fit), 
            by = c("species"))

sp_vec <- unique(dd$species)

for (i in seq_along(sp_vec)) {
  sim_tbl_sub <- sim_tbl %>% filter(species == sp_vec[i])
  fit <- furrr::future_map2(
    sim_tbl_sub$sim_dat, sim_tbl_sub$cyc_fit, 
    function (x, fit) {
      if (sp_vec[i] %in% c("pink")) {
        sdmTMB(
          sim_catch ~ 0 + season_f + day_night + survey_f + scale_dist +
            scale_depth + pink_cycle,
          offset = x$effort,
          data = x,
          mesh =  fit$spde,
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
                rep(1, times = length(unique(x$season_f)))
              )
            )
          ),
          silent = FALSE
        )
      } else {
        sdmTMB(
          sim_catch ~ 0 + season_f + day_night + survey_f + scale_dist +
            scale_depth + sox_cycle,
          offset = x$effort,
          data = x,
          mesh =  fit$spde,
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
                rep(1, times = length(unique(x$season_f)))
              )
            )
          ),
          silent = FALSE
        )
      }
    }
  )
  saveRDS(
    fit,
    here::here("data", "fits", "sim_fit", 
               paste(sp_vec[i], "_nb2_mvrfrw_cycle.rds", sep = ""))
  )
}


sim_fit_list <- purrr::map(
  sp_vec,
  ~ readRDS(
    here::here("data", "fits", "sim_fit", paste(.x, "_nb2_mvrfrw_cycle.rds",
                                                sep = ""))
  )
)

sim_tbl$sim_fit <- do.call(c, sim_fit_list)



# extract fixed and ran effects pars from simulations
sim_tbl$pars <- purrr::map(
  sim_tbl$sim_fit, function (x) {
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
)

sim_pars <- sim_tbl %>% 
  select(species, iter, pars) %>% 
  unnest(cols = c(pars)) %>% 
  mutate(
    iter = as.factor(iter),
    species = abbreviate(species, minlength = 3),
    term = fct_recode(
      as.factor(term), 
      "diel" = "day_nightNIGHT", "depth" = "scale_depth",
      "dist" = "scale_dist",
      "spring_int" = "season_fsp", 
      "summer_int" = "season_fsu", "fall_int" = "season_fwi", 
      "survey_design" = "survey_fipes"
    )
  ) 
dir.create(here::here("data", "preds"), showWarnings = FALSE)
saveRDS(sim_pars, here::here("data", "preds", "sim_pars_mvrfrw.rds"))


# as above but for fitted models
fit_effs <- purrr::map2(
  sim_tbl$fit_cyc, sim_tbl$species, function (x, sp) {
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
      mutate(species = abbreviate(sp, minlength = 3)) %>% 
      # add unique identifier for second range term
      group_by(term) %>% 
      mutate(
        par_id = row_number(),
        term = ifelse(par_id > 1, paste(term, par_id, sep = "_"), term)
      ) %>% 
      ungroup()
  }
) %>% 
  bind_rows() %>% 
  mutate(
    term = fct_recode(
      as.factor(term), "diel" = "day_nightNIGHT", "depth" = "scale_depth",
      "dist" = "scale_dist",
      "spring_int" = "season_fsp", 
      "summer_int" = "season_fsu", "fall_int" = "season_fwi", 
      "survey_design" = "survey_fipes"
    )
  )

ggplot() +
  geom_boxplot(data = sim_pars %>%
                 filter(!grepl("year", term)),
               aes(x = species, y = estimate)) +
  geom_point(data = fit_effs %>% 
               filter(!grepl("year", term)),
             aes(x = species, y = estimate), colour = "red") +
  facet_wrap(~term, scales = "free", ncol = 3) +
  labs(y = "Parameter Estimate", x =  "Species") +
  ggsidekick::theme_sleek() 

