# Simulate and fit unordered dirichlet response models with addition of 
# random intercepts; based on models from following repo:
# https://github.com/carriegu/STAT-840/tree/master/DMRegression

library(dirmult)
library(TMB)
library(tidyverse)
library(DirichletReg)

set.seed(42)

# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
sd_site <- 0.1
# n <- n_sites * n_obs_per_site

# fixed covariate data
J <- 4 #n categories (e.g. stocks)
P <- 3 #n strata (or covariates - 1 if some continuous)
n <- 50 #n observations (e.g. total sampling events)
N <- sample(c(500:1000), n) #sample size per observation

# input data frame and design matrix
dat <- data.frame(strata = sample(1:P, n, replace = T)) %>% 
  mutate(strata_f = as.factor(strata),
         site = sample(1:n_sites, n, replace = T),
         site_f = as.factor(site),
         N = N)
# fixed effects
X <- model.matrix(~strata_f, dat)
beta0 <- matrix(rnorm((ncol(X)) * J), ncol(X), J)
# model matrix for predictions
pred_dat <- dat %>% 
  arrange(strata_f) %>% 
  select(strata_f) %>% 
  distinct()
pred_cov <- model.matrix(~strata_f, pred_dat)

# function to simulate dirichlet data based on parameters and matrix
f_sim <- function(trial = 1) {
  
  # simulate initial beta matrix to generate fixed effects
  fix_eff <- X %*% beta0
  colnames(fix_eff) <- paste("j", seq(1, J, by = 1), sep = "")

  # generate random effects, then combine with fixed
  site_dat <- data.frame(site = seq(1, n_sites, 1),
                         site_mean = rnorm(mean = 0, sd = sd_site, 
                                           n = n_sites))
  dat2 <- dat %>% 
    left_join(., site_dat, by = "site") %>% 
    cbind(., fix_eff)
  
  # reparametrization for data generation
  Gamma = exp(fix_eff + dat2$site_mean) #fixed effects
  Gamma_plus = apply(Gamma, 1, sum) #sum of fixed_effects
  theta = 1/(Gamma_plus + 1)
  pi = apply(Gamma, 2, function(x) {x / theta})
  
  # Simulate the responses Y from package 'dirmult'
  set.seed(123)
  Y = simPop(J = 1, n = N[1], pi = pi[1,], theta = theta[1])$data
  for(jj in c(2:n)){
    Y = rbind(Y, simPop(J = 1, n = N[jj], pi = pi[jj,], theta = theta[jj])$data)
  }
  
  # Switch Y to non-integer
  add_noise <- function(nn) {
    out <- nn
    keep <- which(nn > 1)
    xx <- rnorm(length(nn[keep]), 0, 1)
    noise = xx - mean(xx)
    out[keep] = out[keep] + noise
    return(out)
  }
  Y2 = matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  for (i in seq_len(nrow(Y))) {
    Y2[i, ] <- add_noise(Y[i, ])
  }
  
  # create the observations:
  dat3 <- cbind(dat2, Y)
  
  return(list("trial" = trial, "obs" = Y, "noisey_obs" = Y2, 
              "fixed_fac" = dat3$strata_f, 
              "rand_mean" = site_dat$site_mean, "rand_fac" = dat3$site_f, 
              "full_data" = dat3
  ))
}

# simulate
n_trials <- 40
sims <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  sims[[i]] <- f_sim(trial = i)
}


compile(here::here("src", "dirichlet_randInt.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_randInt")))

#fit models with all data
fit_list_hier <- map(sims, function(sims_in) {
  Y_in <- sims_in$obs #+ runif(length(sims_in$obs), 0, 1)#round(sims_in$obs, 0)
  rfac <- as.numeric(sims_in$rand_fac) - 1 #subtract 1 for indexing in c++
  n_rfac <- length(unique(rfac))
  
  #initial parameter values
  beta_in <- matrix(rnorm((ncol(X)) * J), ncol(X), J)
  rand_int_in <- rep(0, times = length(unique(rfac)))
  
  obj <- MakeADFun(
    data = list(fx_cov = X,
                y_obs = Y_in,
                pred_cov = pred_cov,
                rfac = rfac,
                n_rfac = n_rfac
    ),
    parameters = list(z_ints = beta_in,
                      z_rfac = rand_int_in,
                      log_sigma_rfac = 0
    ),
    random = c("z_rfac"),
    # DLL = "dirichlet_fixInt"
    DLL = "dirichlet_randInt"
  )
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  fix_eff <- data.frame(
    par = "beta",
    par_n = seq(1, length(beta0), by = 1),
    est = ssdr[rownames(ssdr) %in% "z_ints" , 1],
    true = as.vector(beta0)
  )
  rand_eff <- data.frame(
    par = c(rep("alpha", n_rfac), "sigma_a"),
    par_n = c(seq(1, n_rfac, by = 1), 1),
    est = c(ssdr[rownames(ssdr) %in% "z_rfac" , 1], 
            ssdr[rownames(ssdr) %in% "sigma_rfac" , 1]),
    true = c(unique(sims_in$rand_mean), sd_site)
  )
  
  rbind(fix_eff, rand_eff) %>% 
    mutate(tran = sims_in$trans,
           trial = sims_in$trial)  
})

saveRDS(fit_list_hier, here::here("data", "modelFits", "hier_dir_sim_fits.RDS"))
fit_list_hier <- readRDS(here::here("data", "modelFits", "hier_dir_sim_fits.RDS"))



# as above but with noisy data
fit_list_hier_noisey <- map(sims, function(sims_in) {
  Y_in <- sims_in$noisey_obs #+ runif(length(sims_in$obs), 0, 1)#round(sims_in$obs, 0)
  rfac <- as.numeric(sims_in$rand_fac) - 1 #subtract 1 for indexing in c++
  n_rfac <- length(unique(rfac))
  
  #initial parameter values
  beta_in <- matrix(rnorm((ncol(X)) * J), ncol(X), J)
  rand_int_in <- rep(0, times = length(unique(rfac)))
  
  obj <- MakeADFun(
    data = list(fx_cov = X,
                y_obs = Y_in,
                pred_cov = pred_cov,
                rfac = rfac,
                n_rfac = n_rfac
    ),
    parameters = list(z_ints = beta_in,
                      z_rfac = rand_int_in,
                      log_sigma_rfac = 0
    ),
    random = c("z_rfac"),
    # DLL = "dirichlet_fixInt"
    DLL = "dirichlet_randInt"
  )
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  fix_eff <- data.frame(
    par = "beta",
    par_n = seq(1, length(beta0), by = 1),
    est = ssdr[rownames(ssdr) %in% "z_ints" , 1],
    true = as.vector(beta0)
  )
  rand_eff <- data.frame(
    par = c(rep("alpha", n_rfac), "sigma_a"),
    par_n = c(seq(1, n_rfac, by = 1), 1),
    est = c(ssdr[rownames(ssdr) %in% "z_rfac" , 1], 
            ssdr[rownames(ssdr) %in% "sigma_rfac" , 1]),
    true = c(unique(sims_in$rand_mean), sd_site)
  )
  
  rbind(fix_eff, rand_eff) %>% 
    mutate(tran = sims_in$trans,
           trial = sims_in$trial)  
})


# plot coefficient estimates
coef_dat_hier1 <- fit_list_hier %>%
  bind_rows() %>% 
  mutate(
    par_n = case_when(
      par == "sigma_a" ~ 1,
      TRUE ~ par_n),
    sq_err = (true - est)^2
  ) %>% 
  group_by(par, par_n) %>% 
  mutate(rmse = sqrt(mean(sq_err)),
         data_type = "clean") %>% 
  ungroup()
coef_dat_hier2 <- fit_list_hier_noisey %>%
  bind_rows() %>% 
  mutate(
    par_n = case_when(
      par == "sigma_a" ~ 1,
      TRUE ~ par_n),
    sq_err = (true - est)^2
  ) %>% 
  group_by(par, par_n) %>% 
  mutate(rmse = sqrt(mean(sq_err)),
         data_type = "noisy") %>% 
  ungroup()

coef_dat_hier <- rbind(coef_dat_hier1, coef_dat_hier2)

par_est_box <- ggplot(coef_dat_hier) +
  geom_boxplot(aes(x = as.factor(par_n), y = est)) +
  geom_point(aes(x = as.factor(par_n), y = true), color = "red", shape = 21) +
  ggsidekick::theme_sleek() +
  facet_grid(data_type~par, scales = "free")

par_rmse_box <- ggplot(coef_dat_hier) +
  geom_boxplot(aes(x = as.factor(par_n), y = rmse)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~par, scales = "free")

pdf(here::here("figs", "hier_dirich_performance.pdf")) 
par_est_box
par_rmse_box
dev.off()




# look at one sims predictions
# ssdr_out <- fit_list[[1]]$est
# est_beta <- matrix(ssdr_out, nrow = P+1, ncol = K)
# pred_eff <- pred_mm %*% est_beta
# pred_Gamma = exp(pred_eff) #fixed effects
# pred_Gamma_plus = apply(pred_Gamma, 1, sum) #sum of fixed_effects
# pred_theta = 1 / (pred_Gamma_plus + 1)
# pred_pi = apply(pred_Gamma, 2, function(x) {x / pred_theta})
# pred_pi_prop <- pred_pi / rowSums(pred_pi)

ssdr[rownames(ssdr) %in% "pred_pi_prop" , ]
sims_in$full_dat %>%
  pivot_longer(., cols = `1`:`4`, names_to = "cat", values_to = "count") %>% 
  group_by(cat, strata_f) %>%
  summarize(ppn = mean(count) / mean(N)) 
