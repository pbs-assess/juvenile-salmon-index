# Simulate and fit unordered dirichlet response models; based on models from 
# compositional package (previous tests indicate equivalent to 
# DirchletReg)

library(dirmult)
library(TMB)
library(tidyverse)
library(Compositional)

set.seed(42)

# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
sd_site <- 0.5
n <- n_sites * n_obs_per_site
n_groups <- 3

#covdat
reg_dat <- data.frame(reg = as.factor(sample(c(1, 2, 3), size = n, 
                                             replace = T)),
                      fac = as.factor(sample(c(1, 2), size = n, replace = T))
)
# model matrix for fixed effects
fix_mm <- model.matrix(~ (reg + fac), reg_dat)

# fixed intercepts (P x K)
P <- ncol(fix_mm)
K <- n_groups - 1 #groups - 1 for reference category
betas <- matrix(rnorm(P * K), P, K)
log_phi <- 3

# function to simulate dirichlet data based on parameters and matrix
f_sim <- function(trial = 1) {
  
  #calculate fixed effects using parameter matrix and model matrix
  n_fix_cov <- P
  fix_eff <- fix_mm %*% betas
  colnames(fix_eff) <- paste("beta_k", seq(2, K+1, by = 1), sep = "")
  
  # generate random effects, then combine with fixed
  site_mean_a <- rnorm(mean = 0,
                       sd = sd_site,
                       n = n_sites)
  site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                       sd = 0,
                       n = n)
  
  total_eff <- fix_eff + site_obs_ia
  
  site_dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1),
                                              each = n_obs_per_site)),
                         site_mean = rep(site_mean_a, each = n_obs_per_site),
                         site_obs = site_obs_ia,
                         sd_site = sd_site
  ) %>%
    cbind(., reg_dat, total_eff)
  
  
  #calc effects
  mu_eff <- cbind(1, exp(total_eff))
  ma <- mu_eff / rowSums(mu_eff)
  ba <- exp(log_phi) * ma

  # Simulate the responses Y from package 'dirmult'
  Y = rdirichlet(n = 1, ba[1, ])
  for(j in c(2:n)){
    Y = rbind(Y, rdirichlet(n = 1, ba[j, ]))
  }
  
  #add small value to any zeros
  Y_adj <- zeroreplace(Y, a = 2/3)  
  
  # create the observations:
  datf <- cbind(site_dat, Y_adj)
  
  return(list("trial" = trial, "obs" = Y, "obs_adj" = Y_adj,
              "fixed_fac" = datf$reg, "full_data" = datf
  ))
}


#prediction model matrix
pred_dat <- reg_dat %>%
  select(reg, fac) %>%
  distinct() %>%
  arrange(reg, fac)
pred_mm <- model.matrix(~ reg + fac, pred_dat)
rand_fac <- as.numeric(sims[[1]]$full_dat$site) - 1

#initial parameter values
beta_in <- matrix(rnorm(n = P * K, mean = 0), P, K)
rfac_in <- rnorm(length(unique(rand_fac)), 0 , 1)


## One simulation

sims <- f_sim(trial = i)

compile(here::here("src", "dirichlet_randInt_v2.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_randInt_v2")))

obj <- MakeADFun(data=list(fx_cov = fix_mm,
                           y_obs = sims[[1]]$obs,
                           pred_cov = pred_mm,
                           rfac = rand_fac),
                 parameters=list(z_ints = beta_in,
                                 log_phi = runif(1, 1, 10),
                                 z_rfac = rfac_in,
                                 log_sigma_rfac = runif(1, 1, 10)),
                 DLL="dirichlet_randInt_v2")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

# check estimates
matrix(ssdr[rownames(ssdr) %in% "z_ints", "Estimate"], ncol = K)
betas
ssdr[rownames(ssdr) %in% "log_phi", "Estimate"]
log_phi
ssdr[rownames(ssdr) %in% "log_sigma_rfac", "Estimate"]
log(sd_site)

matrix(ssdr[rownames(ssdr) %in% "pred_est", "Estimate"], ncol = n_groups)


# SIMULATE FULL DATASET --------------------------------------------------------

n_trials <- 50
sims <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  sims_dat <- f_sim(trial = i)
  sims[[i]] <- list(trial = i, obs = sims_dat$obs,
                    full_dat = sims_dat$full_data, trans = "raw")
}


ests <- map(sims, function(x) {
  y_obs <- x$obs
  
  rand_fac <- as.numeric(x$full_dat$site) - 1
  n_rfac <- length(unique(rand_fac))
  
  fac_dat <- x$full_dat %>% 
    mutate(facs = as.factor(paste(reg, fac, site, sep = "_")),
           facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
    select(reg, fac, site, facs, facs_n)
  fac_key <- fac_dat %>% 
    distinct() %>% 
    arrange(facs_n)
  
  obj <- MakeADFun(data=list(fx_cov = fix_mm,
                             y_obs = y_obs,
                             pred_cov = pred_mm,
                             rfac = rand_fac),
                   parameters=list(z_ints = beta_in,
                                   log_phi = runif(1, 1, 10),
                                   z_rfac = rfac_in,
                                   log_sigma_rfac = runif(1, 1, 10)),
                   DLL="dirichlet_randInt_v2")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  fix_eff <- data.frame(
    par = "beta",
    par_n = seq(1, length(betas), by = 1),
    est = ssdr[rownames(ssdr) %in% "z_ints" , 1],
    true = as.vector(betas)
  )
  rand_eff <- data.frame(
    par = c(rep("alpha", n_rfac), "sigma_a"),
    par_n = c(seq(1, n_rfac, by = 1), 1),
    est = c(ssdr[rownames(ssdr) %in% "z_rfac" , 1], 
            ssdr[rownames(ssdr) %in% "sigma_rfac" , 1]),
    true = c(unique(x$full_dat$site_mean), sd_site)
  )
  
  rbind(fix_eff, rand_eff) %>% 
    mutate(trial = x$trial) 
}) %>% 
  bind_rows()

coef_dat_hier_m <- ests %>%
  mutate(
    sq_err = (true - est)^2
  ) %>% 
  group_by(par, par_n) %>% 
  mutate(rmse = sqrt(mean(sq_err))) %>% 
  ungroup()

par_est_box <- ggplot(coef_dat_hier_m) +
  geom_boxplot(aes(x = as.factor(par_n), y = est)) +
  geom_point(aes(x = as.factor(par_n), y = true), color = "red", shape = 21) +
  ggsidekick::theme_sleek() +
  facet_wrap(~par, scales = "free")

par_rmse_box <- ggplot(coef_dat_hier_m) +
  geom_boxplot(aes(x = as.factor(par_n), y = rmse)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~par, scales = "free")

pdf(here::here("figs", "hier_dirich2_performance.pdf")) 
par_est_box
par_rmse_box
dev.off()