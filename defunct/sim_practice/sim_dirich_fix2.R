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
  site_dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1),
                                              each = n_obs_per_site)),
                         site_mean = rep(site_mean_a, each = n_obs_per_site),
                         site_obs = site_obs_ia,
                         sd_site = sd_site
  ) %>%
    cbind(., reg_dat, fix_eff)  #bind regional fixed effects data
  
  #calc effects
  mu_eff <- cbind(1, exp(fix_eff))
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

n_trials <- 5
sims <- vector(mode = "list", length = n_trials)
# sims_100 <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  sims_dat <- f_sim(trial = i)
  sims[[i]] <- list(trial = i, obs = sims_dat$obs,
                    full_dat = sims_dat$full_data, trans = "raw")
  # sims_100[[i]] <- list(trial = i, obs = sims_dat$obs_adj,
  #                       full_dat = sims_dat$full_data, trans = "adj100")
}

# sim_list <- c(sims, sims_100)

#initial parameter values
beta_in <- matrix(rnorm(n = P * K, mean = 0), P, K)

#prediction model matrix
pred_dat <- reg_dat %>%
  select(reg, fac) %>%
  distinct() %>%
  arrange(reg, fac)
pred_mm <- model.matrix(~ reg + fac, pred_dat)

# fit compositional package model
temp_x <- sims[[1]]$full_dat %>% 
  select(reg, fac)
mod1 <- Compositional::diri.reg(y = sims[[1]]$obs, x = temp_x,
                                xnew = pred_dat)

compile(here::here("src", "dirichlet_fixInt_v2.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_fixInt_v2")))

obj <- MakeADFun(data=list(fx_cov = fix_mm,
                           y_obs = sims[[1]]$obs,
                           pred_cov = pred_mm),
                 parameters=list(z_ints = beta_in,
                                 log_phi = runif(1, 1, 10)),
                 DLL="dirichlet_fixInt_v2")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

# check estimates
matrix(ssdr[rownames(ssdr) %in% "z_ints", "Estimate"], ncol = K)
mod1$be
betas
ssdr[rownames(ssdr) %in% "log_phi", "Estimate"]
mod1$log.phi
matrix(ssdr[rownames(ssdr) %in% "pred_est", "Estimate"], ncol = n_groups)
mod1$est
