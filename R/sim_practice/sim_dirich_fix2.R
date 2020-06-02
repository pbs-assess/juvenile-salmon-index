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

#covdat
reg_dat <- data.frame(reg = as.factor(sample(c(1, 2, 3), size = n, 
                                             replace = T)),
                      fac = as.factor(sample(c(1, 2), size = n, replace = T))
)
# model matrix for fixed effects
fix_mm <- model.matrix(~ (reg + fac), reg_dat)

# fixed intercepts (P x K)
P <- ncol(fix_mm) - 1
K <- 4
betas <- matrix(rnorm((P+1)*K), (P+1), K)

# function to simulate dirichlet data based on parameters and matrix
f_sim <- function(trial = 1, scalar = 100) {
  
  #calculate fixed effects using parameter matrix and model matrix
  n_fix_cov <- P
  fix_eff <- fix_mm %*% betas
  sum_fix_eff <- apply(fix_eff, 1, sum)
  
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
    cbind(., reg_dat, sum_fix_eff)  #bind regional fixed effects data
  
  #calc gamma
  Gamma = exp(fix_eff) #fixed effects
  Gamma_plus = apply(Gamma, 1, sum) #sum of fixed_effects
  theta = 1 / (Gamma_plus + 1)
  pi = apply(Gamma, 2, function(x) {x / theta})
  
  #number of values to simulate per row of input
  N <- sample(c(1:1000), n, replace = T)
  
  # Simulate the responses Y from package 'dirmult'
  Y = simPop(J = 1, n = N[1], pi = pi[1,], theta = theta[1])$data
  for(jj in c(2:n)){
    Y = rbind(Y, simPop(J = 1, n = N[jj], pi = pi[jj,], theta = theta[jj])$data)
  }
  
  #switch to proportions and add small value
  Y_adj <- Y / rowSums(Y)  
  Y_adj <- Compositional::zeroreplace(Y_adj, a = 2/3)  
  #scalar is necessary to produce realistic SEs
  Y_adj <- Y_adj * scalar
  
  # create the observations:
  datf <- cbind(site_dat, Y_adj)
  
  return(list("trial" = trial, "obs" = Y, "obs_adj" = Y_adj,
              "fixed_fac" = datf$reg, "full_data" = datf
  ))
}

n_trials <- 5
sims <- vector(mode = "list", length = n_trials)
sims_100 <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  sims_dat <- f_sim(trial = i, scalar = 100)
  sims[[i]] <- list(trial = i, obs = sims_dat$obs,
                    full_dat = sims_dat$full_data, trans = "raw")
  sims_100[[i]] <- list(trial = i, obs = sims_dat$obs_adj,
                        full_dat = sims_dat$full_data, trans = "adj100")
}

sim_list <- c(sims, sims_100)

#initial parameter values
beta_in <- matrix(rnorm(n = (P+1)*K, mean = 0), (P+1), K)

#prediction model matrix
pred_dat <- reg_dat %>%
  select(reg, fac) %>%
  distinct() %>%
  arrange(reg, fac)
pred_mm <- model.matrix(~ reg + fac, pred_dat)

compile(here::here("src", "dirichlet_fixInt_v2.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_fixInt_v2")))


Y_in <- sims_in$obs

obj <- MakeADFun(data=list(fx_cov = fix_mm,
                           y_obs = Y_in,
                           pred_cov = pred_mm),
                 parameters=list(z_ints = beta_in),
                 DLL="dirichlet_fixInt")
opt <- nlminb(obj$par, obj$fn, obj$gr)