# Simulate and fit unordered multinomial response models with addition of 
# random intercepts; adjustment of sim_hier_multi_old.
# June 3, 2020

library(tidyverse)
library(TMB)

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

# function to simulate multinomial data and fit hierarchical models
f_sim <- function(trial = 1) {
  
  #calculate fixed effects using parameter matrix and model matrix
  fix_eff <- fix_mm %*% betas
  colnames(fix_eff) <- paste("beta_k", seq(1, K, by = 1), sep = "")
  
  # generate random effects, then combine with fixed
  site_mean_a <- rnorm(mean = 0,
                       sd = sd_site,
                       n = n_sites)
  site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                       sd = 0, 
                       n = n)
  total_eff <- fix_eff + site_obs_ia
  
  datf <- data.frame(site = as.factor(rep(seq(1, 5, by = 1),
                                              each = n_obs_per_site)),
                         site_mean = rep(site_mean_a, each = n_obs_per_site),
                         site_obs = site_obs_ia,
                         sd_site = sd_site) %>%
    cbind(., reg_dat, total_eff) %>% 
    mutate(site = as.factor(site),
           # probability of category 1, 2, and 3:
           p1 = exp(beta_k1) / (1 + exp(beta_k1) + exp(beta_k2)),
           p2 = exp(beta_k2) / (1 + exp(beta_k1) + exp(beta_k2)),
           p3 = 1 - (p1 + p2))
  
  # create the observations:
  y <- numeric(length = nrow(datf))
  p <- datf %>% 
    select(p1:p3) %>%
    as.matrix()
  for (i in seq_along(y)) {
    temp <- rmultinom(1, 1, p[i, ])
    y[i] <- which(temp == 1)
  }
  datf <- cbind(datf, y)

  # matrix version for dmultinom():
  Y <- matrix(ncol = 3, nrow = n, data = 0)
  for (i in seq_along(y)) Y[i, y[i]] <- 1
  
  # true values vector
  tv <- c(as.vector(betas), log(sd_site), site_mean_a, sd_site)
  
  return(list("trial" = trial, "obs" = Y, "fixed_fac" = datf$reg, 
              "rand_fac" = datf$site, "trim_data" = tv, 
              "full_data" = datf))
}

dat_list <- vector(mode = "list", length = 100)
for (i in seq_along(dat_list)) {
  dat_list[[i]] <- f_sim(trial = i)
}


# PREP MODEL INPUTS ------------------------------------------------------------

compile("src/multinomial_randInt_v2.cpp")
dyn.load(dynlib("src/multinomial_randInt_v2"))

## Data and parameters
y_obs <- dat_list[[2]]$obs
#vector of random intercept ids
rfac <- as.numeric(dat_list[[2]]$rand_fac) - 1 #subtract for indexing by 0
fac_dat <- dat_list[[2]]$full_data %>%
  mutate(facs = as.factor(paste(reg, fac, site, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0
  select(reg, fac, site, facs, facs_n)
fac_key <- fac_dat %>%
  distinct() %>%
  arrange(facs_n)

pred_dat <- reg_dat %>%
  select(reg, fac) %>%
  distinct() %>%
  arrange(reg, fac)
pred_mm <- model.matrix(~ reg + fac, pred_dat)

data <- list(y_obs = y_obs,
             rfac = rfac,
             fx_cov = fix_mm, #model matrix from above
             n_rfac = length(unique(rfac)),
             pred_cov = pred_mm
             # ,
             # all_fac = fac_dat$facs_n, # vector of factor combinations
             # fac_key = fac_key$facs_n #ordered unique factor combos in fac_vec
)
parameters <- list(z_rfac = rep(0, times = length(unique(rfac))),
                   z_ints = matrix(0, nrow = ncol(fix_mm),
                                   ncol = ncol(y_obs) - 1),
                   log_sigma_rfac = 0)
obj <- MakeADFun(data, parameters, random = c("z_rfac"),
                 DLL = "multinomial_randInt_v2")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
ssdr <- summary(sdr)


## PLOT PREDICTIONS ------------------------------------------------------------

# make dataframe of true predictions and covariates
dumm <- dat_list[[1]]$full_data 
obs_dat <- dumm %>% 
  pivot_longer(p1:p3, names_to = "cat", names_prefix = "p",
               values_to = "true_prob") %>% 
  select(cat, site, reg, fac, true_prob) %>%
  distinct()

pred_mat <- ssdr[rownames(ssdr) %in% "pred_probs", ]
pred_ci <- data.frame(cat = as.character(rep(1:(K + 1), 
                                             each = nrow(pred_dat))), 
                      pred_prob_est = pred_mat[ , "Estimate"],
                      pred_prob_se =  pred_mat[ , "Std. Error"]) %>% 
  # distinct() %>% 
  mutate(pred_prob_low = pred_prob_est + (qnorm(0.025) * pred_prob_se),
         pred_prob_up = pred_prob_est + (qnorm(0.975) * pred_prob_se),
         reg = rep(pred_dat$reg, times = K + 1),
         fac = rep(pred_dat$fac, times = K + 1)) 

ggplot() +
  geom_boxplot(data = obs_dat,
               aes(x = as.factor(cat), y = true_prob)) +
  geom_pointrange(data = pred_ci, aes(x = as.factor(cat), y = pred_prob_est, 
                                      ymin = pred_prob_low,
                      ymax = pred_prob_up), col = "red") +
  ggsidekick::theme_sleek() +
  facet_grid(reg ~ fac)


## Examine error in coefficient estimates
ests <- map(dat_list, function(x) {
  y_obs <- x$obs
  rfac <- as.numeric(x$rand_fac) - 1 #subtract for indexing by 0
  n_rfac <- length(unique(rfac))
  fac_dat <- x$full_data %>% 
    mutate(facs = as.factor(paste(reg, fac, site, sep = "_")),
           facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
    select(reg, fac, site, facs, facs_n)
  fac_key <- fac_dat %>% 
    distinct() %>% 
    arrange(facs_n)
  
  data <- list(y_obs = y_obs,
               rfac = rfac,
               fx_cov = fix_mm, #model matrix from above
               n_rfac = n_rfac,
               all_fac = fac_dat$facs_n, # vector of factor combinations
               fac_key = fac_key$facs_n #ordered unique factor combos in fac_vec
  )
  parameters <- list(z_rfac = rep(0, times = length(unique(rfac))),
                     z_ints = matrix(0, nrow = ncol(fix_mm),
                                     ncol = k),
                     log_sigma_rfac = 0)
  ## Make a function object
  obj <- MakeADFun(data, parameters, random = c("z_rfac"),
                   DLL = "multinomial_generic_randInt2")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr) 
  
  fix_eff <- data.frame(
    par = "beta",
    par_n = seq(1, length(fix_ints), by = 1),
    est = ssdr[rownames(ssdr) %in% "z_ints" , 1],
    true = as.vector(fix_ints)
  )
  rand_eff <- data.frame(
    par = c(rep("alpha", n_rfac), "sigma_a"),
    par_n = c(seq(1, n_rfac, by = 1), 1),
    est = c(ssdr[rownames(ssdr) %in% "z_rfac" , 1], 
            ssdr[rownames(ssdr) %in% "sigma_rfac" , 1]),
    true = c(unique(x$full_data$site_mean), sd_site)
  )
  
  rbind(fix_eff, rand_eff) %>% 
    mutate(trial = x$trial) 
}) %>% 
  bind_rows()

coef_dat_hier_m <- ests %>%
  bind_rows() %>% 
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

pdf(here::here("figs", "hier_multi_performance.pdf")) 
par_est_box
par_rmse_box
dev.off()
