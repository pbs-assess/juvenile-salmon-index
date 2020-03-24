# Simulate data and fit tweedie model

library(tidyverse)
library(TMB)

set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 250
sd_site <- 0.5
N <- n_sites * n_obs_per_site

fix_dat <- data.frame(reg = as.factor(sample(c(1, 2, 3), size = N, replace = T)),
                      seas = as.factor(sample(c(1, 2), size = N, replace = T)))

# fixed intercepts
fix_ints <- c(0.75, 1, 2, -0.25) 

# model matrix for fixed effects
fix_mm <- model.matrix(~ reg + seas, fix_dat)


# function to simulate hier tweedie data and fit hierarchical models
f_sim <- function(trial = 1) {
  #calculate fixed effects using parameter matrix and model matrix
  sum_fix_eff <- fix_mm %*% fix_ints
  
  # generate random effects, then combine with fixed
  site_mean_a <- rnorm(mean = 0,
                       sd = sd_site,
                       n = n_sites)
  mu_vec <- as.numeric(exp(rep(site_mean_a, each = n_obs_per_site) + 
    sum_fix_eff))
  data.frame(trial = rep(trial, length.out = N),
             site = as.factor(rep(seq(1, 5, by = 1), each = n_obs_per_site)),
             site_mean = rep(site_mean_a, each = n_obs_per_site), 
             sd_site = sd_site
             ) %>%
    cbind(fix_dat) %>%  #bind regional fixed effects data 
    mutate(reg_eff = sum_fix_eff,
           mu = mu_vec,
           site_obs = tweedie::rtweedie(N, mu = mu_vec, power=2, phi=1)) 
}

dat <- f_sim()
hist(dat$site_obs)
ggplot(dat) +
  geom_boxplot(aes(x = reg, y = site_obs)) +
  facet_wrap(~seas)

m2 <- glmmTMB::glmmTMB(site_obs ~ reg + seas, data = dat, 
                       family = glmmTMB::tweedie(link = "log")
                       )
summary(m2)
m1 <- lm(log(site_obs) ~ reg + seas, data = dat)
summary(m1)

#helper function to convert factors 
fct_to_tmb_num <- function(x) {
  as.numeric(as.factor(as.character(x))) - 1
}

fac1k <- fct_to_tmb_num(dat$site)
# empty model matrix for predictions
# mm_pred <- fix_mm[1:(ncol(fix_mm) + 1), ]
# for (i in 1:ncol(mm_pred)) {
#   for (j in 1:nrow(mm_pred)) {
#     mm_pred[j, i] <- 0
#   }}
# mm_pred[,1] <- 1
# for (i in 1:ncol(mm_pred)) {
#   for (j in 1:nrow(mm_pred)) {
#     if (i == j)
#       mm_pred[j, i] <- 1
#   }}
# mm_pred[nrow(mm_pred), ] <- 1

# order matrix based on unique factor levels with most saturated at bottom
mm_pred <- fix_mm %>% 
  unique() 
ord_mat <- mm_pred %>% 
  apply(., 1, sum) %>%
  sort() %>% 
  names()
mm_pred <- mm_pred[ord_mat, ]

data <- list(y1_i = dat$site_obs,
             X1_ij = fix_mm,
             factor1k_i = fac1k,
             nk1 = length(unique(fac1k)),
             X1_pred_ij = mm_pred
             )
parameters = list(
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),#rep(1, ncol(fix_mm)),
  log_phi = log(1.1),
  logit_p = boot::logit(0.8),
  z1_k = rep(0, length(unique(fac1k))),
  log_sigma_zk1 = log(0.25)
  )

## Make a function object
# compile("R/sim_practice/tweediePractice/tweedie_cpue_fe.cpp")
# dyn.load(dynlib("R/sim_practice/tweediePractice/tweedie_cpue_fe"))
# obj <- MakeADFun(data, parameters, DLL = "tweedie_cpue_fe")

compile("R/sim_practice/tweediePractice/tweedie_cpue_1re.cpp")
dyn.load(dynlib("R/sim_practice/tweediePractice/tweedie_cpue_1re"))
obj <- MakeADFun(data, parameters, random = c("z1_k"), 
                 DLL = "tweedie_cpue_1re")

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr

# helper function for inverse logit of p
inv_logit <- function(x) {
  1.0 / (1.0 + exp(-x)) + 1.0
}
