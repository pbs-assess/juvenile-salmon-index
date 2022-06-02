# Simulate data and fit negbinomial model

library(tidyverse)
library(TMB)

set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 250
sd_site <- 0.5
N <- n_sites * n_obs_per_site

fix_dat <- data.frame(eff_z = rnorm(N, 0, 1), #continuous pred
                      reg = as.factor(sample(c(1, 2, 3), size = N, replace = T)),
                      seas = as.factor(sample(c(1, 2), size = N, replace = T)))

# fixed intercepts and slope parameter (last)
fix_pars <- c(0.75, 1, 2, -0.25, 0.5) 

# model matrix for fixed effects
fix_mm <- model.matrix(~ reg + seas + eff_z, fix_dat)


# function to simulate hier negbin data and fit hierarchical models
f_sim <- function(trial = 1) {
  #calculate fixed effects using parameter matrix and model matrix
  sum_fix_eff <- fix_mm %*% fix_pars
  
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
           site_obs_t = tweedie::rtweedie(N, mu = mu_vec, power=2, phi=1),
           site_obs = rnbinom(N, mu = mu_vec, size = 2))
}

dat <- f_sim()
hist(dat$site_obs_t)
hist(dat$site_obs)
ggplot(dat) +
  geom_boxplot(aes(x = reg, y = site_obs)) +
  facet_wrap(~seas)
ggplot(dat) +
  geom_point(aes(x = eff_z, y = site_obs))


m2 <- glmmTMB::glmmTMB(site_obs ~ reg + seas + eff_z, data = dat, 
                       family = glmmTMB::nbinom2(link = "log")
                       )
datb <- dat %>% 
  mutate(eff_z2 = eff_z^2)
m2b <- glmmTMB::glmmTMB(site_obs ~ reg + seas + eff_z + eff_z2, data = datb, 
                        family = glmmTMB::nbinom2(link = "log")
)
summary(m2)

m1 <- lm(log(site_obs + 0.0001) ~ reg + seas + eff_z, data = dat)
summary(m1)

#helper function to convert factors 
fct_to_tmb_num <- function(x) {
  as.numeric(as.factor(as.character(x))) - 1
}

fac1k <- fct_to_tmb_num(dat$site)

# order matrix based on unique factor levels with most saturated at bottom
fac_key <- dat %>% 
  select(reg, seas) %>% 
  distinct() %>% 
  mutate(facs = paste(reg, seas, sep = "_"),
         facs_n = as.numeric(as.factor(facs))) %>% 
  arrange(facs_n)
mm_pred1 <- model.matrix(~ reg + seas, data = fac_key)
mm_pred <- cbind(mm_pred1, eff_z = rep(0, n = nrow(mm_pred1)))

# mm_pred1 <- fix_mm[, -5] %>%
#   unique() 
# ord_mat <- mm_pred1 %>% 
#   apply(., 1, sum) %>%
#   sort() %>% 
#   names()
# mm_pred2 <- mm_pred1[ord_mat, ]
# mm_pred <- cbind(mm_pred2, eff_z = rep(0, n = nrow(mm_pred2)))

data <- list(y1_i = dat$site_obs,
             X1_ij = fix_mm,
             factor1k_i = fac1k,
             nk1 = length(unique(fac1k)),
             X1_pred_ij = mm_pred
             )
parameters = list(
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
  log_phi = log(1.1),
  # logit_p = boot::logit(0.8),
  z1_k = rep(0, length(unique(fac1k))),
  log_sigma_zk1 = log(0.25)
  )


## Make a function object
compile(here::here("src/negbin_1re.cpp"))
dyn.load(dynlib(here::here("src/negbin_1re")))
obj <- MakeADFun(data, parameters, random = c("z1_k"), 
                 DLL = "negbin_1re")

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr


## Make predictions
log_pred_fe <- ssdr[rownames(ssdr) %in% "log_prediction", ]
pred_ci <- data.frame(log_pred_est = log_pred_fe[ , "Estimate"],
                      log_pred_se =  log_pred_fe[ , "Std. Error"]) %>%
  mutate(facs_n = fac_key$facs_n) %>% 
  mutate(log_pred_low = log_pred_est + (qnorm(0.025) * log_pred_se),
         log_pred_up = log_pred_est + (qnorm(0.975) * log_pred_se),
         pred_est = exp(log_pred_est),
         pred_se = exp(log_pred_se),
         pred_low = exp(log_pred_low),
         pred_up = exp(log_pred_up)) %>%
  left_join(., fac_key, by = "facs_n") 

ggplot() +
  geom_boxplot(data = dat, aes(x = reg, y = site_obs), alpha = 0.4) +
  geom_pointrange(data = pred_ci, aes(x = as.factor(reg), y = pred_est,
                                      ymin = pred_low, 
                                      ymax = pred_up)) +
  facet_wrap(~seas)
ggplot() +
  geom_boxplot(data = dat, aes(x = reg, y = log(site_obs)), alpha = 0.4) +
  geom_pointrange(data = pred_ci, aes(x = as.factor(reg), y = log_pred_est,
                                      ymin = log_pred_low, 
                                      ymax = log_pred_up)) +
  facet_wrap(~seas)
