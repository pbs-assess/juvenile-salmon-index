# Simulate abundance (tweedie) and multinomial data then fit combined model

library(tidyverse)
library(TMB)

## Common parameters
n_sites <- 5


## Parameters for tweedie data
n_m1_site <- 40
sd_m1_site <- 0.1
N_m1 <- n_sites * n_m1_site
# Fixed 
fix_ints_m1 <- c(0.75, 1) 
f_dat_m1 <- data.frame(reg = as.factor(sample(c(1, 2), size = N_m1, replace = T)))
fix_mm1 <- model.matrix(~ reg, data = f_dat_m1)


## Parameters for multinomial data
n_m2_site <- 250
sd_m2_site <- 0.5
N_m2 <- n_sites * n_m2_site
# Fixed 
group_ints_m2 <- c(0.3, -1.4)
k <- length(group_ints_m2) #number of groups - 1
f1a_ints_m2 <- c(1.5, -0.25) #i.e. how region 2 differs from reference for each group
# reg3_ints <- c(0.27, 0.49) #i.e. how region 3 differs from reference for each group
# fac2_ints <- c(0, -2) #i.e. how factor 2 differs from reference for each group
fix_ints_m2 <- matrix(c(group_ints_m2, f1a_ints_m2), ncol = k, byrow = T)

f_dat_m2 <- data.frame(reg = as.factor(sample(c(1, 2), size = N_m2, 
                                              replace = T)))
# ,
#                       fac = as.factor(sample(c(1, 2), size = N, replace = T)))
fix_mm2 <- model.matrix(~ reg, data = f_dat_m2)


## Simulate data ---------------------------------------------------------------
# function to simulate multinomial data and fit hierarchical models
f_sim <- function(trial = 1) {
  ## Tweedie abundance data
  fe_m1 <- fix_mm1 %*% fix_ints_m1
  # generate random effects, then combine with fixed
  site_mean_m1 <- rnorm(mean = 0, sd = sd_m1_site, n = n_sites)
  mu_fe_m1 <- as.numeric(exp(rep(site_mean_m1, each = n_m1_site) + 
                             fe_m1))
  dat_m1 <- data.frame(trial = rep(trial, length.out = N_m1),
                       site = as.factor(rep(seq(1, n_sites, by = 1), 
                                            each = n_m1_site)),
                       site_mean = rep(site_mean_m1, each = n_m1_site)) %>%
    cbind(f_dat_m1) %>%  #bind regional fixed effects data 
    mutate(fe = as.vector(fe_m1),
           mu = mu_fe_m1,
           site_obs = tweedie::rtweedie(N_m1, mu = mu_fe_m1, 
                                        power = 2, phi = 1))
  
  
  ## Multinomial composition data
  #Calculate multinomial FE
  fe_m2_kk <- matrix(NA, nrow = N_m2, ncol = k)
  for (kk in 1:k) {
    fe_m2_kk[ , kk] <- fix_mm2 %*% fix_ints_m2[, kk] 
  }
  
  # generate random effects, then combine with fixed
  site_mean_m2 <- rnorm(mean = 0,
                       sd = sd_m2_site,
                       n = n_sites)
  site_dat_m2 <- data.frame(site = as.factor(rep(seq(1, n_sites, by = 1), 
                                              each = n_m2_site)),
                         site_m2 = rep(site_mean_m2, each = n_m2_site)) %>% 
    cbind(f_dat_m2) %>%  #bind regional fixed effects data 
    mutate(k1_fe = fe_m2_kk[, 1],
           k2_fe = fe_m2_kk[, 2])
  
  # generate log_odds based on fixed and random intercepts
  dat_m2 <- site_dat_m2 %>% 
    mutate(trial = trial,
           site = as.factor(site),
           log_odds_1_3 = k1_fe + site_m2,
           log_odds_2_3 = k2_fe + site_m2,
           # probability of category 1, 2, and 3:
           p1 = exp(log_odds_1_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
           p2 = exp(log_odds_2_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
           p3 = 1 - (p1 + p2)
    )
  
  # create the observations:
  y <- numeric(length = nrow(dat_m2))
  p <- dat_m2 %>% 
    select(p1:p3) %>%
    as.matrix()
  for (i in seq_along(y)) {
    temp <- rmultinom(1, 1, p[i, ])
    y[i] <- which(temp == 1)
  }
  dat_m2 <- cbind(dat_m2, y)
  
  # matrix version for dmultinom():
  Y_m2 <- matrix(ncol = 3, nrow = N_m2, data = 0)
  for (i in seq_along(y)) Y_m2[i, y[i]] <- 1
  
  return(list(dat_m1 = dat_m1, dat_m2 = dat_m2, obs_m2 = Y_m2))
}

sim_dat <- f_sim()

## Check simulated data
# Abundance data
ggplot(sim_dat$dat_m1) +
  geom_histogram(aes(x = site_obs))
ggplot(sim_dat$dat_m1) +
  geom_boxplot(aes(x = reg, y = site_obs))

# Multinomial data
ggplot(sim_dat$dat_m2) +
    geom_histogram(aes(x = as.factor(y)), stat = "count") + 
    facet_wrap(~reg, nrow = 2)
