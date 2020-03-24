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
m1_dat <- sim_dat$dat_m1
m2_dat <- sim_dat$dat_m2
m2_obs <- sim_dat$obs_m2

## Check simulated data
# Abundance data
ggplot(m1_dat) +
  geom_histogram(aes(x = site_obs))
ggplot(m1_dat) +
  geom_boxplot(aes(x = reg, y = site_obs))

# Multinomial data
ggplot(m2_dat) +
    geom_histogram(aes(x = as.factor(y)), stat = "count") + 
    facet_wrap(~reg, nrow = 2)



## Prep data to pass to TMB
# order matrix based on unique factor levels with most saturated at bottom
mm_pred <- fix_mm1 %>% 
  unique() 
ord_mat <- mm_pred %>% 
  apply(., 1, sum) %>%
  sort() %>% 
  names()
mm_pred <- mm_pred[ord_mat, ]

# fit dummy model to speed up tweedie abundance estimates
m1 <- lm(log(site_obs) ~ reg, data = m1_dat)

#helper function to convert factors 
fct_to_tmb_num <- function(x) {
  as.numeric(as.factor(as.character(x))) - 1
}
fac1k <- fct_to_tmb_num(m1_dat$site)
fac2k <- fct_to_tmb_num(m2_dat$site)

# generate factor key for multinomial data
fac_dat <- m2_dat %>% 
  mutate(facs = reg,#as.factor(paste(reg, sep = "_")),
         facs_n = fct_to_tmb_num(facs)) %>% #subtract for indexing by 0 
  select(reg, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)

# combine
data <- list(
  #abundance data
  y1_i = m1_dat$site_obs,
  X1_ij = fix_mm1,
  factor1k_i = fac1k,
  nk1 = length(unique(fac1k)),
  X1_pred_ij = mm_pred,
  #composition data
  y2_ig = m2_obs,
  X2_ij = fix_mm2,
  factor2k_i = fac2k,
  nk2 = length(unique(fac2k)),
  m2_all_fac = fac_dat$facs_n,
  m2_fac_key = fac_key$facs_n
)

parameters = list(
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),#rep(1, ncol(fix_mm)),
  log_phi = log(1.1),
  logit_p = boot::logit(0.8),
  z1_k = rep(0, length(unique(fac1k))),
  log_sigma_zk1 = log(0.25),
  b2_jg = matrix(0, nrow = ncol(fix_mm2), ncol = ncol(m2_obs) - 1),
  z2_k = rep(0, times = length(unique(fac2k))),
  log_sigma_zk2 = log(0.25)
)

## Make a function object
compile("R/sim_practice/tweediePractice/tweedie_multinomial_1re.cpp")
dyn.load(dynlib("R/sim_practice/tweediePractice/tweedie_multinomial_1re"))
obj <- MakeADFun(data, parameters, random = c("z1_k", "z2_k"), 
                 DLL = "tweedie_multinomial_1re")

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr
