# Simulate and fit unordered multinomial response models with addition of 
# random intercepts

library(tidyverse)

set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
N <- n_sites * n_obs_per_site
sd_site <- 0.5
sd_global <- 1

# function to simulate multinomial data and fit hierarchical models
f_sim <- function(b0 = 0.3, b2 = -1.4) {
  site_mean_a <- rnorm(mean = 0, sd = sd_global, n = n_sites)
  site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                       sd = sd_site, 
                       n = N)
  dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1), 
                                         each = n_obs_per_site)),
                    g1_int = b0, # intercept describing log odds of category 1 vs. 3
                    g2_int = b2, # intercept describing log odds of category 2 vs. 3
                    site_mean = rep(site_mean_a, each = n_obs_per_site),
                    site_obs = site_obs_ia,
                    sd_site = sd_site) 
  # ggplot(dat, aes(x = site, y = site_obs)) +
  #   geom_jitter(width = 0.1, height = 0, alpha = 0.25) +
  #   geom_hline(aes(yintercept = site_mean), color = "red") +
  #   ggsidekick::theme_sleek() +
  #   facet_wrap(~site, scales = "free_x")
  
  # generate log_odds based on fixed and random intercepts
  datf <- dat %>% 
    mutate(site = as.factor(site),
           log_odds_1_3 = g1_int + site_obs,
           log_odds_2_3 = g2_int + site_obs,
           # probability of category 1, 2, and 3:
           p1 = exp(log_odds_1_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
           p2 = exp(log_odds_2_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
           p3 = 1 - (p1 + p2)
    )
  
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
  # ggplot(datf, aes(x = y)) +
  #   geom_histogram() + #negative values in site 2 inc. probability of third cat
  #   ggsidekick::theme_sleek() +
  #   facet_wrap(~site)
  
  # matrix version for dmultinom():
  Y <- matrix(ncol = 3, nrow = N, data = 0)
  for (i in seq_along(y)) Y[i, y[i]] <- 1
  
  return(list("obs" = Y, "rand_fac" = datf$site, "full_data" = datf))
}

dat_list <- vector(mode = "list", length = 100)
for (i in 1:20) {
  dat_list[[i]] <- f_sim()
}


library(TMB)
compile("R/multinomialPractice/multinomial_generic_randInt.cpp")
dyn.load(dynlib("R/multinomialPractice/multinomial_generic_randInt"))

## Data and parameters
y_obs <- dat_list[[1]]$obs
fac1 <- as.numeric(dat_list[[1]]$rand_fac)
data <- list(y_obs = y_obs,
             fac1 = fac1)
parameters <- list(beta1 = rep(0, times = ncol(y_obs) - 1),
                   z_fac1 = rep(0, times = length(unique(fac1))),
                   log_sigma = 0)#,
                   # log_sigma_fac1 = 0)

## Make a function object
obj <- MakeADFun(data, parameters, #random = c("z_fac1"),
                 DLL = "multinomial_generic_randInt")
 

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

ssdr <- summary(sdr)
ssdr

r <- obj$report()
r$probs
r$log_odds
r$logit_probs





trans_sig <- dat_out1 %>%
  filter(var %in% c("sigma", "sigma_site")) %>%
  group_by(var) %>% 
  mutate(est = exp(est)) %>% 
  ungroup()

dat_out1 %>% 
  filter(var %in% c("int1", "int2")) %>%
  rbind(., trans_sig) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = var, y = est)) +
  geom_point(aes(x = var, y = true), colour = "red")


## Don't worry about following code ##

rand_eff <- dat_out1 %>% 
  filter(!var %in% c("int1", "int2", "sigma")) %>% 
  pivot_longer(-var, names_to = "type", values_to = "value")

ggplot(rand_eff) +
  geom_boxplot(aes(x = type, y = value)) +
  facet_wrap(~var)
# looks pretty good


# Predictions ------------------------------------------------------------------
# calculate our predictions:
fitted_log_odds_1_3 <- m2$par[1] + m2$par[2] * X$x
fitted_log_odds_2_3 <- m2$par[3] + m2$par[4] * X$x

fitted_p1 <- exp(fitted_log_odds_1_3) /
  (1 + exp(fitted_log_odds_1_3) + exp(fitted_log_odds_2_3))
fitted_p2 <- exp(fitted_log_odds_2_3) /
  (1 + exp(fitted_log_odds_1_3) + exp(fitted_log_odds_2_3))
fitted_p3 <- 1 - (fitted_p1 + fitted_p2)

fitted <- cbind(fitted_p1, fitted_p2, fitted_p3)
image(t(fitted))


# Second generic attempt -------------------------------------------------------

k <- 3 #number of groups

# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
N <- n_sites * n_obs_per_site
sd_site <- 0.5
sd_global <- 1

# fixed effects for groups (i.e. years)
yr_means <- c(0, 0)
n_yrs <- length(yr_means)

site_mean_a <- rnorm(mean = 0, sd = sd_global, n = n_sites)
site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                     sd = sd_site, 
                     n = N)
yr_obs <- rnorm(mean = rep(yr_means, each = N / length(unique(yr_means))),
                sd = 0.5, 
                n = N)

dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1), 
                                       each = n_obs_per_site)),
                  g1_int = 0.3, # intercept describing log odds of category 1 vs. 3
                  g2_int = -1.4, # intercept describing log odds of category 2 vs. 3
                  site_mean = rep(site_mean_a, each = n_obs_per_site),
                  site_obs = site_obs_ia,
                  #fixed effect
                  yr = sample(c("yr1", "yr2"), size = N, replace = T)) %>% 
  mutate(
    yr_mean = case_when(
      yr == "yr1" ~ yr_means[1],
      yr == "yr2" ~ yr_means[2]
    ),
    yr_obs = 0,#rnorm(mean = yr_mean, sd = 0.5, n = N),
    log_odds_1_3 = g1_int + site_obs + yr_obs,
    log_odds_2_3 = g2_int + site_obs + yr_obs,
    # probability of category 1, 2, and 3:
    p1 = exp(log_odds_1_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
    p2 = exp(log_odds_2_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
    p3 = 1 - (p1 + p2)
  )

# create the observations:
y <- numeric(length = nrow(dat))
p <- dat %>% 
  select(p1:p3) %>%
  as.matrix()

# function to calc
f2 <- function() {
  for (i in seq_along(y)) {
    temp <- rmultinom(1, 1, p[i, ])
    y[i] <- which(temp == 1)
  }
  datf <- cbind(dat, y)
  
  # matrix version for dmultinom():
  Y <- matrix(ncol = 3, nrow = N, data = 0)
  for (i in seq_along(y)) Y[i, y[i]] <- 1
  
  # X <- model.matrix(~ site + yr - 1, dat)
  fix_X <- model.matrix(~ yr - 1, dat)
  rand_X <- model.matrix(~ site - 1, dat)
  
  nll2 <- function(par) {
    #calc fixed effs
    fix_pars = par[3:4]
    fix_e = matrix(NA, nrow = N, ncol = ncol(fix_X))
    for(j in seq_along(fix_pars)) {
      fix_e[ , j] = fix_pars[j] * fix_X[ , j]
    }
    sum_fix_e <- apply(fix_e, 1, sum)
    
    #calc rand effs
    log_sigma <- par[10] #log of random variance
    rand_pars = par[5:9]
    rand_e = matrix(NA, nrow = N, ncol = ncol(rand_X))
    for(j in seq_along(rand_pars)) {
      rand_e[ , j] = rand_pars[j] * rand_X[ , j]
    }
    sum_rand_e <- apply(rand_e, 1, sum) * log_sigma
    
    .log_odds_1_3 <- par[1] + sum_fix_e + sum_rand_e
    .log_odds_2_3 <- par[2] + sum_fix_e + sum_rand_e
    
    .p1 <- exp(.log_odds_1_3) /
      (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))
    .p2 <- exp(.log_odds_2_3) /
      (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))
    .p3 <- 1 - (.p1 + .p2)
    
    nll <- vector(length = length(y))
    for (i in seq_along(y)) { # not vectorized
      nll[i] <- -dmultinom(Y[i,], size = 1,
                           prob = c(.p1[i], .p2[i], .p3[i]), log = TRUE)
    }
    sum(nll)
  }
  
  par_in <- c(rep(0, k-1), #number of groups
              rep(0, n_yrs), #number of fixed effects
              rep(0, n_sites), #number of random effects
              -1) #variance term
  m2 <- nlminb(par_in, nll2)
 
  dat_out <- data.frame(var = c("int1", "int2", levels(datf$yr),
                                levels(datf$site), "sigma"),
                        est = m2$par,
                        true = c(unique(datf$g1_int), unique(datf$g2_int), 
                                 yr_means, site_mean_a, log(sd_global)))
  return(dat_out) 
}
f2()

dat_list <- vector(mode = "list", length = 25)
for (i in 1:25) {
  dat_list[[i]] <- f2()
}
dat_out <- dat_list %>% 
  bind_rows() 

#biased estimates
ggplot(dat_out) +
  geom_boxplot(aes(x = var, y = est)) +
  geom_point(aes(x = var, y = true), colour = "red")

rand_eff <- dat_out %>% 
  filter(!var %in% c("int1", "int2", "sigma")) %>% 
  pivot_longer(-var, names_to = "type", values_to = "value")

ggplot(rand_eff) +
  geom_boxplot(aes(x = type, y = value)) +
  facet_wrap(~var)










## Code chunk for first generic version
##negative log likelihood function using dmultinom:
nll2 <- function(par) {
  log_sigma <- par[8]
  log_sigma_r <- par[9]

  # multiply random coefficients by model matrix to get linear predictor
  z1_k <- par[3:7]
  z1_i <- XX %*% z1_k
  .log_odds_1_3b <- par[1] + (z1_i * exp(log_sigma))
  .log_odds_2_3b <- par[2] + (z1_i * exp(log_sigma))

  .p1 <- exp(.log_odds_1_3) /
    (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))
  .p2 <- exp(.log_odds_2_3) /
    (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))
  .p3 <- 1 - (.p1 + .p2)

  nll <- vector(length = length(y))
  for (i in seq_along(y)) { # not vectorized
    nll[i] <- -dmultinom(Y[i,], size = 1,
                         prob = c(.p1[i], .p2[i], .p3[i]), log = TRUE)
  }
  # probability of random coefficients
  # sum(nll)
  nll_r <- vector(length = length(z1_k))
  for (k in seq_along(z1_k)) {
    nll_r[k] <- -dnorm(z1_k[k], 0, exp(log_sigma_r), log = TRUE)
  }
  sum(nll, nll_r)
}

par_in <- c(0, 0, rep(0, n_sites), -0.75, -0.75)
m2 <- nlminb(par_in, nll2)

dat_out <- data.frame(var = c("int1", "int2", unique(datf$site), "sigma",
                              "sigma_site"),
                      est = m2$par,
                      true = c(b0, b2, site_mean_a, sd_global, sd_site))
return(dat_out)