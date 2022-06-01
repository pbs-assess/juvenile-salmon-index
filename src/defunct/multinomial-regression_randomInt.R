# Simulate and fit unordered multinomial response models with addition of 
# random intercepts

library(tidyverse)

set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
N <- n_sites * n_obs_per_site
sd_site <- 0.5
# sd_global <- 1

# function to simulate multinomial data and fit hierarchical models
f_sim <- function(trial = trial, b0 = 0.3, b2 = -1.4) {
  site_mean_a <- runif(n = n_sites, min = -1, max = 1)
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
  
  # matrix version for dmultinom():
  Y <- matrix(ncol = 3, nrow = N, data = 0)
  for (i in seq_along(y)) Y[i, y[i]] <- 1
  
  # true values vector
  tv <- c(b0, b2, log(sd_site), site_mean_a, sd_site) 
  
  return(list("trial" = trial, "obs" = Y, "rand_fac" = datf$site, 
              "trim_data" = tv, "full_data" = datf))
}

dat_list <- vector(mode = "list", length = 100)
for (i in 1:200) {
  dat_list[[i]] <- f_sim(trial = i)
}


library(TMB)
compile("R/multinomialPractice/multinomial_generic_randInt.cpp")
dyn.load(dynlib("R/multinomialPractice/multinomial_generic_randInt"))

## Data and parameters
y_obs <- dat_list[[1]]$obs
fac1 <- as.numeric(dat_list[[1]]$rand_fac) - 1 #subtract for indexing by 0
data <- list(y_obs = y_obs,
             fac1 = fac1,
             n_fac = length(unique(fac1)))
parameters <- list(beta1 = rep(0, times = ncol(y_obs) - 1),
                   z_fac1 = rep(0, times = length(unique(fac1))),
                   #log_sigma = 0,
                   log_sigma_fac1 = 0)

## Make a function object
obj <- MakeADFun(data, parameters, random = c("z_fac1"),
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

ests <- map(dat_list, function(x) {
  y_obs <- x$obs
  fac1 <- as.numeric(x$rand_fac) - 1 #subtract for indexing by 0
  data <- list(y_obs = y_obs,
               fac1 = fac1,
               n_fac = length(unique(fac1)))
  parameters <- list(beta1 = rep(0, times = ncol(y_obs) - 1),
                     z_fac1 = rep(0, times = length(unique(fac1))),
                     # log_sigma = 0,
                     log_sigma_fac1 = 0)
  ## Make a function object
  obj <- MakeADFun(data, parameters, random = c("z_fac1"),
                   DLL = "multinomial_generic_randInt")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  fac_seq <- paste("z_fac1", seq(1, length(unique(fac1)), by = 1),
                   sep = "_")
  as.data.frame(ssdr) %>% 
    mutate(par = c("beta1", "beta2", "log_sigma_fac1", 
                   fac_seq,
                   "sigma_fac1"),
           trial = x$trial,
           vals = x$trim_data)
}) %>% 
  bind_rows()

ests %>% 
  filter(par %in% c("beta1", "beta2", "sigma_fac1")) %>%
  ggplot(.) +
  geom_boxplot(aes(x = par, y = Estimate)) +
  geom_point(aes(x = par, y = vals), colour = "red")

r_effs <- ests %>% 
  filter(grepl("z_fac1", par)) %>%
  mutate(p_err = (Estimate - vals) / vals) 

r_effs %>% 
  select(par, est = Estimate, true = vals) %>% 
  pivot_longer(-par, names_to = "type", values_to = "value") %>% 
  ggplot(.) +
  geom_boxplot(aes(x = type, y = value)) +
  facet_wrap(~par)

ggplot(r_effs) +
  geom_boxplot(aes(x = par, y = p_err))

