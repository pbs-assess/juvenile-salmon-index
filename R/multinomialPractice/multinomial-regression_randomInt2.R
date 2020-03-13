# Simulate and fit unordered multinomial response models with addition of 
# random intercepts; extension of multinomial-regression_randomInt that also
# includes two additional fixed effects

library(tidyverse)

set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 250
sd_site <- 0.5
N <- n_sites * n_obs_per_site

# fixed intercepts
group_ints <- c(0.3, -1.4)
k <- length(group_ints) #number of groups - 1
reg2_ints <- c(1.5, -0.25) #i.e. how region 2 differs from reference for each group
fac2_ints <- c(0, -0.75) #i.e. how region 2 differs from reference for each group
fix_ints <- matrix(c(reg2_ints, fac2_ints
                 ), ncol = k, byrow = T)

reg_dat <- data.frame(reg = as.factor(sample(c(1, 2), size = N, replace = T)),
                      fac = as.factor(sample(c(1, 2), size = N, replace = T)))

# model matrix for fixed effects
fix_mm <- model.matrix(~ reg:fac - 1, reg_dat)

# function to simulate multinomial data and fit hierarchical models
f_sim <- function(trial = 1) {
  
  #calculate fixed effects using parameter matrix and model matrix
  n_fix_cov <- nrow(ints)
  sum_fix_eff <- matrix(NA, nrow = N, ncol = k)
  for (kk in 1:k) {
    fix_eff <- matrix(nrow = N, ncol = n_fix_cov)
    for (i in 1:N) {
      for (j in 1:n_fix_cov) {
        fix_eff[i, j] <- group_ints[kk] + (ints[j, kk] * fix_mm[i, j])  
      }
    }
    sum_fix_eff[ , kk] <- apply(fix_eff, 1, sum)
  }

  # generate random effects, then combine with fixed
  site_mean_a <- runif(n = n_sites, min = -1, max = 1)
  site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                       sd = sd_site, 
                       n = N)
  site_dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1), 
                                         each = n_obs_per_site)),
                         site_mean = rep(site_mean_a, each = n_obs_per_site),
                         site_obs = site_obs_ia,
                         sd_site = sd_site
                    ) %>% 
    cbind(reg_dat) %>%  #bind regional fixed effects data 
    # arrange(site, reg) %>% 
    mutate(g1_fe = sum_fix_eff[, 1],
           g2_fe = sum_fix_eff[, 2])
  
  # generate log_odds based on fixed and random intercepts
  datf <- site_dat %>% 
    mutate(site = as.factor(site),
           log_odds_1_3 = g1_fe + site_obs,
           log_odds_2_3 = g2_fe + site_obs,
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
  #   # facet_wrap(~reg, nrow = 2)
  # facet_wrap(reg~site, nrow = 2)
  
  # matrix version for dmultinom():
  Y <- matrix(ncol = 3, nrow = N, data = 0)
  for (i in seq_along(y)) Y[i, y[i]] <- 1
  
  # true values vector
  tv <- c(as.vector(ints), log(sd_site), site_mean_a, sd_site)
  
  return(list("trial" = trial, "obs" = Y, "fixed_fac" = site_dat$reg, 
              "rand_fac" = datf$site, "trim_data" = tv, 
              "full_data" = datf))
}

dat_list <- vector(mode = "list", length = 100)
for (i in seq_along(dat_list)) {
  dat_list[[i]] <- f_sim(trial = i)
}


library(TMB)
compile("R/multinomialPractice/multinomial_generic_randInt_fixInt.cpp")
dyn.load(dynlib("R/multinomialPractice/multinomial_generic_randInt_fixInt"))

## Data and parameters
y_obs <- dat_list[[1]]$obs
rfac <- as.numeric(dat_list[[1]]$rand_fac) - 1 #subtract for indexing by 0
data <- list(y_obs = y_obs,
             rfac = rfac,
             fx_cov = fix_mm, #model matrix from above
             n_rfac = length(unique(rfac))
             )
parameters <- list(z_rfac = rep(0, times = length(unique(rfac))),
                   z_ints = matrix(0, nrow = ncol(fix_mm),
                                   ncol = ncol(y_obs) - 1),
                   log_sigma_rfac = 0)

## Make a function object
obj <- MakeADFun(data, parameters, random = c("z_rfac"),
                 DLL = "multinomial_generic_randInt_fixInt")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

ssdr <- summary(sdr)
ssdr


ests <- map(dat_list, function(x) {
  y_obs <- x$obs
  rfac <- as.numeric(x$rand_fac) - 1 #subtract for indexing by 0
  data <- list(y_obs = y_obs,
               rfac = rfac,
               fx_cov = fix_mm, #model matrix from above
               n_rfac = length(unique(rfac))
  )
  parameters <- list(z_rfac = rep(0, times = length(unique(rfac))),
                     z_ints = matrix(0, nrow = ncol(fix_mm),
                                     ncol = k),
                     log_sigma_rfac = 0)
  ## Make a function object
  obj <- MakeADFun(data, parameters, random = c("z_rfac"),
                   DLL = "multinomial_generic_randInt_fixInt")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  fac_seq <- paste("z_rfac", seq(1, length(unique(rfac)), by = 1),
                   sep = "_")
  as.data.frame(ssdr) %>% 
    mutate(par = c("g1_int", "g2_int", "reg2_g1", "reg2_g2", "log_sigma_rfac", 
                   fac_seq,
                   "sigma_fac1"),
           trial = x$trial,
           vals = x$trim_data)
}) %>% 
  bind_rows()


ests %>% 
  filter(par %in% c("g1_int", "g2_int", "reg2_g1", "reg2_g2", "sigma_fac1")) %>%
  ggplot(.) +
  geom_boxplot(aes(x = par, y = Estimate)) +
  geom_point(aes(x = par, y = vals), colour = "red")

r_effs <- ests %>% 
  filter(grepl("z_rfac", par)) %>%
  mutate(p_err = (Estimate - vals) / vals) 

r_effs %>% 
  select(par, est = Estimate, true = vals) %>% 
  pivot_longer(-par, names_to = "type", values_to = "value") %>% 
  ggplot(.) +
  geom_boxplot(aes(x = type, y = value)) +
  facet_wrap(~par)

ggplot(r_effs) +
  geom_boxplot(aes(x = par, y = p_err))

