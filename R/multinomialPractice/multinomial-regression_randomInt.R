# Simulate and fit unordered multinomial response models with addition of 
# random intercepts

# e.g. https://www.bristol.ac.uk/media-library/sites/cmm/migrated/documents/unordered-multi-r-models.pdf
# but note a number of typos in the equations in those slides!

library(tidyverse)

set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
N <- n_sites * n_obs_per_site
sd_site <- 0.5
sd_global <- 1

# function to simulate multinomial data and fit hierarchical models
f <- function(b0 = 0.3, b2 = -1.4) {
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
                    sd_global = sd_global) 
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
  
  # model matrix treating RE as FE
  XX <- model.matrix(~ site - 1, datf) %>% 
    as.data.frame()
  
  # negative log likelihood function using dmultinom:
  nll2 <- function(par) {
    log_sigma <- log(par[8])
    .log_odds_1_3 <- par[1] + (par[3] * XX$site1 + par[4] * XX$site2 + 
                                 par[5] * XX$site3 + par[6] * XX$site4 + par[7] * XX$site5) * log_sigma
    .log_odds_2_3 <- par[2] + (par[3] * XX$site1 + par[4] * XX$site2 + 
                                 par[5] * XX$site3 + par[6] * XX$site4 + par[7] * XX$site5) * log_sigma
    
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
  
  par_in <- c(0 , 0, rep(0, n_sites), 1.5)
  m2 <- nlminb(par_in, nll2)
  
  dat_out <- data.frame(var = c("int1", "int2", unique(datf$site), "sigma"),
                        est = m2$par,
                        true = c(b0, b2, site_mean_a, sd_site))
  return(dat_out)
}

dat_list <- vector(mode = "list", length = 100)
for (i in 1:100) {
  dat_list[[i]] <- f()
}
dat_out <- dat_list %>% 
  bind_rows() 

trans_sig <- dat_out %>% 
  filter(var == "sigma") %>% 
  mutate(var = fct_recode(var, log_sigma = "sigma"),
         est = log(est))

dat_trim <- dat_out %>% 
  filter(var %in% c("int1", "int2", "sigma")) %>% 
  rbind(., trans_sig)

ggplot(dat_trim) +
  geom_boxplot(aes(x = var, y = est)) +
  geom_point(aes(x = var, y = true), colour = "red")

rand_eff <- dat_out %>% 
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

site_mean_a <- rnorm(mean = 0, sd = sd_global, n = n_sites)
site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                     sd = sd_site, 
                     n = N)
dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1), 
                                       each = n_obs_per_site)),
                  # g1_int = 0.3, # intercept describing log odds of category 1 vs. 3
                  # g2_int = -1.4, # intercept describing log odds of category 2 vs. 3
                  site_mean = rep(site_mean_a, each = n_obs_per_site),
                  site_obs = site_obs_ia,
                  #fixed effect
                  yr = sample(c("yr1", "yr2"), size = N, replace = T)) 

X <- model.matrix(~ yr:site - 1, dat)

int <- c(0.3, -1.4, 0.3)
beta <- c(-0.25, 0.75, -0.5)

# append with reference category
#combine intercepts and betas into one matrix that can be multiplied by model
#matrix; columns are parameters for each group (excluding ref. cat.)
betas <- matrix(c(int, beta), nrow = ncol(X), ncol = k - 1, byrow = T)
m <- nrow(betas) #number of covariates (i.e. 1 int. + slopes)

log_odds <- matrix(NA, nrow = N, ncol = k)
for (h in 1:k) {
  if (h < k) {
    cov_effects <- matrix(NA, nrow = N, ncol = m)
    for (i in seq_len(N)) {
      for (j in seq_len(m)) {
        cov_effects[i, j] <- betas[j , h] * X[i, j]  
      }
      log_odds[i, h] <- sum(cov_effects[i, ])
    }
  } else if (h == k) {
    log_odds[ , h] <- NA
  }
}
exp_log_odds <- exp(log_odds)

denominator <- rep(NA, length.out = N)
for (i in 1:N) {
  denominator[i] <- 1 + sum(exp_log_odds[i, 1:(k-1)])
}  

probs <- matrix(NA, nrow = N, ncol = k)
for (i in 1:k) {
  if (i < k) {
    probs[ , i] <- exp_log_odds[ , i] / denominator
  } 
  else if (i == k) {
    for (h in 1:N) {
      probs[h, i] <- 1 - sum(probs[h, 1:(k-1)]) 
    }
  }
}

#generate observations
y <- numeric(length = N)
for (i in seq_along(y)) {
  temp <- rmultinom(1, 1, probs[i, ])
  y[i] <- which(temp == 1)
}
plot(X[ , 2], jitter(y, 0.3))
