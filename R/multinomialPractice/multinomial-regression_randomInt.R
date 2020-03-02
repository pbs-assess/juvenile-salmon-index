# Simulate and fit unordered multinomial response models

# e.g. https://www.bristol.ac.uk/media-library/sites/cmm/migrated/documents/unordered-multi-r-models.pdf
# but note a number of typos in the equations in those slides!

library(tidyverse)
set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
N <- n_sites * n_obs_per_site
sd_site <- 1
sd_global <- 5

site_mean_a <- rnorm(mean = 0, sd = sd_global, n = n_sites)
site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                     sd = sd_site, 
                     n = N)
dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1), 
                                       each = n_obs_per_site)),
                  b0 = 0.3, # intercept describing log odds of category 1 vs. 3
                  b2 = -1.4, # intercept describing log odds of category 2 vs. 3
                  site_mean = rep(site_mean_a, each = n_obs_per_site),
                  site_obs = site_obs_ia,
                  sd_site = sd_site) 
ggplot(dat, aes(x = site, y = site_obs)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.25) +
  geom_hline(aes(yintercept = site_mean), color = "red") +
  ggsidekick::theme_sleek() +
  facet_wrap(~site, scales = "free_x")

# add fixed effects (log_odds relative to reference category)
datf <- dat %>% 
  mutate(site = as.factor(site),
         log_odds_1_3 = b0 + site_obs,
         log_odds_2_3 = b2 + site_obs,
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
ggplot(datf, aes(x = y)) +
  geom_histogram() + #negative values in site 2 inc. probability of third cat
  ggsidekick::theme_sleek() +
  facet_wrap(~site)

# create temporary vectors for fitting:
y_1 <- as.numeric(y == 1)
y_2 <- as.numeric(y == 2)

# matrix version for dmultinom():
Y <- matrix(ncol = 3, nrow = N, data = 0)
for (i in seq_along(y)) Y[i, y[i]] <- 1

# model matrix treating RE as FE
XX <- model.matrix(~ site + b0 - 1, datf) %>% 
  as.data.frame()


# negative log likelihood function using dbinom:
# nll <- function(par) {
#   .log_odds_1_3 <- par[1] + par[2] * X$x
#   .log_odds_2_3 <- par[3] + par[4] * X$x
# 
#   .p1 <- exp(.log_odds_1_3) /
#     (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))
#   .p2 <- exp(.log_odds_2_3) /
#     (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))
# 
#   -sum(dbinom(y_1, size = 1, prob = .p1, log = TRUE) +
#     dbinom(y_2, size = 1, prob = .p2, log = TRUE))
# }
# 
# m <- nlminb(rep(0, 4), nll)
# m

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
# This should be the same but is slightly different(!?).
par_in <- c(0 , 0, rep(0, n_sites), 0.5)
m2 <- nlminb(par_in, nll2)
m2

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

