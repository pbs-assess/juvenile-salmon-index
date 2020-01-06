## Visualize random fields simulation (included in 
# generic_multinomial-regression-tmb.R) 

library(sdmTMB)
library(ggplot2)
library(dplyr)

.rf_sim <- function(model, x, y) {
  out <- sdmTMB:::rf_sim(model, x, y)
  out - mean(out)
}

N <- 300 # number of observations
k <- 4 #number of groups

log_odds <- matrix(NA, nrow = N, ncol = (k - 1))
sp_list <- vector(mode = "list", length = k)
coords <- expand.grid(x = seq(0, 1, length.out = round(1.25*sqrt(N), 0)), 
                      y = seq(0, 1, length.out = round(1.25*sqrt(N), 0))) %>% 
  sample_n(., size = N)
rf_pars <- list(sig = 0.9, kappa = 1.7)
# rf_pars <- NULL
betas <- c(0, 0, 0) #slopes
ints <- c(0, 0, 0) #slopes
X <- runif(N)

set.seed(123)
for (h in 1:(k - 1)) {
  # normal calculation of log_odds w/out random fields
  if (is.null(rf_pars)) {
    log_odds[ , h] <- ints[h] + betas[h] * X
  }
  
  # with simulated draws from random field 
  if (!is.null(rf_pars)) {
    rf_model <- RandomFields::RMmatern(nu = 1, var = rf_pars$sig^2, 
                                       scale = 1/rf_pars$kappa)
    #assumes relatively coarse grid size currently
    rf_dat <- coords %>% 
      mutate(rf_effect = .rf_sim(model = rf_model, coords$x, coords$y)) %>% 
      sample_n(., size = N)
    
    #### TEMP LIST TO LOOK AT SPATIAL OBSERVATIONS ####
    sp_list[[h]] <- rf_dat      
    
    log_odds[ , h] <- ints[h] + betas[h] * X + rf_dat$rf_effect
  }
}
exp_log_odds <- exp(log_odds)

denominator <- rep(NA, length.out = N)
for (i in 1:N) {
  denominator[i] <- 1 + sum(exp_log_odds[i, ])
}

probs <- matrix(NA, nrow = N, ncol = k)
for (h in 1:k) {
  if (h < k) {
    probs[ , h] <- exp_log_odds[ , h] / denominator
  } 
  else if (h == k) {
    for (i in 1:N) {
      probs[i, h] <- 1 - sum(probs[i, 1:(k-1)]) 
    }
  }
}

#generate observations
y <- numeric(length = N)
for (i in seq_along(y)) {
  temp <- rmultinom(1, 1, probs[i, ])
  y[i] <- which(temp == 1)
}

#matrix of observations
y_obs <- matrix(ncol = k, nrow = N, data = 0)
for (i in seq_along(y)) {
  y_obs[i, y[i]] <- 1
}

## temporary vis
temp <- NULL
for (i in 1:3) {
  dum <- cbind(det = y_obs[ , i], sp_list[[i]][, 1:2])
  dum$k <- i
  dum$rf_effect <- sp_list[[i]]$rf_effect
  temp <- rbind(temp, dum)
}

rf <- ggplot(temp, aes(x, y, fill = rf_effect)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~k, ncol = 1)
obs_det <- ggplot(temp, aes(x, y, fill = det)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~k, ncol = 1)

ggpubr::ggarrange(rf, obs_det, ncol = 2)

# ------------------------------------------------------------------------------

## Random fields sim from S. Anderson
set.seed(122)

d <- expand.grid(X = seq(0, 1, length.out = 100), 
                 Y = seq(0, 1, length.out = 100), year = 1:5)

.rf_sim <- function(model, x, y) {
  out <- sdmTMB:::rf_sim(model, x, y)
  out - mean(out)
}

sigma_O <- 0.3
sigma_Z <- 0.3
sigma_E <- 0.15
kappa <- 1.7
# as kappa shrinks less granularity

rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1/kappa)
d$omega_s <- .rf_sim(model = rf_omega, d$X, d$Y)

ggplot(dplyr::filter(d, year == 1), aes(X, Y, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2()

rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1/kappa)

d <- d %>% group_by(year) %>%
  group_split() %>%
  purrr::map_dfr(
    ~tibble(X = .x$X, Y = .x$Y, omega_s = .x$omega_s, year = .x$year,
            eps_st = .rf_sim(model = rf_epsilon, .x$X, .x$Y))) %>%
  ungroup()

ggplot(d, aes(X, Y, fill = eps_st)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_grid(cols = vars(year))

ggplot(d, aes(X, Y, fill = eps_st + omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_grid(cols = vars(year))

rf_zeta <- RandomFields::RMmatern(nu = 1, var = sigma_Z^2, scale = 1/kappa)
d$zeta_s <- .rf_sim(model = rf_zeta, d$X, d$Y)
d$year_cent <- d$year - mean(d$year)
d$`With spatially\nvarying trend` <- d$omega_s + d$eps_st + d$zeta_s * d$year_cent
d$`Without spatially\nvarying trend` <- (d$omega_s + d$eps_st) * 1.5 # approx. equalize variance

g <- reshape2::melt(d, id.vars = c("X", "Y", "year"),
                    measure.vars = c("With spatially\nvarying trend", 
                                     "Without spatially\nvarying trend")) %>%
  ggplot(aes(X, Y, fill = value)) +
  geom_raster() +
  # scale_fill_gradient2(high = scales::muted("red"), mid = "grey94",
  #   low = scales::muted("blue")) +
  scale_fill_viridis_c(option = "A") +
  facet_grid(cols = vars(year), rows = vars(variable)) +
  ggsidekick::theme_sleek() +
  theme(panel.border = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        panel.spacing.x = unit(-0.15, "lines"), panel.spacing.y = unit(-0.25, "lines"),
        axis.ticks = element_blank()) +
  guides(fill = FALSE)
print(g)


