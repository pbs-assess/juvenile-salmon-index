# TMB ------------------------------------------------------

set.seed(1)
N <- 2000 # number of observations
k <- 4 #number of groups
# X <- runif(N) # cont. predictor
#fixed predictor representing two cats
yrCov <- data.frame(yr = sample(c("yr1", "yr2", "yr3"), size = N, replace = T)) 
X <- model.matrix(~yr, yrCov)

int <- c(0.3, -1.4, 0.3)
beta_yr2 <- c(-0.25, 0.75, -0.5)
beta_yr3 <- c(0.5, -1, 0.1)

# append with reference category
#combine intercepts and betas into one matrix that can be multiplied by model
#matrix; columns are parameters for each group (excluding ref. cat.)
betas <- matrix(c(int, beta_yr2, beta_yr3), nrow = ncol(X), ncol = k - 1, 
                byrow = T)
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
# plot(X[ , 2], jitter(y, 0.3))

#matrix of observations
y_obs <- matrix(ncol = k, nrow = N, data = 0)
for (i in seq_along(y)) {
  y_obs[i, y[i]] <- 1
}

library(TMB)
compile("R/sim_practice/multinomialPractice/multinomial_generic.cpp")
dyn.load(dynlib("R/sim_practice/multinomialPractice/multinomial_generic"))

## Data and parameters
# .X <- cbind(1, X) #predictor with intercept
.X <-  model.matrix(~yr, yrCov) #intercept already bound
data <- list(cov = .X, y_obs = y_obs)
parameters <- list(betas = matrix(data = 0, nrow = ncol(.X), 
                                  ncol = ncol(y_obs) - 1))

## Make a function object
obj <- MakeADFun(data, parameters, DLL="multinomial_generic")

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


# compare with vanilla multinomial model
long_dat <- y_obs %>% 
  cbind(yrCov, .) %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = `1`:`4`) %>% 
  filter(!value == 0) %>% 
  select(-value) %>% 
  mutate(name = as.factor(name))
m1 <- nnet::multinom(name ~ yr, data = long_dat)
summary(m1)
exp(coef(m1))
head(pp <- fitted(m1))


# replace all category 4 covariates to yr1
cat_4_cov <- which(y_obs[, 4] == 1)
yrCov[cat_4_cov, ] <- "yr2"
yrCov[cat_4_cov, ]

#rerun basic model
long_dat <- y_obs %>% 
  cbind(yrCov, .) %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = `1`:`4`) %>% 
  filter(!value == 0) %>% 
  select(-value) %>% 
  mutate(name = as.factor(name))
m2 <- nnet::multinom(name ~ yr, data = long_dat)
summary(m2)
exp(coef(m2))

#rerun tmb model
obj <- MakeADFun(data, parameters, DLL="multinomial_generic")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
# SEs blow up on parameter estimates, but predictions are surprisingly stable

head(pp <- fitted(m2))
r <- obj$report()
head(r$probs)
