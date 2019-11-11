## Simulate and fit unordered multinomial response models with k categories
# Based on multinomial-regression.R

set.seed(42)
N <- 100 # number of observations
X <- data.frame(x = runif(N)) # our predictor
# for the following, category 3 will be our reference category
b0 <- 0.3 # intercept describing log odds of category 1 vs. 3
b1 <- -3 # effect on log odds of category 1 vs. 3 for one unit change in `x`
b2 <- -1.4 # intercept describing log odds of category 2 vs. 3
b3 <- 4 # effect on log odds of category 2 vs. 3 for one unit change in `x`

k <- 3 #number of groups
refCat <- "3"

library(tidyverse)

# structure inputs
set.seed(42)
N <- 100 # number of observations
k <- 3 #number of groups

### Gen function 1 - simulate data from variable number of groups

#Gen function 1 inputs 
X <- runif(N) #predictor
ints <- c(0.3, -1.4)
betas <- c(-3, 4)
k <- 3

# Function guts
cat <- seq(1, k, by = 1)
int <- c(ints, NaN) 
beta <- c(betas, NaN)

log_odds <- matrix(NA, nrow = N, ncol = k)
for (i in 1:k) {
  if (i < k) {
    log_odds[ , i] <- int[i] + beta[i] * X
  } else if (i == k) {
    log_odds[ , i] <- NA
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
      #not sure if we can vectorize this without apply...
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
plot(X, jitter(y, 0.3))
head(cbind(x = X, y))

#matrix of observations
Y <- matrix(ncol = k, nrow = N, data = 0)
for (i in seq_along(y)) {
  Y[i, y[i]] <- 1
}


### Generic function 2 - estimate probabilities from input matrix
#Gen function 2 inputs 
ints <- c(0.3, -1.4)
betas <- c(-3, 4)
k <- 3 #number of groups

# pars <- c(0.3, -1.4, -3, 4, 1, 1, 2, 2)
pars <- c(0.3, -1.4, -3, 4)
X <- matrix(data = runif(N), nrow = N, ncol = 1) #predictor

nll2 <- function(pars) {
  # ideally k and the covariate matrix should be arguments within the function,
  # but not sure how to pass that to nlminb 
  k = 3
  cov = X
  
  # define intercepts and betas based on pars vector and add third reference 
  # category to vector of parameters
  numBetas <- (length(pars) / (k-1)) - 1
  int <- c(pars[1:(k-1)], NA)
  betasVec <- pars[k:length(pars)]
  beta <- matrix(NA, nrow = numBetas, ncol = k) 
  for (j in 1:numBetas) {
    beta[j , ] <- c(betasVec[(j * 2 - 1):(j * 2)], NA)
  }
  N <- nrow(cov)
  
  ### adjusted to allow for multiple covariate effects, but not very elegant...
  # calculate log-odds for each category 
  log_odds <- matrix(NA, nrow = N, ncol = k)
  for (h in 1:k) {
    if (h < k) {
      #calculate cumulative covariate effects
      covEffects <- matrix(NA, nrow = N, ncol = numBetas)
      for (j in 1:numBetas) {
        covEffects[ , j] <- beta[j, h] * cov[ , j]
      }
      for (i in 1:N) {
        log_odds[i , h] <- int[h] + sum(covEffects[i, ])
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
      for (i in 1:N) {
        #not sure if we can vectorize this without apply...
        probs[i, h] <- 1 - sum(probs[i, 1:(k-1)]) 
      }
    }
  }
  
  #generate matrix of observations
  y <- numeric(length = N)
  for (i in 1:N) {
    temp <- rmultinom(1, 1, probs[i, ])
    y[i] <- which(temp == 1)
  }
  Y <- matrix(ncol = k, nrow = N, data = 0)
  for (i in 1:N) {
    Y[i, y[i]] <- 1
  }
  
  #log-likelihood
  nll <- vector(length = length(y))
  for (i in seq_along(y)) { # not vectorized
    nll[i] <- -dmultinom(Y[i, ], size = 1, prob = probs[i, ], log = TRUE)
  }
  sum(nll)
}
# Fit with initial values
m2 <- nlminb(rep(0, 4), nll2)
m2
