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
    i} else if (i == k) {
      for (h in 1:N) {
        probs[h, i] <- 1 - sum(probs[h, 1:(k-1)]) 
      }
    }
}


#probability of each category
probs <- matrix(NA, nrow = N, ncol = k)
colnames(probs) <- unique(effects$cat)
for (i in 1:k) {
  if (i < k) {
    dum <- effects %>% 
      filter(cat == i)
    probs[ , i] <- dum$exp_log_odds / denominator
  i} else if (i == k) {
    probs[ , i] <- apply(probs[ , 1:(k-1)], 1, function(x) 1 - sum(x))
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
X <- runif(N) #predictor
ints <- c(0.3, -1.4)
betas <- c(-3, 4)
k <- 3


nll2 <- function(ints, betas, k, X) {
  # .log_odds_1_3 <- par[1] + par[2] * X$x
  # .log_odds_2_3 <- par[3] + par[4] * X$x
  cat <- seq(1, k, by = 1) %>% 
    rep(., each = length(X))
  int <- c(ints, NaN) %>% 
    rep(., each = length(X))
  beta <- c(betas, NaN) %>% 
    rep(., each = length(X))
  effects <- data.frame(cat,
                        int,
                        beta,
                        x = rep(X, times = k)) %>% 
    mutate(log_odds = case_when(
      cat != refCat ~ int + beta * x,
      cat == refCat ~ NaN),
      exp_log_odds = exp(log_odds))
  
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
m2 <- nlminb(rep(0, 4), nll2)
m2