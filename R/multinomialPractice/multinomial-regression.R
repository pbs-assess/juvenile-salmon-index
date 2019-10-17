# Simulate and fit unordered multinomial response models

# e.g. https://www.bristol.ac.uk/media-library/sites/cmm/migrated/documents/unordered-multi-r-models.pdf
# but note a number of typos in the equations in those slides!

set.seed(42)
N <- 100 # number of observations
X <- data.frame(x = runif(N)) # our predictor
# for the following, category 3 will be our reference category
b0 <- 0.3 # intercept describing log odds of category 1 vs. 3
b1 <- -3 # effect on log odds of category 1 vs. 3 for one unit change in `x`
b2 <- -1.4 # intercept describing log odds of category 2 vs. 3
b3 <- 4 # effect on log odds of category 2 vs. 3 for one unit change in `x`

log_odds_1_3 <- b0 + b1 * X$x
log_odds_2_3 <- b2 + b3 * X$x

# probability of category 1, 2, and 3:
p1 <- exp(log_odds_1_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3))
p2 <- exp(log_odds_2_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3))
p3 <- 1 - (p1 + p2)

# create the observations:
y <- numeric(length = nrow(X))
for (i in seq_along(y)) {
  temp <- rmultinom(1, 1, c(p1[i], p2[i], p3[i]))
  y[i] <- which(temp == 1)
}
plot(X$x, jitter(y, 0.3))
head(cbind(x = X$x, y))

# create temporary vectors for fitting:
y_1 <- as.numeric(y == 1)
y_2 <- as.numeric(y == 2)

# matrix version for dmultinom():
Y <- matrix(ncol = 3, nrow = N, data = 0)
for (i in seq_along(y)) Y[i,y[i]] <- 1

# negative log likelihood function using dbinom:
nll <- function(par) {
  .log_odds_1_3 <- par[1] + par[2] * X$x
  .log_odds_2_3 <- par[3] + par[4] * X$x

  .p1 <- exp(.log_odds_1_3) /
    (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))
  .p2 <- exp(.log_odds_2_3) /
    (1 + exp(.log_odds_1_3) + exp(.log_odds_2_3))

  -sum(dbinom(y_1, size = 1, prob = .p1, log = TRUE) +
    dbinom(y_2, size = 1, prob = .p2, log = TRUE))
}

m <- nlminb(rep(0, 4), nll)
m

# negative log likelihood function using dmultinom:
nll2 <- function(par) {
  .log_odds_1_3 <- par[1] + par[2] * X$x
  .log_odds_2_3 <- par[3] + par[4] * X$x

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

