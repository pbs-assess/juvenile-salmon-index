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

# structure inputs
effects <- data.frame(cat = seq(1, k, by = 1),
                      int = c(0.3, -1.4, NaN),
                      beta = c(-3, 4, NaN)) %>% 
  expand(nesting(cat, int, beta), x = runif(N)) %>%
  mutate(log_odds = case_when(
    cat != refCat ~ int + beta * x,
    cat == refCat ~ NaN),
    exp_log_odds = exp(log_odds))

##denominator for calculating probabilities
denominator <- effects %>% 
  group_by(x) %>% 
  summarize(den = 1 + sum(exp_log_odds, na.rm = T)) %>% 
  pull(den) 

probs <- matrix(NA, nrow = N, ncol = k)
for (i in 1:k) {
  if (i < k) {
    dum <- effects %>% 
      filter(cat == i)
    probs[ , i] <- dum$exp_log_odds / denominator
  i} else if (i == k) {
    probs[ , i] <- apply(probs[ , 1:(k-1)], 1, function(x) 1 - sum(x))
  }
}

apply(probs[ , 1:(k-1)], 1, function(x) 1 - sum(x))

#if else doesn't work with mix of numeric and logical vectors
# for (i in 1:k) {
#   dum <- effects %>% 
#     filter(cat == i)
#   dum <- if (k < 1) {
#      dum %>% 
#       mutate(log_odds = int + beta * x)
#   } else {
#     dum %>% 
#       mutate(log_odds = NA)
#   }
# }

prob <- vector(mode = "list", length = k)
for (i in 1:k) {
  if (i < k) {
    prob[[i]] <- exp_log_odds[[i]] / (1 + exp_log_odds[[i]]) 
  }
}

tt <- exp_log_odds %>% 
  sapply(., c) %>% 
  # t() %>% 
  data.frame() %>% 
  glimpse()
apply(tt, 1, sum) + 1

data.frame(t(sapply(mylistlist, c)))

