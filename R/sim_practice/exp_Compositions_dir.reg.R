library(dirmult)
library(TMB)
library(tidyverse)
library(Compositional)

xx <- as.vector(iris[, 4:5])
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
pred_x <- iris[, 4:5] %>% 
  distinct() %>% 
  arrange(Petal.Width) 

mod1 <- diri.reg(y, xx, xnew = pred_x)
# mod2 <- diri.reg2(y, xx)

# guts of Compositions::dir.reg function
dm <- dim(y)
n <- dm[1]
x <- model.matrix(y ~ ., data.frame(xx))
d <- dm[2] - 1
z <- Rfast::Log(y)

runtime <- proc.time()
rla <- z[, -1] - z[, 1]
ini <- as.vector(solve(crossprod(x), crossprod(x, rla)))

param <- c(3, ini)
# dirireg <- function(param, z, x, n, d) {
  phi <- exp(param[1])
  para <- param[-1]
  be <- matrix(para, ncol = d)
  mu1 <- cbind(1, exp(x %*% be))
  ma <- mu1/rowSums(mu1)
  ba <- phi * ma
  -n * lgamma(phi) + sum(lgamma(ba)) - sum(z * (ba - 1))
# }

kk <- ncol(y)
mu1b <- matrix(NA, nrow = n, ncol = kk)
exp_eff <- exp(x %*% be)
for (i in 1:n) {
  mu1b[i, 1] <- 1
  for (k in 1:d) {
    mu1b[i, (k + 1)] <- exp_eff[i, k]
  }
}
  
compile(here::here("src", "dirichlet_fixInt_v2.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_fixInt_v2")))

cov_mat <- model.matrix(~., data = data.frame(xx))
obj <- MakeADFun(data=list(fx_cov = cov_mat,
                           y_obs = y,
                           pred_cov = model.matrix(~., data = pred_x)
                           ),
                 parameters=list(z_ints = matrix(ini, ncol = d),
                                 log_phi = 3),
                 DLL="dirichlet_fixInt_v2")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

pp <- matrix(ssdr[rownames(ssdr) %in% "pred_est", ], ncol = d + 1)

matrix<Type> pred_mu1(n_preds, (n_fix_cov + 1));
matrix<Type> pred_est(n_preds, (n_fix_cov + 1));

matrix<Type> pred_eff = pred_cov * z_ints;
matrix<Type> exp_pred_eff = exp(pred_eff.array());

for(int i = 0; i < n_obs; i++) {
  pred_mu1(i, 0) = 1;
  for(int m = 0; m < n_fix_cov; m++) {
    pred_mu1(i, (m + 1)) = exp_pred_eff(i, m);
  }
}

func1 <- 'NumericMatrix mmult1(NumericMatrix exp_pred_eff) {
  int n_obs = exp_fx_eff.nrow();
  int n_fix_cov = exp_fx_eff.ncol():
  NumericMatrix mu1 = no_init_matrix(n_obs, 3);
  
for(int i = 0; i < n_obs; i++) {
  mu1(i, 0) = 1;
  for(int m = 0; m < n_fix_cov; m++) {
    mu1(i, (m + 1)) = exp_pred_eff(i, m);
  }
}
  return(pred_mu1)
}'
