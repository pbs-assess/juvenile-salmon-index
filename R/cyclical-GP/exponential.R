# GP regression example

GP_sq_exp <- function(x, b0 = 0, sigma = 0.1, gp_sigma = 0.1, gp_theta = 2) {
  Sigma <- sapply(x, function(x1) {
    sapply(x, function(x2) {
      dist_sq <- (x2 - x1)^2
      gp_sigma*gp_sigma * exp(- 1 / (2 * gp_theta * gp_theta) * dist_sq)
    })
  })
  eta <- MASS::mvrnorm(mu = rep(0, length(x)), Sigma = Sigma)
  y <- rnorm(length(x), mean = b0 + eta, sd = sigma)
  list(x = x, y = y, eta = eta)
}

library(TMB)
compile("exponential.cpp")
dyn.load(dynlib("exponential"))

set.seed(108)
out <- GP_sq_exp(1:40, b0 = 0, sigma = 0.2, gp_sigma = 0.4, gp_theta = 1.5)
plot(out$x, out$eta, type = "l")
points(out$x, out$y)

data <- list(x = out$x, y = out$y, flag = 1)
parameters <- list(
  b0 = 0,
  log_gp_sigma = log(0.5), 
  log_gp_theta = log(1),
  log_sigma = log(0.1), 
  eta = rep(0, length(out$x))
)
obj <- TMB::MakeADFun(data, parameters, DLL="exponential",
  random = "eta")
obj <- normalize(obj, flag="flag")
opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
sdr <- TMB::sdreport(obj)
sdr

r <- obj$report()
print(r)

