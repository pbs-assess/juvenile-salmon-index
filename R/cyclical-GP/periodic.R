# periodic GP regression example

# See:
# Roberts, S., Osborne, M., Ebden, M., Reece, S., Gibson, N., and Aigrain, S. 2013. Gaussian processes for time-series modelling. Proc. R. Soc. A 371(1984): 20110550. doi:10.1098/rsta.2011.0550.
# Eq. 20

# Not able to get TMB model to work!

GP_periodic <- function(x, b0 = 0, sigma = 0.1,
  gp_sigma = 0.2, gp_theta = 1.5, gp_T = 12) {
  Sigma <- sapply(x, function(x1) {
    sapply(x, function(x2) {
      gp_part1 <- pi * abs((x2 - x1) / gp_T)
      gp_part2 <- sin(gp_part1) ^ 2
      gp_sigma * gp_sigma * exp((- 1 / (2 * gp_theta * gp_theta)) * gp_part2)
    })
  })
  eta <- MASS::mvrnorm(mu = rep(0, length(x)), Sigma = Sigma)
  y <- rnorm(length(x), mean = b0 + eta, sd = sigma)
  list(x = x, y = y, eta = eta, Sigma = Sigma)
}

library(TMB)
compile("periodic.cpp")
dyn.load(dynlib("periodic"))

set.seed(1018)
out <- GP_periodic(1:20, sigma = 0.1, gp_sigma = 0.5, gp_theta = 0.1, gp_T = 6)
plot(out$x, out$eta, type = "l", ylim = range(union(out$y, out$eta)))
points(out$x, out$y)
abline(v = seq(1, 36, 6))

data <- list(x = out$x, y = out$y, gp_T = 6, flag = 1, Sigma = out$Sigma)
parameters <- list(
  b0 = 0,
  log_gp_sigma = log(0.5), 
  log_gp_theta = log(0.9),
  log_sigma = log(0.1), 
  eta = rep(0, length(out$x))
)
obj <- TMB::MakeADFun(data, parameters, DLL="periodic",
  random = "eta")
# obj <- normalize(obj, flag="flag")
opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
sdr <- TMB::sdreport(obj)
sdr

r <- obj$report()
print(r)

