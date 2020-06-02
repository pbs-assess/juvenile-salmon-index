x <- as.vector(iris[, 4])
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
mod1 <- diri.reg(y, x)
mod2 <- diri.reg2(y, x)



# guts of Compositions::dir.reg function
dm <- dim(y)
n <- dm[1]
x <- model.matrix(y ~ ., data.frame(x))
d <- dm[2] - 1
z <- Rfast::Log(y)

runtime <- proc.time()
rla <- z[, -1] - z[, 1]
ini <- as.vector(solve(crossprod(x), crossprod(x, rla)))
el <- NULL
oop <- options(warn = -1)
on.exit(options(oop))

param <- c(3, ini)
dirireg <- function(param, z, x, n, d) {
  phi <- exp(param[1])
  para <- param[-1]
  be <- matrix(para, ncol = d)
  mu1 <- cbind(1, exp(x %*% be))
  ma <- mu1/rowSums(mu1)
  ba <- phi * ma
  -n * lgamma(phi) + sum(lgamma(ba)) - sum(z * (ba - 1))
}

qa <- nlm(dirireg, c(3, ini), z = z, x = x, n = n, d = d)
el1 <- -qa$minimum
qa <- nlm(dirireg, qa$estimate, z = z, x = x, n = n, d = d)
el2 <- -qa$minimum
vim <- 2
while (el2 - el1 > 1e-06) {
  el1 <- el2
  qa <- nlm(dirireg, qa$estimate, z = z, x = x, n = n, 
            d = d)
  el2 <- -qa$minimum
}
qa <- optim(qa$estimate, dirireg, z = z, x = x, n = n, d = d, 
            hessian = TRUE)
log.phi <- qa$par[1]
be <- matrix(qa$par[-1], ncol = d)
colnames(be) <- colnames(y[, -1])
seb <- sqrt(diag(solve(qa$hessian)))
std.logphi <- seb[1]
seb <- matrix(seb[-1], ncol = d)
if (!is.null(colnames(y))) {
  colnames(seb) <- colnames(y[, -1])
}
else colnames(seb) <- paste("Y", 1:d, sep = "")
if (!is.null(xnew)) {
  xnew <- model.matrix(~., data.frame(xnew))
  mu <- cbind(1, exp(xnew %*% be))
  est <- mu/Rfast::rowsums(mu)
  lev <- NULL
}
else {
  mu <- cbind(1, exp(x %*% be))
  est <- mu/Rfast::rowsums(mu)
  lev <- (exp(log.phi) + 1) * Rfast::rowsums((y - est)^2/mu)
}
runtime <- proc.time() - runtime
rownames(be) <- colnames(x)
if (!is.null(seb)) 
  rownames(seb) <- colnames(x)
list(runtime = runtime, loglik = -qa$value, phi = exp(log.phi), 
     log.phi = log.phi, std.logphi = std.logphi, be = be, 
     seb = seb, lev = lev, est = est)


compile(here::here("src", "dirichlet_fixInt_v2.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_fixInt_v2")))

obj <- MakeADFun(data=list(fx_cov = matrix(x, ncol = 2),
                           y_obs = y),
                 parameters=list(z_ints = matrix(ini, ncol = 2),
                                 log_phi = 3),
                 DLL="dirichlet_fixInt_v2")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

