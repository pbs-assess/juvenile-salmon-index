library(mgcv)
library(tidyverse)

set.seed(123)
n <- 12
x <- seq(1/n, 1, length.out = n)
x2 <- rnorm(12, 0, 1)#seq(1, 12, length.out = n)
y1 <- sin(x*2*pi) + rnorm(n)*.2 + (0.1 * x2)
y2 <- sin(x*2*pi) + rnorm(n)*.2 + (0.1 * x2)
d <- data.frame(
  x = rep(x * n, 2), 
  x2 = rep(x2, 2),
  y = c(y1, y2), 
  f = rep(c("a", "b"), each = 12)
  ) %>% 
  mutate(y_2 = sin(x*2*pi) + rnorm(n * 2)*.2)
d <- d[-c(16, 18, 20), ]
d$f <- as.factor(d$f)
d_sub <- d[d$f == "b", ]

ggplot(d) +
  geom_point(aes(x = x, y = y, col = f))

kk <- 5

# One factor level
# use predict, type = lpmatrix, not gam fit = false
# m1 <- gam(rep(0, length.out = nrow(nd)) ~ s(x, bs = "tp", k = 6),
#               knots = list(x = c(1, 12)), 
#               data=nd, 
#           fit = FALSE)
m1_fit <- gam(y ~ s(x, bs = "tp", k = 6) + x2,
              knots = list(x = c(1, 12)), 
              data=d_sub)
nd <- expand.grid(x = seq(1, 12, length.out = 50),
                  x2 = mean(d$x2))

Xp <- predict(m1_fit, nd, type = "lpmatrix")
p <- Xp %*% coef(m1_fit)
p2 <- predict(m1_fit, nd)
# p3 <- m1$X %*% coef(m1_fit)
plot(d_sub$x, d_sub$y)
lines(nd$x, p)
lines(nd$x, p2, col = "red")
# lines(nd$x, p3, col = "blue")


# Multiple factor levels
m2_fit <- gam(y ~ s(x, bs = "tp", k = 6, by = f) + x2,
              knots = list(x = c(1, 12)), 
              data=d)
nd2 <- expand.grid(
  x = seq(1, 12, length.out = 50),
  x2 = mean(d$x2),
  f = unique(d$f)
  )

Xp <- predict(m2_fit, nd2, type = "lpmatrix")
nd2$p <- as.vector(Xp %*% coef(m1_fit))
nd2$p2 <- as.vector(predict(m2_fit, nd2))

ggplot() +
  geom_point(data = d, aes(x = x, y = y, colour = f)) +
  geom_line(data = nd2, aes(x = x, y = p, colour = f)) +
  geom_line(data = nd2, aes(x = x, y = p2, colour = f), linetype = 2)
  


# Use smoothcon and list (INCOMPLETE AFTER SWITCHING BACK TO GAM)
sm <- smoothCon(s(x, bs = "tp", k = 4), knots = list(x = c(1, 12)), 
                data=d_sub)[[1]]
gam_sm <- gam(y ~ s(x, bs = "tp", k = 4) + x2, 
              knots = list(x = c(1, 12)), 
              data = d_sub)
y <- d_sub$y
beta <- coef(lm(y ~ sm$X + d_sub$x2 + 0))
nd <- data.frame(x = seq(1, 12, length.out = 200),
                 x2 = 0)
Xp <- PredictMat(sm, nd)
p <- Xp %*% beta[1:4] + (nd$x2 * beta[5]) 
p2 <- predict(gam_sm, nd)

plot(d$x, d$y)
lines(nd$x, p)
lines(nd$x, p2, col = "red")

# incorporate factors
d_sub2 <- d[-c(3, 11), ]
y <- d_sub2$y

sm <- smoothCon(s(x, bs = "tp", k = 7, by = f), knots = list(x = c(1, 12)), 
                data=d_sub2)

# combine 
tt <- map(sm, function(x) {
  dum <- x$X
  labs <- paste(x$label, seq(1, x$bs.dim, by = 1), sep = "_")
  colnames(dum) <- labs
  dum
}) %>% 
  do.call(cbind, .)
betaT <- coef(lm(y ~ tt + 0))

nd <- data.frame(x = rep(seq(1, 12, length.out = 100), 2),
                 f = rep(c("a", "b"), each = 100))
nd_list <- split(nd, nd$f)

## combine each 
out <- vector(2, mode = "list")
#number of splines (necessary to pad matrix)
n_spline <- sm[[1]]$bs.dim
n_fac <- length(sm)
#can't use map because predictmat requires the full smoothcon list, but map 
#requires lists of equal length
for (i in seq_len(n_fac)) {
  #make empty list that will represent 0s and pred matrix for one fac level
  int_list <- vector(n_fac, mode = "list")
  for (j in seq_len(n_fac)) {
    if (i == j) {
      # add prediction matrix
      int_list[[j]] <- PredictMat(sm[[i]], nd_list[[i]])  
    } else {
      # otherwise pad with 0s
      int_list[[j]] <- matrix(0, nrow = nrow(nd_list[[i]]), ncol = n_spline)
    }
    colnames(int_list[[j]]) <- paste(sm[[i]]$label, 
                                     seq(1, n_spline, by = 1), sep = "_")
  }
  out[[i]] <- do.call(cbind, int_list)
}

ttt <- do.call(rbind, out) 
p <- ttt %*% betaT
plot(d_sub2$x, d_sub2$y)
lines(nd$x[1:100], p[1:100])
lines(nd$x[101:200], p[101:200])

## Unsure how to incorporate covariates...

