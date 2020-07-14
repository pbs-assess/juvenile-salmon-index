library(mgcv)
library(tidyverse)

n <- 12
x <- seq(1/n, 1, length.out = n)
y1 <- sin(x*2*pi) + rnorm(n)*.2
y2 <- sin(x*2*pi) + rnorm(n)*.2
d <- data.frame(x = rep(x * n, 2), y = c(y1, y2), f = rep(c("a", "b"), each = 12)) %>% 
  mutate(y_2 = sin(x*2*pi) + rnorm(n * 2)*.2)
d <- d[-c(16, 18, 20), ]
d$f <- as.factor(d$f)
d_sub <- d[d$f == "b", ]

kk <- 6

## 1) Are splines sensitive to response?
m1 <- gam(y ~ s(x, bs = "tp", by = f, k = kk), data=d, fit = F)
m2 <- gam(y_2 ~ s(x, bs = "tp", by = f, k = kk), data=d, fit = F)

identical(m1$X, m2$X)

# Not sensitive


## 2) Are splines sensitive to explanatory?
n_2 <- n * 2
x_2 <- seq(1/n_2, 1, length.out = n_2)
y1_2 <- sin(x_2*2*pi) + rnorm(n_2)*.2
y2_2 <- sin(x_2*2*pi) + rnorm(n_2)*.2
d2 <-  data.frame(x = rep(x_2 * n, 2),
                  y = c(y1_2, y2_2),
                  f = rep(c("a", "b"), each = n_2))

m1 <- gam(y ~ s(x, bs = "tp", by = f, k = kk), 
          knots = list(x = c(1, 12)), 
          data=d, fit = F)
m2 <- gam(y ~ s(x, bs = "tp",  by = f, k = kk), 
          knots = list(x = c(1, 12)), 
          data=d2, fit = F)

m1$X[1, ]
m2$X[1, ]
# different

m1_fit <- gam(y ~ s(x, bs = "tp", by = f, 
                    k = kk), knots = list(x = c(1, 12)), 
              data = d)
nd1 <- d %>% select(-y)
nd1$pred_y <- m1$X %*% coef(m1_fit) #predict.gam(m1_fit, nd1)
# nd2 <- d2 %>% select(-y)
# nd2$pred_y <- m2$X %*% coef(m1_fit) #predict.gam(m1_fit, nd2)

ggplot() +
  geom_point(data = d, aes(x, y, colour = f)) + 
  geom_line(data = nd1, aes(x, pred_y, colour = f))
# ggplot() +
#   geom_point(data = d2, aes(x, y, colour = f)) + 
#   geom_line(data = nd2, aes(x, pred_y, colour = f))

# multiplying GAM model matrix by predictors non-sensical with gams

# try predicting across gaps
nd3 <- expand.grid(x = seq(1, 12, length.out = 50), 
                   f = unique(d$f), 
                   y = 0)
m3 <- gam(y ~ s(x, bs = "tp", by = f, k = kk), knots = list(x = c(1, 12)), 
          data = nd3, fit = F)
nd3$pred_y <- m3$X %*% coef(m1_fit)

ggplot() +
  geom_point(data = d, aes(x, y, colour = f)) + 
  geom_line(data = nd3, aes(x, pred_y, colour = f))


# as above but with smoothcon and list
sm <- smoothCon(s(x, bs = "tp", k = 4), knots = list(x = c(1, 12)), 
                data=d_sub)[[1]]
y <- d_sub$y
beta <- coef(lm(y ~ sm$X + 0))
nd <- data.frame(x = seq(1, 12, length.out = 200))
Xp <- PredictMat(sm, nd)
p <- Xp %*% beta
# par(mfrow = c(2, 1))
# matplot(sm$X, type = "l", lty = 1)
plot(d$x, d$y)
lines(nd$x, p)


# incorporate factors
d_sub2 <- d[-c(3, 11), ]

sm <- smoothCon(s(x, bs = "tp", k = 4, by = f), knots = list(x = c(1, 12)), 
                data=d_sub2)

tt <- map(sm, function(x) {
  dum <- x$X
  labs <- paste(x$label, seq(1, x$bs.dim, by = 1), sep = "_")
  colnames(dum) <- labs
  as.data.frame(dum)
}) %>% 
  bind_cols() %>% 
  as.matrix()
coef(lm(y ~ tt + 0))


nd <- data.frame(x = rep(seq(1, 12, length.out = 100), 2),
                 f = rep(c("a", "b"), each = 100))
nd_list <- split(nd, nd$f)
ttt <- map2(sm, nd_list, function (x, y) {
  PredictMat(x, y)
  }) %>% 
  bind_rows()

y <- d_sub2$y
beta <- coef(lm(y ~ sm[[2]]$X + 0))
Xp <- PredictMat(sm[[2]], nd)
p <- Xp %*% beta
# par(mfrow = c(2, 1))
# matplot(sm$X, type = "l", lty = 1)
plot(d_sub2$x, d_sub2$y)
lines(nd$x, p)
