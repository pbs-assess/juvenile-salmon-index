library(mgcv)
library(ggplot2)

n <- 12
x <- seq(1/n, 1, length.out = n)
y1 <- sin(x*2*pi) + rnorm(n)*.2
y2 <- sin(x*2*pi) + rnorm(n)*.2
d <- data.frame(x = rep(x * n, 2), y = c(y1, y2), f = rep(c("a", "b"), each = 12))
d <- d[-c(13:14, 24), ]
d2 <- d[d$f == 2, ]
d$f <- as.factor(d$f)

ggplot(d, aes(x, y, colour = f)) + geom_point() +
  facet_wrap(~f)

# cc with all data
sm <- smoothCon(s(x, bs = "cc"), data=d)[[1]]
y <- d$y
beta <- coef(lm(y ~ sm$X + 0))
nd <- data.frame(x = seq(1, 12, length.out = 200))
Xp <- PredictMat(sm, nd)
p <- Xp %*% beta
par(mfrow = c(2, 1))
matplot(sm$X, type = "l", lty = 1)
plot(d$x, d$y)
lines(nd$x, p)

# cc with 1 and 12 knots; missing some months
sm <- smoothCon(s(x, bs = "cc"), knots = list(x = c(1, 12)), data=d2)[[1]]
y <- d2$y
beta <- coef(lm(y ~ sm$X + 0))
nd <- data.frame(x = seq(1, 12, length.out = 200))
Xp <- PredictMat(sm, nd)
p <- Xp %*% beta
par(mfrow = c(2, 1))
matplot(sm$X, type = "l", lty = 1)
plot(d2$x, d2$y, xlim = c(1, 12))
lines(nd$x, p)

# tp; missing some months
sm <- smoothCon(s(x, bs = "tp", k = 5), data=d2)[[1]]
y <- d2$y
beta <- coef(lm(y ~ sm$X + 0))
nd <- data.frame(x = seq(1, 12, length.out = 200))
Xp <- PredictMat(sm, nd)
p <- Xp %*% beta
par(mfrow = c(2, 1))
matplot(sm$X, type = "l", lty = 1)
plot(d2$x, d2$y, xlim = c(1, 12))
lines(nd$x, p)

# cc; by; all data
sm <- smoothCon(s(x, bs = "tp", by = f, k = 5), data=d)[[1]]
sm$X
y <- d$y
beta <- coef(lm(y ~ sm$X + 0))
nd1 <- data.frame(x = seq(1, 12, length.out = 200), f = "a")
nd2 <- data.frame(x = seq(1, 12, length.out = 200), f = "b")
nd <- rbind(nd1, nd2)
nd$f <- as.factor(nd$f)

Xp <- PredictMat(sm, nd)
nd$p <- Xp %*% beta
ggplot(nd, aes(x, p)) + geom_point() + facet_wrap(~f)


