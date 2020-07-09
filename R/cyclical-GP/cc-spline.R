library(mgcv)
library(MASS) ## load for mcycle data.


CET <- url("http://www.metoffice.gov.uk/hadobs/hadcet/cetml1659on.dat")
# writeLines(readLines(CET, n = 10))
cet <- read.table(CET, sep = "", skip = 6, header = TRUE,
  fill = TRUE, na.string = c(-99.99, -99.9))
names(cet) <- c(month.abb, "Annual")
## remove last row of incomplete data
cet <- cet[-nrow(cet), ]
## get rid of the annual too - store for plotting
rn <- as.numeric(rownames(cet))
Years <- rn[1]:rn[length(rn)]
annCET <- data.frame(Temperature = cet[, ncol(cet)],
  Year = Years)
cet <- cet[, -ncol(cet)]
cet <- stack(cet)[,2:1]
names(cet) <- c("Month","Temperature")
cet <- transform(cet, Year = (Year <- rep(Years, times = 12)),
  nMonth = rep(1:12, each = length(Years)),
  Date = as.Date(paste(Year, Month, "15", sep = "-"),
    format = "%Y-%b-%d"))
cet <- cet[with(cet, order(Date)), ]
cet <- transform(cet, Time = as.numeric(Date) / 1000)
plot(Temperature ~ Year, data = annCET, type = "l",
  main = "CET")

m <- gam(Temperature ~ s(nMonth, bs = "cc", k = 12),
  data = cet)
plot(m)



sm <- smoothCon(s(nMonth,bs = "cc", k=12), data=cet,knots=NULL)[[1]]

## use it to fit a regression spline model...
beta <- coef(lm(cet$Temperature~sm$X-1))

nMonths <- seq(1,12, length=200)  ## create prediction times
## Get matrix mapping beta to spline prediction at 'times'
Xp <- PredictMat(sm, data.frame(nMonth = nMonths))
plot(nMonths, Xp%*%beta, type = "l") ## add smooth to plot

bs <- list(
  smoothCon(s(nMonth), data = cet),
  smoothCon(s(Time), data = cet)
)

beta <- c


require(mgcv)
n <- 15
x <- runif(n)
y <- sin(x*2*pi) + rnorm(n)*.2
mod <- gam(y~s(x,bs="cc",k=6),knots=list(x=seq(0,1,length=6)), fit = FALSE)
mm <- model.matrix(mod)
mm <- as.data.frame(mm)
names(mm) <- c("b0", "b1", "b2", "b3", "b4")

m <- lm(y ~ 0 + b0 + b1 + b2 + b3 + b4, data = mm)


x <- seq(0, 1, length.out = 100)
y_fake <- rep(1, length(x))

mod <- gam(y_fake~s(x,bs="cc",k=6),knots=list(x=seq(0,1,length=6)), fit = TRUE)
mm <- model.matrix(mod)
mm <- as.data.frame(mm)
names(mm) <- c("b0", "b1", "b2", "b3", "b4")

p <- predict(m, newdata = mm)
plot(x, p)

n <- 15
x <- runif(n)
y <- sin(x*2*pi) + rnorm(n)*.2
jd <- mgcv::jagam(formula = y~s(x,bs="cc",k=12), family = gaussian(),
  file = tempfile(fileext = ".jags"),
  diagonalize = TRUE)
jd$jags.data$X

