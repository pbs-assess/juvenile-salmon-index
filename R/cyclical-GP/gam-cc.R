# based on https://github.com/skaug/tmb-case-studies

library(TMB)
library(mgcv) #Use gam
library(Matrix) #Use sparse matrices
# 
# set.seed(2) ## simulate some data...
# dat <- gamSim(1, n = 400, dist="normal", scale=2)
# plot(dat$x2, dat$y)
# 
# m <- gam(y ~ s(x2, bs = "cs"), data = dat)
# plot(m)
# 
# dat$x2 <- dat$x2 * 100
# m <- gam(y ~ s(x2, bs = "cc", k = 12), data = dat)
# plot(m)

# https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
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


# Tutorial: sparse matrixes in R and TMB
M = matrix(1:4,2,2)   # ordinary 2x2 matrix
M_block_diag = .bdiag(list(M,M)) # Block diagonal (sparse) matrix
data.class(M_block_diag)     # Check data.class
print(M_block_diag)    # dots means 0 value

compile("pSplines.cpp")
dyn.load(dynlib("pSplines"))
#----------------------------------------

#Set up spline structure by using mgcv---
gam_setup = gam(Temperature ~ s(nMonth, bs = "cs"),
  data = cet, fit=FALSE)

gam_setup = gam(Temperature ~ s(nMonth, bs = "cs"),
  data = cet, fit=FALSE)

S1 = gam_setup$smooth[[1]]$S[[1]]
S_list = list(S1)
S_combined = .bdiag(S_list)         # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S

x2 = seq(min(cet$nMonth),max(cet$nMonth), length.out = 100)

x2Report = PredictMat(gam_setup$smooth[[1]],data = data.frame(x2))

designMatrixForReport = list(x2Report)

data = list(Y = dat$y, # Response
  X = gam_setup$X[,-1],  # Design matrix, without intercept
  S = S_combined,      # Combined penalty matrix
  Sdims = Sdims,
  designMatrixForReport = .bdiag(designMatrixForReport),
  flag = 1)

par = list(
  beta0 = 0,  # Intercept
  beta = rep(0,sum(Sdims)),  # Spline coefficients
  log_lambda = rep(rep(0,length(Sdims))), #Log spline penalization coefficients
  log_sigma = 0
)

obj = MakeADFun(data = data, parameters = par,random="beta",DLL = "pSplines")
opt = nlminb(obj$par,obj$fn,obj$gr)
rep = sdreport(obj)

rep

muSpline = rep$value[names(rep$value)=="splineForReport"]
sdSpline<-rep$sd[names(rep$value)=="splineForReport"]

plot(x2, muSpline, lty=1,type = 'l',ylim = range(dat$y - mean(dat$y)),ylab = "f(x2)",main = "Spline for x2")
lines(x2, muSpline- 1.96*sdSpline, lty=2)
lines(x2, muSpline+ 1.96*sdSpline, lty=2)
points(dat$x2, dat$y - mean(dat$y), col = "#00000070")
abline(h = 0)
