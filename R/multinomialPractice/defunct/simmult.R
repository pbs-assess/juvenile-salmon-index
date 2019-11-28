#============================================================================
#Simulation test of TMB algorithm 
#Catarina Wor
# October 2019
#============================================================================
library(ggplot2)
library(TMB)


compile("multinom_test.cpp")
dyn.load(dynlib("multinom_test"))

pk<-c(0.05, 0.05, 0.01, 0.003, 0.4, 0.005, 0.002, 0.2, 0.01, 0.1, 0.16, 0.01)
#c(0.1, 0.1, 0.1, 0.1, 0.4, 0.1, 0.1)

simfunc <- function(pk, n){

	x <- matrix(NA,ncol=length(pk),nrow=n)

	N <- rpois(n,20)

	for(i in seq_len(n)){
		x[i,]<-rmultinom(1,N[i],pk)[,1]
	}
	
	return(x)

}



parameters <- list(
  p=log(rep(1/length(pk),length(pk)))
  )


#==================================
#simeval

nsim<-1000
parest<-matrix(NA,nrow=nsim,ncol=length(pk))

for(a in seq_len(nsim)){


	data <- list(X=simfunc(pk, 30),nobs=30)
	obj <- MakeADFun(data,parameters,DLL="multinom_test")	
	opt<-nlminb(obj$par,obj$fn,obj$gr)

	#
	parest[a,] <- exp(opt$par)/sum(exp(opt$par))
}

df<-reshape::melt(parest)

dftrue<-data.frame(value=pk,
	X2=seq_along(pk))
 
p <-ggplot(df)
p <- p + geom_density(aes(value))
p <- p + facet_wrap(~ X2, scales="free")
p <- p + theme_bw()
p <- p + geom_vline(data=dftrue,aes(xintercept=value))
p

#seems







