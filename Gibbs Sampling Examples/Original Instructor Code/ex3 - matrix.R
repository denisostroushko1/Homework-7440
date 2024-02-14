rm(list=ls())
library(mvtnorm)
############
#Get Data
############
set.seed(525)
beta.true=c(0,1); n=30
x=rnorm(n)
X=cbind(1,x)
sig.true=1
Y=cbind(1,x)%*%beta.true+rnorm(n,0,sqrt(sig.true))

############
#Set up priors
############
a=0.001; b=0.001  #common noninformative prior for sig2

############
#Initialize Gibbs Sampler
############
nits=10000
beta=array(dim=c(2,nits))
sig2=rep(NA,nits)
beta[,1]=0  #why not?
sig2[1]=var(Y)

######. #######
#Gibbs Sampler
#######. ######
for(it in 2:nits){
  #########
  #update beta
  ########
  muvar = sig2[it-1] * solve(t(X)%*%X)
  mumean=solve(t(X)%*%X)%*%t(X)%*%Y
  beta[,it]=rmvnorm(1,mumean,muvar)

  #########
  #update sig2
  #########
  sig2[it]=1/rgamma(1,n/2+a,sum( (Y-X%*%beta[,it])^2 )/2 + b)
}

dev.new()
fit=lm(Y~x)
par(mfrow=c(2,3))
plot(beta[1,],type='l')
abline(h=beta.true[1],col=2); abline(h=fit$coef[1],col=3)
plot(beta[2,],type='l')
abline(h=beta.true[2],col=2); abline(h=fit$coef[2],col=3)
plot(sig2,type='l')
abline(h=sig.true,col=2); abline(h=summary(fit)$sigma^2,col=3)
hist(beta[1,],breaks=100)
abline(v=beta.true[1],col=2); abline(v=fit$coef[1],col=3)
hist(beta[2,],breaks=100)
abline(v=beta.true[2],col=2); abline(v=fit$coef[2],col=3)
hist(sig2,breaks=100)
abline(v=sig.true,col=2); abline(v=summary(fit)$sigma^2,col=3)
