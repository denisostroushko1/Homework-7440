rm(list=ls())
############
#Get Data
############
set.seed(525)
beta.true=c(0,1); n=30
x=rnorm(n)
sig.true=1
Y=cbind(1,x)%*%beta.true+rnorm(n,0,sqrt(sig.true))

############
#Set up priors
############
tau=10000 #big prior variances for beta
a=0.001; b=0.001  #common noninformative prior for sig2

############
#Initialize Gibbs Sampler
############
nits=10000
beta0=beta1=sig2=rep(NA,nits)
beta0[1]=0  #why not? -- assume that there is no relationship 
beta1[1]=0  #ditto
sig2[1]=var(Y) # bayes data driven approach to initiate the algorithm 
 
#############
#Gibbs Sampler
#############
for(it in 2:nits){
  #########
  #update beta0
  ########
  mumean=tau/(tau+sig2[it-1]/n) * mean(Y-x*beta1[it-1]) ## this is a mean for the beta 0 sampling distribution 
  muvar =(tau * sig2[it-1]/n) / (tau + sig2[it-1]/n) ## this is a variance for the beta 0 sampling distribution 
  beta0[it]=rnorm(1,mumean,sqrt(muvar)) # this is one draw from the updated sampling distribution of beta 

  #########
  #update beta1
  ########
  mumean=tau/(tau+sig2[it-1]/sum(x^2)) * sum(x*(Y-beta0[it]))/sum(x^2)
  muvar =(tau * sig2[it-1]/sum(x^2)) / (tau + sig2[it-1]/sum(x^2))
  beta1[it]=rnorm(1,mumean,sqrt(muvar))

  #########
  #update sig2
  #########
  sig2[it]=1/rgamma(1,n/2+a,sum( (Y-beta0[it]-x*beta1[it])^2 )/2 + b)
}

fit=lm(Y~x)
par(mfrow=c(2,3))
plot(beta0,type='l')
abline(h=beta.true[1],col=2); abline(h=fit$coef[1],col=3)
plot(beta1,type='l')
abline(h=beta.true[2],col=2); abline(h=fit$coef[2],col=3)
plot(sig2,type='l')
abline(h=sig.true,col=2); abline(h=summary(fit)$sigma^2,col=3)
hist(beta0,breaks=100)
abline(v=beta.true[1],col=2); abline(v=fit$coef[1],col=3)
hist(beta1,breaks=100)
abline(v=beta.true[2],col=2); abline(v=fit$coef[2],col=3)
hist(sig2,breaks=100)
abline(v=sig.true,col=2); abline(v=summary(fit)$sigma^2,col=3)
