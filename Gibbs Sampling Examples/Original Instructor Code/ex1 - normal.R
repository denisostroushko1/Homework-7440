rm(list=ls())
############
#Get Data
############
set.seed(525)
mu.true=0; sig2.true=1
n=30; Y=rnorm(n,mu.true,sqrt(sig2.true))

############
#Set up priors
############
theta=0; tau2=100000 #proper, but not very informative prior for mu
a=0.001; b=0.001  #common noninformative prior for sig2

############
#Initialize Gibbs Sampler
############
nits=10000
mu=sig2=rep(NA,nits)
mu[1]=10000 #mean(Y)    #seems like logical initial values      
sig2[1]=200 #var(Y)

#############
#Gibbs Sampler
#############
for(it in 2:nits){
  #########
  #update mu
  ########
  mumean=(tau2*mean(Y) + sig2[it-1]/n * theta)/(tau2 + sig2[it-1]/n)
  muvar =(tau2 * sig2[it-1]/n) / (tau2 + sig2[it-1]/n)
  mu[it]=rnorm(1,mumean,sqrt(muvar))

  #########
  #update sig2
  #########
  sig2[it]=1/rgamma(1,n/2+a,sum( (Y-mu[it])^2 )/2 + b)
}

burnin=1:100

par(mfrow=c(2,3))
plot(exp(mu[-burnin]),type='l'); abline(h=exp(mu.true),col=2)
plot(exp(sqrt(sig2[-burnin])),type='l'); abline(h=exp(sqrt(sig2.true)),col=2)
plot(exp(mu[-burnin]+sig2[-burnin]/2),type='l'); abline(h=exp(mu.true+sig2.true/2),col=2)
hist(exp(mu[-burnin]),breaks=100,freq=FALSE)
hist(exp(sqrt(sig2[-burnin])),breaks=100,freq=FALSE)
hist(exp(mu[-burnin]+sig2[-burnin]/2),breaks=100,freq=FALSE)