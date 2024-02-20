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

q=c(.6,.6,2)

############
#Initialize Gibbs Sampler
############
nits=100000
beta0=beta1=sig2=rep(NA,nits)
beta0[1]=0  #why not?
beta1[1]=0  #ditto
sig2[1]=var(Y)

#############
#Gibbs Sampler
#############
for(it in 2:nits){
  #########
  #update beta0
  ########
  b0s=rnorm(1,beta0[it-1],q[1])
  ra=sum( (Y-b0s-beta1[it-1]*x)^2 ) - sum( (Y-beta0[it-1]-beta1[it-1]*x)^2 )
  rb=b0s^2 - beta0[it-1]^2
  rc=(beta0[it-1] - b0s)^2 - (b0s - beta0[it-1])^2
  r=exp(-ra/(2*sig2[it-1]) - rb/(2*tau) - rc/(2*q[1]))
  beta0[it]=ifelse(r>runif(1),b0s,beta0[it-1])

  #########
  #update beta1
  ########
  b1s=rnorm(1,beta1[it-1],q[2])
  ra=sum( (Y-beta0[it]-b1s*x)^2 ) - sum( (Y-beta0[it]-beta1[it-1]*x)^2 )
  rb=b1s^2 - beta1[it-1]^2
  rc=(beta1[it-1] - b1s)^2 - (b1s - beta1[it-1])^2
  r=exp(-ra/(2*sig2[it-1]) - rb/(2*tau) - rc/(2*q[2]))
  beta1[it]=ifelse(r>runif(1),b1s,beta1[it-1])

  #########
  #update sig2
  #########
  sigs=1/rgamma(1,q[3],q[3]*sig2[it-1])
  ra=(sigs/sig2[it-1])^(-(n/2+a)-1 + 2*q[3]+1)
  rb=sum( (Y-beta0[it]-beta1[it]*x)^2 ) * (1/sigs - 1/sig2[it-1])
  rc=b*(1/sigs - 1/sig2[it-1])
  rd=q[3]*(sigs/sig2[it-1] - sig2[it-1]/sigs)
  r=ra * exp(-rb/2 - rc - rd)
  sig2[it]=ifelse(r>runif(1),sigs,sig2[it-1])
}


fit=lm(Y~x)
dev.new()
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

apply(cbind(beta0,beta1,sig2)[-1,]!=cbind(beta0,beta1,sig2)[-nits,],2,mean)
apply(cbind(beta0,beta1,sig2),2,quantile,c(.5,.025,.975))
#From Gibbs
#            beta0     beta1      sig2
#50%    0.02552855 1.4499046 1.5441631
#2.5%  -0.43361556 0.9792423 0.9542619
#97.5%  0.48034502 1.9189622 2.7611493
