rm(list=ls())
set.seed(1234)
I=10; J=5
beta.true=c(-1,1)
x=array(rnorm(I*J),dim=c(I,J))
sig=1
z.true=rnorm(I,0,sqrt(sig))
thetaij=beta.true[1] + beta.true[2]*x + rep(z.true,times=J)

Y=array(rbinom(I*J,1,exp(thetaij)/(1+exp(thetaij))),dim=c(I,J))

#############
#priors
tau=c(10000,10000)
#a=0.001; b=0.001
a=2/2; b=2/2 * 1

#############
#candidate variances
qb=c(0.9,0.9)
qz=rep(0.9,I)

#############
#initial values
nits=50000
beta0=beta1=sig2=rep(NA,nits)
beta0[1]=0
beta1[1]=0
sig2[1]=1
z=array(dim=c(I,nits))
z[,1]=0

##############
#Gibbs sampler
for(it in 2:nits){
  #beta0
  b0s=rnorm(1,beta0[it-1],qb[1])
  ra=sum(Y)*(b0s-beta0[it-1])
  rb=prod( (1+exp(beta0[it-1]+beta1[it-1]*x+rep(z[,it-1],times=J)))/
           (1+exp(b0s+beta1[it-1]*x+rep(z[,it-1],times=J))) )
  rc=(b0s^2-beta0[it-1]^2)
  r=rb * exp(ra - rc/(2*tau[1]))
  beta0[it]=ifelse(r>runif(1),b0s,beta0[it-1])

  #beta1
  b1s=rnorm(1,beta1[it-1],qb[2])
  ra=sum(Y * x) * (b1s - beta1[it-1])
  rb=prod( (1+exp(beta0[it]+beta1[it-1]*x+rep(z[,it-1],times=J)))/
           (1+exp(beta0[it]+b1s*x+rep(z[,it-1],times=J))) )
  rc=(b1s^2-beta1[it-1]^2)
  r=rb * exp(ra - rc/(2*tau[2]))
  beta1[it]=ifelse(r>runif(1),b1s,beta1[it-1])

  #z
  for(i in 1:I){
    zs=rnorm(1,z[i,it-1],qz[i])
    ra=sum(Y[i,])*(zs-z[i,it-1])
    rb=prod( (1+exp(beta0[it-1]+beta1[it]*x[i,]+z[i,it-1]))/
             (1+exp(beta0[it-1]+beta1[it]*x[i,]+zs)) )
    rc=(zs^2-z[i,it-1]^2)
    r=rb * exp(ra - rc/(2*sig2[it-1]))
    z[i,it]=ifelse(r>runif(1),zs,z[i,it-1])
  }

  #sig2
  sig2[it]=1/rgamma(1,I/2+a,sum(z[,it]^2)/2+b)
}

pars=cbind(beta0,beta1,t(z),sig2)
apply(pars[-1,]!=pars[-nits,],2,mean)
apply(pars,2,quantile,c(.5,.025,.975))

par(mfrow=c(2,6))
plot(beta0,type='l')
plot(beta1,type='l')
for(i in 1:I){
  plot(z[i,],type='l')
}

par(mfrow=c(2,6))
hist(beta0,breaks=100)
hist(beta1,breaks=100)
for(i in 1:I){
  hist(z[i,],breaks=100); abline(v=z.true[i],col=2)
}

