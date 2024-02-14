rm(list=ls())
############
#Get Data
############
set.seed(525)
mu.true=0 ## these are true parameters that we need to estimate using data Y 
sig2.true=1

n=30 
Y=rnorm(n = n,
        mean = mu.true,
        sd = sqrt(sig2.true)
        )

############
#Set up priors
############
theta=0
tau2=100000 #proper, but not very informative prior for mu; we assume that mean is distributed around 
  # zero, but with huge variance 
  # variance is huge => prior is non informative 

a=0.001
b=0.001  #common noninformative prior for sig2, this will be some gamma 


############
#Initialize Gibbs Sampler
############
nits=10000 # will use 10,000 steps/iterations to get posterior distributions of Mean and Variance of Distribution wehre Y came from 

mu=sig2=rep(NA,nits)

### these are extremely wrong, but because the 
### problem is relatively simple, so this does not matter in the end 

## these are just our point estimate~guesses at what the true parameters are 
mu[1]=10000 #mean(Y)    #seems like logical initial values      
sig2[1]=200 #var(Y)


############
# this is the sampler now 

for(it in 2:nits){
  #########
  #update mu
  ########
  mumean=(tau2*mean(Y) + sig2[it-1]/n * theta)/(tau2 + sig2[it-1]/n) # calculate updated mean of mu distribution 
    # recall: mu is the mean of distribution that generated values of Y 
    # but, in Bayes inference, we need to get distribution of MU - true mean of Y 
    # we need to get distribution of values of MU 
    # this parameter is the mean of distribution of MU 
  
    # how we calculated it?? we use proportional argument and derive by hand what 
    # parameters distribution of mu must have 
    # it depends of observed mean of Y, 
    # prior assumed variance (standard error) of MU distribution,  
    # prior assumed mean of MU distribution, 
  
    # most recent estimate of the variance of Y
  
  muvar =(tau2 * sig2[it-1]/n) / (tau2 + sig2[it-1]/n)
    # now we need to estimate variance of distirbution of MU 
    # again, it is expressed as some formula that we derived by hand 
  
  mu[it]=rnorm(1,mumean,sqrt(muvar))
    # now using updated values of parameters of distribution of MU, get one value of MU, and save it 
    # these values will eventually form posterior sampling distributions of MU, which we will do inference on 
  
  #########
  #update sig2
  #########
  
  # same idea, no reason to comment on this 
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

dev.off()

## after filtering out bad burn in values, we get a symmetric normal distribution of the parameter MU: 
#     which goverens distirubiton of values of Y 
mu_f <- mu[-burnin]
hist(mu_f)  
summary(mu_f)
sd(mu_f)
mean(mu_f)
median(mu_f)

sig2_f <- sig2[-burnin]
hist(sig2_f)  
summary(sig2_f)
mean(sig2_f)
median(sig2_f)

### recall: true values of MU and SIGMA for distirbution of Y were 0 and 1 respectively 
### this is what we get here, which is great. 

##### In slides, week2_lecture1: page 22 shows how we get formulas for updating mu and sigma values 

