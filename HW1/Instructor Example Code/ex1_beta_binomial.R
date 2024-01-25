
#           tab_save = tab

tab_ex = tab_save
tab_ex=tab_ex[tab_ex$Notes!='Total',]
tab_ex <- tab_ex[!is.na(tab_ex$County.Code), ]

Nr=2
Ns=67

Y=matrix(tab_ex$Deaths,Ns,Nr,byrow=TRUE)
N=matrix(tab_ex$Population,Ns,Nr,byrow=TRUE)

state_rate=sum(Y)/sum(N)

rates=Y/N

high=rates>state_rate
black=rep(1:0,each=Ns)

y=c(sum(rates[,1]>state_rate,na.rm=TRUE), sum(rates[,2]>state_rate,na.rm=TRUE))
n=c(sum(!is.na(rates[,1])), sum(!is.na(rates[,2])))

a=1; b=1

set.seed(1234)
theta_B=rbeta(10000,y[1]+a,n[1]-y[1]+b)
theta_W=rbeta(10000,y[2]+a,n[2]-y[2]+b)

mod=glm(c(high)~c(black),family=binomial)


differ=theta_B-theta_W
logOR=log( theta_B/(1-theta_B) / (theta_W/(1-theta_W)) )

#png('HD_betab_exinomial.png',height=2*5*360,width=2*5*540,res=600)
par(mfrow=c(2,2),cex=1,mar=c(4,4,2,1)+.1)
lims=0:1
hist(theta_B,breaks=seq(lims[1],lims[2],length.out=100),freq=FALSE)
curve(dbeta(x,a,b),from=lims[1],to=lims[2],col=2,add=TRUE)
legend('topright',legend='Prior',lty=1,col=2,cex=1,bty='n')
hist(theta_W,breaks=seq(lims[1],lims[2],length.out=100),freq=FALSE)
curve(dbeta(x,a,b),from=lims[1],to=lims[2],col=2,add=TRUE)
legend('topright',legend='Prior',lty=1,col=2,cex=1,bty='n')
hist(differ,breaks=100,freq=FALSE,main='Histogram of theta_B - theta_W')
hist(logOR,breaks=100,freq=FALSE)
curve(dnorm(x,summary(mod)$coef[2,1],summary(mod)$coef[2,2]),col=2,add=TRUE)
legend('topright',legend='Freq',lty=1,col=2,cex=1,bty='n')
#dev.off()
