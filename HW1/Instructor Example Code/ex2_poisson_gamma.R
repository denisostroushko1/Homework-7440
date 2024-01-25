
#           tab_save = tab

tab_ex = tab_save
tab_ex=tab_ex[tab_ex$Notes!='Total',]
tab_ex <- tab_ex[!is.na(tab_ex$County.Code), ]

Nr=2
Ns=67

Y=matrix(tab$Deaths,Ns,Nr,byrow=TRUE)
N=matrix(tab$Population,Ns,Nr,byrow=TRUE)

black=rep(1:0,each=Ns)

mod=glm(c(Y)~c(black),family=poisson,offset=log(c(N)),subset=N!=0)

y=c(sum(Y[,1]),sum(Y[,2]))
n=c(sum(N[,1]),sum(N[,2]))

lambda0=sum(Y)/sum(N)
n0=1000; y0=lambda0*n0
#This corresponds to a prior worth 1,000 people

set.seed(1234)
lambda_B=rgamma(10000,y[1]+y0,n[1]+n0)
lambda_W=rgamma(10000,y[2]+y0,n[2]+n0)

differ=lambda_B-lambda_W
logRR=log( lambda_B/lambda_W )

#png('HD_poissongamma.png',height=2*5*360,width=2*5*540,res=600)
par(mfrow=c(2,2),cex=1,mar=c(4,4,2,1)+.1)
lims=quantile(c(lambda_B,lambda_W),0:1)
hist(lambda_B,breaks=seq(lims[1],lims[2],length.out=100),freq=FALSE)
curve(dgamma(x,y0,n0),from=lims[1],to=lims[2],col=2,add=TRUE)
legend('topright',legend='Prior',lty=1,col=2,cex=1,bty='n')
hist(lambda_W,breaks=seq(lims[1],lims[2],length.out=100),freq=FALSE)
curve(dgamma(x,y0,n0),from=lims[1],to=lims[2],col=2,add=TRUE)
legend('topright',legend='Prior',lty=1,col=2,cex=1,bty='n')
hist(differ,breaks=100,freq=FALSE,main='Histogram of lambda_B - lambda_W')
hist(logRR,breaks=100,freq=FALSE)
curve(dnorm(x,summary(mod)$coef[2,1],summary(mod)$coef[2,2]),col=2,add=TRUE)
legend('topright',legend='Freq',lty=1,col=2,cex=1,bty='n')
#dev.off()
