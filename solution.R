rm(list=ls())
#First we read in the data and define a few things...
load(file='midterm_data.rdata')
list2env(mdata,envir=sys.frame(sys.nframe()))

#################
#Y is a 67x2 array of the number of low-weight births
#  in PA counties by White mothers vs. Black mothers
hist(Y)
#################
#n is a 67x2 array of the number of total births
#  in PA counties by White mothers vs. Black mothers
hist(n)

Ns=dim(Y)[1] #67 counties
Ng=dim(Y)[2] #2 races being considered

###################
#prior specifications
###################
tau2=10000  #beta0~N(0,tau2)
as=bs=0.001   #sig2~IG(as,bs)

###################
#candidate density variances
###################
qt=array(1,dim=c(Ns,Ng))

nsims=10000
beta0=sig2=array(dim=c(Ng,nsims))
theta=array(dim=c(Ns,Ng,nsims))
for(a in 1:Ng){
  beta0[a,1]= 0 #prior mean might not be wise here...
#  beta0[a,1]=log( sum(Y[a,dY[a,]])/sum(n[a,dY[a,]]) )  
  sig2[a,1]=1  #bs/as = 1, so this seems neutral
  for(i in 1:Ns){
    theta[i,a,1]=beta0[a,1] #0
  }
}

for(it in 2:nsims){
  for(a in 1:Ng){

    ############
    #update beta0
    bvar=(Ns/sig2[a,it-1] + 1/tau2)^(-1)
    bmean=bvar * (sum(theta[,a,it-1])/sig2[a,it-1])
    beta0[a,it]=rnorm(1,bmean,sqrt(bvar))

    ############
    #update theta
    for(i in 1:Ns){
      ts=rnorm(1,theta[i,a,it-1],qt[i,a])
      ra=Y[i,a] * (ts - theta[i,a,it-1])
      rb=(1+exp(theta[i,a,it-1]))/(1+exp(ts))
      rc=((ts-beta0[a,it])^2 - (theta[i,a,it-1]-beta0[a,it])^2)
      r=exp(ra + n[i,a]*log(rb) - rc/(2*sig2[a,it-1]))
      theta[i,a,it]=ifelse(r>runif(1),ts,theta[i,a,it-1])
    }

    ############
    #update sig2
    sig2[a,it]=1/rgamma(1,Ns/2 + as, sum((theta[,a,it]-beta0[a,it])^2)/2 + bs)
  }
  if(it/100 == floor(it/100)){
    cat('Iteration:',it,'\n')
    for(a in 1:Ng){
      for(i in 1:Ns){
        trate=mean(theta[i,a,it-100+1:99]!=theta[i,a,it-100+2:100])
        trate=ifelse(trate>.75,.75,ifelse(trate<.2,.2,trate))
        qt[i,a]=qt[i,a]*trate/0.44
      }
    }
  }
}

##################
#Problem 3
##################
par(mfrow=c(1,2))
matplot(t(beta0),type='l')
matplot(t(sig2),type='l')
cat("Some burn-in is required;
I'm going to use 2,000 because that's what I used before...\n")
burnin=1:2000

##################
#Problem 4
##################
logOR=beta0[2,]-beta0[1,]
par(mfrow=c(1,1))
hist(logOR[-burnin],breaks=100)
cat("The entire histogram above is greater than 0,
so we have strong evidence of a racial disparity\n")

##################
#Problem 5
##################
pii=exp(theta)/(1+exp(theta))
pci=apply(pii[,,-burnin],1:2,quantile,c(.5,.025,.975))

load('penn.rdata')
#install.packages(c('maptools','RColorBrewer'))
library(maptools)
library(RColorBrewer)
ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]

for(a in 1:Ng){
tcuts=quantile(pci[1,,a],1:(ncols-1)/ncols)*100
tcolb=array(rep(pci[1,,a]*100,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png(paste('PAmap_',a,'.png',sep=''),height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
    plot(penn,col=cols[tcol],border='lightgray',lwd=.5)
    legend('right',inset=c(-.15,0),xpd=TRUE,
           legend=c(paste(
           c('Below',round(tcuts[-(ncols-1)],2),'Over'),
           c(' ',rep( ' - ',ncols-2),' '),
           c(round(tcuts,2),round(tcuts[ncols-1],2)),sep='')),
           fill=cols,title='LW Births per 100',bty='n',cex=1.5,
           border='lightgray')
dev.off()
}

rat=pii[,2,]/pii[,1,]
rci=apply(rat[,-burnin],1,quantile,c(.5,.025,.975))
tcuts=quantile(rci[1,],1:(ncols-1)/ncols)
tcolb=array(rep(rci[1,],each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('PAmap_BW_disparity.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
    plot(penn,col=cols[tcol],border='lightgray',lwd=.5)
    legend('right',inset=c(-.15,0),xpd=TRUE,
           legend=c(paste(
           c('Below',round(tcuts[-(ncols-1)],2),'Over'),
           c(' ',rep( ' - ',ncols-2),' '),
           c(round(tcuts,2),round(tcuts[ncols-1],2)),sep='')),
           fill=cols,title='B/W Disparity',bty='n',cex=1.5,
           border='lightgray')
dev.off()

##################
#Problem 6
##################
par(mfrow=c(1,2))
for(i in c(51,57)){
  hist(rat[i,-burnin],breaks=100,main=paste(penn$NAME[i],'County'),
       xlab='B/W Disparity')
  abline(v=(Y[i,2]/n[i,2])/(Y[i,1]/n[i,1]),col=2)
  abline(v=(sum(Y[,2])/sum(n[,2]))/(sum(Y[,1])/sum(n[,1])),col=4)
}

##################
#Bonus / Model Informativeness
##################
b.med=apply(beta0[,-burnin],1,median)
p0=exp(b.med)/(1+exp(b.med))

info=(1+exp(beta0))/sig2-exp(beta0)/(1+exp(beta0))
ici=apply(info[,-burnin],1,quantile,c(.5,.025,.975))

par(mfrow=c(1,2))
for(k in 1:Ng){
plot(pci[1,,k]/(pci[3,,k]-pci[2,,k]),x=Y[,k])
#alpha=0.5; a=b=0
a=ici[1,k]
b=a*(1-p0[k])/p0[k]
curve(qbeta(.5,x+a,x*(1/p0[k]-1)+b)/
     (qbeta(.975,x+a,x*(1/p0[k]-1)+b)-
      qbeta(.025,x+a,x*(1/p0[k]-1)+b)),
        from=1,to=max(Y[,k]),col=2,lty=1,add=TRUE)

}