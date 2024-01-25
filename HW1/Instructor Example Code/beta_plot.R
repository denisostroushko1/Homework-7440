rm(list=ls())

a=c(1/2,1,1/2,1,2)
b=c(1/2,1,1,1/2,2)
library(RColorBrewer)
cols=brewer.pal(5,'Dark2')

png('beta.png',height=5*480,width=5*540,res=600)
par(cex=1,mar=c(4,4,1,1)+.1)
curve(dbeta(x,a[1],b[1]),from=0,to=1,col=cols[1],ylim=c(0,3),
      xlab='theta ~ Beta(a,b)',ylab='density')
for(k in 2:5){
  curve(dbeta(x,a[k],b[k]),from=0,to=1,col=cols[k],add=TRUE,lty=k)
}
legend('top',legend=paste('a=',format(a,nsmall=1),
                        ', b=',format(b,nsmall=1),sep=''),lty=1:5,
       col=cols,bty='n',cex=2/3)
dev.off()