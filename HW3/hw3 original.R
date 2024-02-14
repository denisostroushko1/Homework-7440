#https://wonder.cdc.gov/controller/saved/D140/D34F844
rm(list=ls())
#First we read in the data and define a few things...
stroke=read.table('2016_PA_stroke_total.txt',sep='\t',
                  stringsAsFactors=FALSE,header=TRUE)
Ng=3        #three age groups
alabs=unique(stroke$Age.Group.Code)
Ns=67     #67 counties
clabs=unique(stroke$County)

#Next we organize things a bit...
Y=array(stroke$Deaths,dim=c(Ng,Ns))
n=array(stroke$Population,dim=c(Ng,Ns))

###################
###################
#per part 4, all Y's below 10 have been suppressed
###################
thres=10      #Suppression threshold
dY=!is.na(Y)  #0 for suppressed, 1 for observed
nsupp=apply(!dY,1,sum) #how many suppressed per age
#note: we do not know the true Y's
#      CDC did the suppression, not me

###################
###################
#insert your prior info here
###################
lam0=c(75,250,1000)/100000
n0=      ###############
Y0=      ###############

###################
###################
#initialize your Gibbs sampler here
nsims=10000
lami=array(dim=c(Ng,Ns,nsims))
for(a in 1:Ng){
  lami[a,,1]=    ################
  Y[k,!dY[a,]]=  ################
  #Note: the preceding line assumes we don't care
  #      what the posterior dist of the missing
  #      Y's looks like -- I'm just plugging the
  #      current guesses directly into my data vector
}

for(it in 2:nsims){
  for(a in 1:Ng){
    ###################
    ###################
    #ADDRESS SUPPRESSED Y HERE
    ###################
    ###################

    ###################
    ###################
    #ESTIMATE LAMBDA_{ia} HERE
    ###################
    ###################
  }
}

#################
#################
#Get posterior samples
#of the age-adjusted rates
#################
aalami=array(dim=c(Ns,nsims))
for(i in 1:Ns){
  aalami[i,]=  ##################
}

##################
##################
#calculate the posterior medians
# of the age-adjusted rates
aa.med= ##################

##################
##################
#THE BELOW CODE SHOULD BE LEFT AS-IS!
#IT ASSUMES YOU NAMED
#THE POSTERIOR MEDIANS OF THE AGE-ADJUSTED RATES
#"aamed" USING THE CODE ABOVE,
#AND WILL CREATE A MAP "PAmap.png"
#THAT WILL BE SAVED TO YOUR CURRENT DIRECTORY
##################

load('penn.rdata')
install.packages(c('maptools','RColorBrewer'))
library(maptools)
library(RColorBrewer)
ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(aa.med*100000,1:(ncols-1)/ncols)
tcolb=array(rep(aa.med*100000,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('PAmap.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
    plot(penn,col=cols[tcol],border='lightgray',lwd=.5)
    legend('right',inset=c(-.15,0),xpd=TRUE,
           legend=c(paste(
           c('Below',round(tcuts[-(ncols-1)],0),'Over'),
           c(' ',rep( ' - ',ncols-2),' '),
           c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
           fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
           border='lightgray')
dev.off()