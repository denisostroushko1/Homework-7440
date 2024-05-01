
rm(list = ls())

                        ##############################
                        # setup_mn.R file contains 
                        ##############################

#####################
#Set Things Up
#####################

                        ##############################
                        # setup_mn.R file contains 
                        ##############################
# source('./Final/HQ Code/Code/Setup/setup_packages.r')
                        
if (!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if (!require(mvtnorm)){
  install.packages("mvtnorm")
  library(mvtnorm)
}

if (!require(MCMCpack)){
  install.packages("MCMCpack")
  library(MCMCpack)
}

if (!require(spdep)){
  install.packages("spdep")
  library(spdep)
}

if (!require(sp)){
  install.packages("sp")
  library(sp)
}

if (!require(arm)){
  install.packages("arm")
  library(arm)
}

if (!require(raster)){
  install.packages("raster")
  library(raster)
}

if (!require(maps)){
    install.packages("maps")
    library(maps)
}

#if (!require(rgeos)){
#  install.packages("rgeos")
#  library(rgeos)
#}

print("")


library(spdep)      # spatial dependencies functions
# library(rgdal)      # manages spatial projections
library(maps)   # general functions for map/spatial manipulation
library(sp)         # foundational definition of spatial classes in R
library(arm)        # Loads R2WinBugs, Lattice, Matrix, others
library(raster)
# library(rgeos)
library(tidyverse)
getwd()

load(file='./Final/HQ Code/Code/Shapefiles/counties.rdata')
load(file='./Final/HQ Code/Code/Shapefiles/states.rdata')

us0=spTransform(us,CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 
     +x_0=0.0 +y_0=0 +k=1.0 +units=m +no_defs"))

mn=us0[us0$STATE_FIPS=='27',]

## Creating adjacency matrix
mn.wt1<-spdep::poly2nb(mn,snap=1000)
#snap=1000 is used because *for this shapefile* some neighbor pairs
#  were not being recognized; this value of "snap" seems to fix that
#By default, poly2nb uses a "queen" definition of contiguity
#  (i.e., sharing a single *point* is sufficient to be a neighbor);
#  this can be changed by using the argument "queen=FALSE"
#  (which should require neighbors to share an *edge*)

mn.WBwt<-spdep::nb2WB(mn.wt1)
adj<-mn.WBwt$adj
num<-mn.WBwt$num; sum(num==0)

## Define Variables
Ns =length(mn)

neigh=list()
count=0
for(i in 1:Ns){
  if(num[i]>0){
    neigh[[i]]=adj[count+1:num[i]]
    count=count+num[i]
  }
}
m=num
                                            
                                            #####################
                                            # hq_final.R contents 
                                            #####################

library(mvtnorm)
library(MCMCpack)


tab = read.table(file = "./Final/HQ Code/Code/Data/mn_suicide.txt", sep='\t',header=TRUE, stringsAsFactors=FALSE)


##################
#Get race, gender, and county labels
##################
rlabs=unique(tab$Race)
R=length(rlabs)-1; rlabs=rlabs[1:R]
slabs=unique(tab$Gender)
S=length(slabs)-1; slabs=slabs[1:S]
clabs=unique(tab$County)
I=length(clabs)-1; clabs=clabs[1:I]

##################
#Throw out "extra" total rows
##################
tab=tab[-dim(tab)[1],]		#last row is the overall total
tab=tab[tab$Gender!="",]	#throw out the race totals

##################
#Set up Y, n, Ytot, Yobs
##################
Y=array(as.numeric(tab$Deaths),dim=c(I+1,S,R))
n=array(as.numeric(tab$Population),dim=c(I+1,S,R))
Ytot=Y[I+1,,]; ntot=n[I+1,,]
Y=Y[-(I+1),,]; n=n[-(I+1),,]
Yobs=apply(Y,2:3,sum,na.rm=TRUE)
nYmiss=Ytot-Yobs

###################
#Missing Y's
###################
dY=!is.na(Y)
nsupp=apply(!dY,2:3,sum)
Ythres=c(1,9)

###################
#prior specifications
###################
gam2=10000		#beta0~N(0,gam2)
as=1; bs=1/7	#sig2~IG(as,bs)
at=1; bt=1/100	#tau2~IG(at,bt)

###################
#candidate density variances
###################


qt=array(1,dim=dim(Y))

nsims=50000


tab %>% 
  mutate(D = as.numeric(Deaths)) %>% 
  group_by(Gender, Race) %>% 
  summarise(mean = log(mean(D/Population, na.rm = T))) %>% 
  ungroup() %>% 
  select(mean) %>% unlist() -> crude_data_beta

set.seed(1234)
beta0=sig2=tau2=array(dim=c(S,R,nsims))
lami=z=theta=array(dim=c(I,S,R,nsims))
beta0[,,1]=	 0 #############
Ymiss=list()

if(T == F){
  
  s = 1; r= 1
  s = 2; r = 1
  
}

  
system.time({
  count = 1
  
  for(s in 1:S){
    Ymiss[[s]]=list()
    for(r in 1:R){
      Ymiss[[s]][[r]] = array(dim=c(nsupp[s,r],nsims))
      Ymiss[[s]][[r]][,1] =	crude_data_beta[count] #############
      Y[!dY[,s,r],s,r]=Ymiss[[s]][[r]][,1]
      
      count = count + 1
    }
  }
  
  sig2[,,1]=	10 #############
  tau2[,,1]=	2 #############
  z[,,,1]=	0 #############
  for(i in 1:I){
    theta[i,,,1]=beta0[,,1] + z[i,,,1]
    lami[i,,,1]=exp(theta[i,,,1])
  }
  
  for(it in 2:nsims){
    for(s in 1:S){
      for(r in 1:R){
          ############
          #Update Y; account for state total and bounds
          if(nsupp[s,r]>0){
            nlam=n[!dY[,s,r],s,r]*lami[!dY[,s,r],s,r,it-1]
            good=FALSE; attempt=0
            while(!good){
            remain=nYmiss[s,r]
      
      if(attempt%%2==0){
            Yord=order(nlam)
            for(eye in 1:(nsupp[s,r]-1)){
              pie=nlam[Yord[eye]]/sum(nlam[Yord[eye:nsupp[s,r]]])
              probs=dbinom(Ythres[1]:Ythres[2],remain,pie)
              Ymiss[[s]][[r]][Yord[eye],it]=ifelse(remain==0,0,sample(Ythres[1]:Ythres[2],1,prob=probs))
              remain=remain-Ymiss[[s]][[r]][Yord[eye],it]
            }
            Ymiss[[s]][[r]][Yord[nsupp[s,r]],it]=remain
      }else{
            Yord=1:nsupp[s,r]
            for(eye in 1:(nsupp[s,r]-1)){
              pie=nlam[eye]/sum(nlam[eye:nsupp[s,r]])
              probs=dbinom(Ythres[1]:Ythres[2],remain,pie)
              Ymiss[[s]][[r]][eye,it]=ifelse(remain==0,0,sample(Ythres[1]:Ythres[2],1,prob=probs))
              remain=remain-Ymiss[[s]][[r]][eye,it]
            }
            Ymiss[[s]][[r]][nsupp[s,r],it]=remain
      }
            good=min(Ymiss[[s]][[r]][,it])>=Ythres[1] & 
                 max(Ymiss[[s]][[r]][,it])<=Ythres[2]
            attempt=attempt+1
            if(!good & attempt>=10){cat('attempt',attempt,'\n')}
            }
            Y[!dY[,s,r],s,r]=Ymiss[[s]][[r]][,it]
          }
      
          ############
          #update beta0
          bvar=(I/tau2[s,r,it-1] + 1/gam2)^(-1)
          bmean=bvar * (sum(theta[,s,r,it-1]-z[,s,r,it-1])/tau2[s,r,it-1])
          beta0[s,r,it]=rnorm(1,bmean,sqrt(bvar))
      
          ############
          #update z
          z[,s,r,it]=z[,s,r,it-1]
          for(i in 1:I){
            mui=mean(z[neigh[[i]],s,r,it])
            sigi=sig2[s,r,it-1]/m[i]
      
            zvar=1/(1/tau2[s,r,it-1] + 1/sigi)
            zmean=zvar * ( (theta[i,s,r,it-1]-beta0[s,r,it])/tau2[s,r,it-1] + mui/sigi )
            z[i,s,r,it]=rnorm(1,zmean,sqrt(zvar))
          }
          z[,s,r,it]=z[,s,r,it]-mean(z[,s,r,it])  #sum-to-zero constraint
      
          ############
          #update theta
          for(i in 1:I){
            ts=rnorm(1,theta[i,s,r,it-1],qt[i,s,r])
            ra=Y[i,s,r] * (ts - theta[i,s,r,it-1])
            rb=n[i,s,r] * (exp(ts) - exp(theta[i,s,r,it-1]))
            rc=((ts-beta0[s,r,it]-z[i,s,r,it])^2 - 
                (theta[i,s,r,it-1]-beta0[s,r,it]-z[i,s,r,it])^2)
            r0=exp(ra - rb - rc/(2*tau2[s,r,it-1]))
            theta[i,s,r,it]=ifelse(r0>runif(1),ts,theta[i,s,r,it-1])
      
            lami[i,s,r,it]=exp(theta[i,s,r,it])
          }
      
          ############
          #update tau2
          tau2[s,r,it]=1/rgamma(1,I/2 + at, 
                       sum((theta[,s,r,it]-beta0[s,r,it]-z[,s,r,it])^2)/2 + bt)
      
          ############
          #update sig2
          bss=0
          for(i in 1:I){
            bss=bss+z[i,s,r,it]^2*m[i] - z[i,s,r,it]*sum(z[neigh[[i]],s,r,it])
          }
          sig2[s,r,it]=1/rgamma(1,(I-1)/2 + as, bss/2 + bs)
        }
      }
    if(it/100 == floor(it/100)){
      print(it)
      cat('Iteration:',it,'\n')
        for(s in 1:S){
        for(r in 1:R){
        for(i in 1:I){
          trate=mean(theta[i,s,r,it-100+1:99]!=theta[i,s,r,it-100+2:100])
          trate=ifelse(trate>.75,.75,ifelse(trate<.2,.2,trate))
          qt[i,s,r]=qt[i,s,r]*trate/0.44
        }
        }
        }
    }
  }
  
  
  })

#    user  system elapsed 
#    951.690  31.459 992.511

beepr::beep(3)

##################
# export resutls 

final_resutls <- list(
  beta0,
  lami,
  z,
  theta,
  sig2,
  tau2 
)

names(final_resutls) <-c(
  'beta0',
  'lami',
  'z',
  'theta',
  'sig2',
  'tau2' 
)
  
write_rds(final_resutls, "./Final/res.rds")
