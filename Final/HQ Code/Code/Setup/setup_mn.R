#####################
#Set Things Up
#####################
source('./Final/HQ Code/Code/Setup/setup_packages.r')
#library(spdep)      # spatial dependencies functions
#library(rgdal)      # manages spatial projections
#library(maptools)   # general functions for map/spatial manipulation
#library(sp)         # foundational definition of spatial classes in R
#library(arm)        # Loads R2WinBugs, Lattice, Matrix, others
#library(raster)
#library(rgeos)

load(file='Shapefiles/counties.rdata')
load(file='Shapefiles/states.rdata')

us0=spTransform(us,CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 
     +x_0=0.0 +y_0=0 +k=1.0 +units=m +no_defs"))

mn=us0[us0$STATE_FIPS=='27',]

## Creating adjacency matrix
mn.wt1<-poly2nb(mn,snap=1000)
#snap=1000 is used because *for this shapefile* some neighbor pairs
#  were not being recognized; this value of "snap" seems to fix that
#By default, poly2nb uses a "queen" definition of contiguity
#  (i.e., sharing a single *point* is sufficient to be a neighbor);
#  this can be changed by using the argument "queen=FALSE"
#  (which should require neighbors to share an *edge*)

mn.WBwt<-nb2WB(mn.wt1)
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

