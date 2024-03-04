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

