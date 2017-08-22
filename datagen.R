####################################
####1. Generating a simple dataset##
####################################
rm(list=ls())
library(data.table)
library(splines)
library(FSA)
library(zoo)
library(boot)
library(plyr)
# Generate Data -  binary outcome Y, binary treatment indicator A
set.seed(1234) 
# number of subjects
N <- 50
# number of time points
Ndat <-6
# number of replicates
test1<-data.frame()
for(s in 1:N){
  for (j in 1:Ndat){
    if(j == 1){
      L<-Y<-A<-C<-t<-t0<-NULL
      t[j]<-j
      t0[j]<-j-1
      #Baseline covariate, does not vary within each subject, e.g. gender
      X <- rbinom(1,1,0.5)  
      #X <- rnorm(1, 0, 0.2)
      #Baseline L0, binary
      L[j] <- rbinom(1, 1, 0.5)
      #Baseline exposure level, binary
      A[j] <- rbinom(1, 1, plogis(0.5))
      #Baseline outcome at t0, binary, assumed to be all 0 at this time
      Y[j] <- 0
      #Censoring at t0
      C[j] <- 0
    }
    else if(j >= 2){
      t[j]<-j
      t0[j]<-j-1
      X[j]<-X[j-1]
      L[j]<- rbinom(1, 1, plogis(2*A[j-1]+0.3*X[j]))
      A[j]<- rbinom(1, 1, plogis(0.5+0.03*L[j]+0.3*X[j]))
      Y[j]<- rbinom(1, 1, plogis(-2-0.8*A[j]-2*L[j-1]+1.5*A[j-1]+0.3*X[j]))
      C[j]<- rbinom(1, 1, plogis(-0.5-0.5*L[j-1]-0.1*A[j-1]+X[j]))
    }
    {if (C[j]==1){Y[j]<-NA}}
    {if (Y[j]==1|C[j]==1){break}}
  }
  L.csum<- cumsum(L)
  A.csum<- cumsum(A)
  lag.A<-as.vector(stats::lag(zoo(A),-1,na.pad=TRUE))
  lag.L<-as.vector(stats::lag(zoo(L),-1,na.pad=TRUE))
  lag.A.csum<-as.vector(stats::lag(zoo(A.csum),-1,na.pad=TRUE))
  lag.L.csum<-as.vector(stats::lag(zoo(L.csum),-1,na.pad=TRUE))
  avglag.L.csum<-lag.L.csum/t0
  test1 <- rbind(test1,cbind(id=s,t,t0,X,L,A,Y,C,lag.A,lag.L,A.csum,L.csum,lag.A.csum,lag.L.csum,avglag.L.csum))
}

#Write the fake dataset into csv file
write.table(test1, 
                   file="C:/Users/zzhang/Desktop/Gformula Causal Inference for Longitudinal Data/test1.csv",
                   sep=",")
