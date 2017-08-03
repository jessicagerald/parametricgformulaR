###Simulation test code for g-formula in observational study
#May 24th, 2017
rm(list=ls())
library(data.table)
library(splines)
library(FSA)
library(zoo)
library(boot)
library(plyr)
library(dplyr)
####################################
####1. Generating a simple dataset
####################################
# Generate Data -  binary outcome Y, binary treatment indicator A
set.seed(1234) 
# number of subjects
N <- 20
# number of time points
Ndat <-5
# number of replicates
#e = list()
test1<-data.frame()
for(s in 1:N){
  #print("...............")
  #print(paste("Patient:",s,"..............."))
  #Normal random variable U with mean 0 and sd 0.2
  for (j in 1:Ndat){
    if(j == 1){
    L<-Y<-A<-C<-p<-p0<-NULL
    p[j]<-j
    p0[j]<-j-1
    #Baseline covariate, does not vary within each subject, e.g. gender
    X <- rbinom(1,1,0.5)  
    #U <- rnorm(1, 0, 0.2)
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
      p[j]<-j
      p0[j]<-j-1
      X[j]<-X[j-1]
      L[j]<- rbinom(1, 1, plogis(2*A[j-1]+0.3*X[j]))
      A[j]<- rbinom(1, 1, plogis(0.5+0.03*L[j]+0.3*X[j]))
      Y[j]<- rbinom(1, 1, plogis(-2-0.8*A[j]-2*L[j-1]+1.5*A[j-1]+0.3*X[j]))
      C[j]<- rbinom(1, 1, plogis(-0.5-0.5*L[j-1]-0.1*A[j-1]+X[j]))
    }
    {if (Y[j]==1){break}}
    
    }
  #e[[s]]<- data.frame(id=s,p=j,X=X,L=L,A=A,Y=Y)
  L.csum<- cumsum(L)
  A.csum<- cumsum(A)
  lag.A<-as.vector(lag(zoo(A),-1,na.pad=TRUE))
  lag.L<-as.vector(lag(zoo(L),-1,na.pad=TRUE))
  lag.A.csum<-as.vector(lag(zoo(A.csum),-1,na.pad=TRUE))
  lag.L.csum<-as.vector(lag(zoo(L.csum),-1,na.pad=TRUE))
  avglag.L.csum<-lag.L.csum/p0
  test1 <- rbind(test1,cbind(id=s,p,p0,X,L,A,Y,C,lag.A,lag.L,A.csum,L.csum,lag.A.csum,lag.L.csum,avglag.L.csum))
}
#test1[1:10,]

# #We use "boot" to get bootstrapping results
# Lpred <-function(formula,data,indices){
# d     <-data[indices,]
# fitL.1<-glm(formula, family = binomial, d)
# #return(fitL.1)
# return(fitL.1$coefficients[1])
# }
# 
# cptime<-proc.time()
# Lboot <- boot(data=test1, R=3, formula=L ~ avglag.L.csum + lag.A.csum + X + ns(p0,2))
# proc.time()
# 
# #We use "boot" to get bootstrapping results
# Lpred <-function(formula,data){
#   #d     <-data[indices,]
#   fitL.1<-glm(formula, family = binomial, data)
#   return(fitL.1)
#   #return(fitL.1$coefficients[1])
# }
# Lpred(L ~ avglag.L.csum + lag.A.csum + X + ns(p0,2),1)
# 


#manipulate the dataset
#test1<-do.call(rbind, lapply(e, data.frame, stringsAsFactors=FALSE))
#table(test1$id)
#test1$p<-ave(test1$id,test1$id,FUN=seq_along) #create index number by subject
#test1$p0<-test1$p-1

#Generate cumulative L
#test1$L.csum <- ave(test1$L,test1$id, FUN=cumsum)
#test1$A.csum <- ave(test1$A,test1$id, FUN=cumsum)
#test1<-data.table(test1)
# original<-c("A","L")
# lags <-paste("lag",original,sep=".")
# test1[, lag.A1 := shift(A,1), by=id]
# test1[, lag.A2 := shift(A,2), by=id]
# test1[, lag.A3 := shift(A,3), by=id]
# 
# test1[, lag.L1 := shift(L,1), by=id]
# test1[, lag.L2 := shift(L,2), by=id]
# test1[, lag.L3 := shift(L,3), by=id]

#test1[, lag.L.csum:= shift(L.csum,1), by=id]
#test1[, lag.A.csum:= shift(A.csum,1), by=id]
#test1[, lag.A:= shift(A,1), by=id]
#test1[, lag.L:= shift(L,1), by=id]
#test1$avglag.L.csum<-test1$lag.L.csum/test1$p0
#for (i in seq_along(test1)) set(test1, i=which(is.na(test1[[i]])), j=i, value=0)


########################                                                                            
####2. Analyses Steps###
########################
#####Step I:
#Ltp <- data.frame(alpha0=rep(NA,p0),alpha1=rep(NA,p0),alpha2=rep(NA,p0),alphaX=rep(NA,p0))
tp <- Ndat-1
for (ii in 1:tp){
  #fitL<-list()
  dat<-subset(test1,p0==ii)
  #Model for each pooled time point using previous L and A with lag 1
  fitL[[ii]]<-glm(L ~ lag.L + lag.A + X , family = binomial, data=dat)
}

# Model using average cumulative information of L and A, baseline non-vary X and spline function of time 
fitL.1<-glm(L ~ avglag.L.csum + lag.A.csum + X + ns(p0,2) , family = binomial, data=test1)


##fit parametric model for Y:
fitY.1<-glm(Y ~ lag.A.csum + lag.L + X + p0, family = binomial, data=test1) 

#####Step II: Monte Carlo simulation for each subject
##Static treatment: always treat vs. always not treat
#Static treatment: always treat 
set.seed(12345) 
new<-data.frame()
for (jj in 1:N){
  subj<-subset(test1,id==jj)
  X<-subj$X
  Astar<-rep(1,length(subj$p))
  lag.A<-as.vector(lag(zoo(Astar),-1,na.pad=TRUE))
  lag.cumAstar<-pcumsum(Astar)
  t0<-seq(0,length(subj$p)-1,by=1)
  t <-seq(1,length(subj$p),by=1)

  for (kk in 1:length(subj$p)){

    if (kk==1){
      L<-lag.L<-Y<-Py<-prodp0<-prodp1<-prodp2<-NULL
      L[kk]<-subj$L[kk]
      lag.L[kk]<- Y[kk]<-lag.A[kk]<-0
      d0<-cbind(lag.cumAstar[kk],lag.L[kk],X[kk],t0[kk])
      colnames(d0)<-c("lag.A.csum","lag.L","X","p0")
      Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
      prodp0[kk]<- Py[kk]
      prodp1[kk]<- Py[kk]
    }
    else if (kk>1){
      lag.L[kk]<-L[kk-1]
      d1<-cbind(lag.L[kk],lag.A[kk],X[kk])
      colnames(d1)<-c("lag.L","lag.A","X")
      L[kk]<-rbinom(1,1,predict(fitL[[kk-1]],type="response",newdata=data.frame(d1)))
      d0<-cbind(lag.cumAstar[kk],lag.L[kk],X[kk],t0[kk])
      colnames(d0)<-c("lag.A.csum","lag.L","X","p0")
      Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
      prodp0[kk]<- 1-Py[kk]
      prodp1[kk]<-cumprod(prodp0[1:kk])[kk]
    }
    #poprisk <- cumsum(prodp1)[length(subj$p)]
    poprisk <- cumsum(prodp1)
  }
  new<- rbind(new,cbind(id=jj,t,t0, Astar, X, L, lag.L,Py,prodp0,prodp1,poprisk))
}
#new[1:5,]
#Static treatment: always not treat
set.seed(12345) 
new0<-data.frame()
for (jj in 1:N){
  subj<-subset(test1,id==jj)
  X<-subj$X
  Astar<-rep(0,length(subj$p))
  lag.A<-as.vector(lag(zoo(Astar),-1,na.pad=TRUE))
  lag.cumAstar<-pcumsum(Astar)
  t0<-seq(0,length(subj$p)-1,by=1)
  t <-seq(1,length(subj$p),by=1)
  
  for (kk in 1:length(subj$p)){
    
    if (kk==1){
      L<-lag.L<-Y<-Py<-prodp0<-prodp1<-prodp2<-NULL
      L[kk]<-subj$L[kk]
      lag.L[kk]<- Y[kk]<-lag.A[kk]<-0
      d0<-cbind(lag.cumAstar[kk],lag.L[kk],X[kk],t0[kk])
      colnames(d0)<-c("lag.A.csum","lag.L","X","p0")
      Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
      prodp0[kk]<- Py[kk]
      prodp1[kk]<- Py[kk]
    }
    else if (kk>1){
      lag.L[kk]<-L[kk-1]
      d1<-cbind(lag.L[kk],lag.A[kk],X[kk])
      colnames(d1)<-c("lag.L","lag.A","X")
      L[kk]<-rbinom(1,1,predict(fitL[[kk-1]],type="response",newdata=data.frame(d1)))
      d0<-cbind(lag.cumAstar[kk],lag.L[kk],X[kk],t0[kk])
      colnames(d0)<-c("lag.A.csum","lag.L","X","p0")
      Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
      prodp0[kk]<- 1-Py[kk]
      prodp1[kk]<-cumprod(prodp0[1:kk])[kk]
    }
    #poprisk <- cumsum(prodp1)[length(subj$p)]
    poprisk <- cumsum(prodp1)
  }
  new0<- rbind(new0,cbind(id=jj,t,t0, Astar, X, L, lag.L,Py,prodp0,prodp1,poprisk))
}
new0[1:5,]


##Using pooled model to simulate L and Y
#Always treat

set.seed(12345) 
pool<-data.frame()
for (jj in 1:N){
  subj<-subset(test1,id==jj)
  X<-subj$X
  Astar<-rep(1,length(subj$p))
  lag.A<-as.vector(lag(zoo(Astar),-1,na.pad=TRUE))
  lag.cumAstar<-pcumsum(Astar)

  t0<-seq(0,length(subj$p)-1,by=1)
  t <-seq(1,length(subj$p),by=1)
  
  for (kk in 1:length(subj$p)){
    
    if (kk==1){
      L<-lag.L<-Y<-Py<-prodp0<-prodp1<-prodp2<-lag.L.csum<-L.csum<-avglag.L.csum<-NULL
      L[kk]<-subj$L[kk]
      lag.L[kk]<- Y[kk]<-lag.A[kk]<-lag.L.csum[kk]<-avglag.L.csum[kk]<-L.csum[kk]<-0
      L.csum[kk]<-L[kk]
      
      d0<-cbind(lag.cumAstar[kk],lag.L[kk],X[kk],t0[kk])
      colnames(d0)<-c("lag.A.csum","lag.L","X","p0")
      Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
      prodp0[kk]<- Py[kk]
      prodp1[kk]<- Py[kk]
    }
    else if (kk>1){
      lag.L[kk]<-L[kk-1]
      #L.csum[kk]<-cumsum(L[1:kk])[kk]
      L.csum[kk]<-cumsum(L)
      lag.L.csum[kk]<-L.csum[kk-1]
      avglag.L.csum[kk]<-lag.L.csum[kk]/t0[kk]
      d1<-cbind(avglag.L.csum[kk],lag.cumAstar[kk],X[kk],t0[kk])
      colnames(d1)<-c("avglag.L.csum","lag.A.csum","X","p0")
      L[kk]<-rbinom(1,1,predict(fitL.1,type="response",newdata=data.frame(d1)))
      
      d0<-cbind(lag.cumAstar[kk],lag.L[kk],X[kk],t0[kk])
      colnames(d0)<-c("lag.A.csum","lag.L","X","p0")
      Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
      prodp0[kk]<- 1-Py[kk]
      prodp1[kk]<-cumprod(prodp0[1:kk])[kk]
    }
    poprisk <- cumsum(prodp1)
  }
 pool<- rbind(pool,cbind(id=jj,t,t0, Astar, X, L, lag.L,L.csum,lag.L.csum,lag.A.csum,Py,prodp0,prodp1,poprisk))
}

#III. Sampling

ids<-data.frame(sample(test1$id,20,replace = TRUE))
colnames(ids)<-c("id")
#sample1<- test1[test1$id %in% ids, ]
#sample1<-merge(id, test1, by="id",all.x = TRUE)
test1a <- as.data.table(test1)
setkey(test1a, "id")
# create the new data set
sample1 <- test1a[J(ids), allow.cartesian = TRUE]
