###R program for g-formula in observational study
#Created on May 24th, 2017
#directory of the GitHub folder: C:\Users\zzhang\Documents\GitHub\parametricgformulaR
rm(list=ls())
library(data.table)
library(splines)
library(FSA)
library(zoo)
library(boot)
library(plyr)
library(SuperLearner)

########################                                                                            
####1.Read in dataset###
########################
test1<-read.table(file="C:/Users/zzhang/Desktop/Gformula Causal Inference for Longitudinal Data/test1.csv",sep=",")

##########################                                                                            
####2. Running Gformula###
##########################

#####Step I: Model fitting
##NOTE: If we read in external dataset, need to find out the number of time points and number of subjects first.
pred.func.L<-function(model,dat){
  fitL<-glm(model, family = binomial, data=dat)
  return(fitL)
}
# fit the parametric model using average cumulative information of L and A, baseline non-vary baseline 
# covariate X and spline function of time points 
#fitL.1<-pred.func.L(L ~ avglag.L.csum + lag.A.csum + X + ns(t0,2), subset(test1,t0>0)) ##Limit to time>0
#fitL.1<-glm(L ~ avglag.L.csum + lag.A.csum + X + ns(t0,2) , family = binomial, data=test1)

pred.func.Y<-function(model,dat){
  fitY<-glm(model, family = binomial, data=dat)
  return(fitY)
}
##fit parametric outcome model using pooled records (added the current observations A and L):
#fitY.1<-pred.func.Y(Y ~ A + L + lag.A + lag.L + X + t0,subset(test1,t0>0))
#fitY.1<-glm(Y ~ A + L + lag.A + lag.L + X + t0, family = binomial, data=test1) 


#####Step II: Monte Carlo simulation for each subject
##Static treatment: always treat vs. always not treat
startTime <-proc.time()
set.seed(9375)
N <- 100         # number of subjects
Ndat <-6        # number of time points
tp <-Ndat-1
no.sample <- 100 # Number of subjects, for each bootstrapping. Same as the number of subjects in the real data
Bsample <- 20    # Number of bootstrapping, make this number small first
Regimes <- c(0:1)# Never treat-->0 or/and always treat-->1
Result_AT <-matrix(NA, nrow=Bsample, ncol=Ndat)
Result_NT <-matrix(NA, nrow=Bsample, ncol=Ndat)
colnames(Result_AT) <-paste("TimePoint",0:tp,sep="")
colnames(Result_NT) <-paste("TimePoint",0:tp,sep="")
RiskRatio <-list(AT=Result_AT,NT=Result_NT)

#Function to simulate using pooled over time information.
gformula_pool_sim <- function(trt){
  
  for (jj in 1:N){
    subj<-subset(sample,id==jj)
    X <- rep(subj$X[1],Ndat)
    Astar<-rep(trt,Ndat)
    lag.A<-as.vector(stats::lag(zoo(Astar),-1,na.pad=TRUE))
    lag.cumAstar<-pcumsum(Astar)
    
    t0<-seq(0,Ndat-1,by=1)
    t <-seq(1,Ndat,by=1)
    
    for (kk in 1:Ndat){
      if (kk==1){
        L<-lag.L<-Y<-Py<-prodp0<-prodp1<-prodp2<-lag.L.csum<-L.csum<-avglag.L.csum<-NULL
        L[kk]<-subj$L[kk]
        lag.L[kk] <- Y[kk]<-lag.A[kk]<-lag.L.csum[kk]<-avglag.L.csum[kk]<-L.csum[kk]<-0
        L.csum[kk]<-L[kk]
        d0<-cbind(Astar[kk],L[kk],lag.A[kk],lag.L[kk],X[kk],t0[kk])
        colnames(d0)<-c("A","L","lag.A","lag.L","X","t0")
        Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
        prodp0[kk]<- 1-Py[kk]
        prodp1[kk]<- Py[kk]
      }
      else if (kk>1){
        lag.L.csum[kk]<-L.csum[kk-1]
        avglag.L.csum[kk]<-lag.L.csum[kk]/t0[kk]
        d1<-cbind(avglag.L.csum[kk],lag.cumAstar[kk],X[kk],t0[kk])
        colnames(d1)<-c("avglag.L.csum","lag.A.csum","X","t0")
        L[kk]<-rbinom(1,1,predict(fitL.1,type="response",newdata=data.frame(d1)))
        L.csum[kk]<-cumsum(L)[kk]
        lag.L[kk]<-L[kk-1]
        d0<-cbind(Astar[kk],L[kk],lag.A[kk],lag.L[kk],X[kk],t0[kk])
        colnames(d0)<-c("A","L","lag.A","lag.L","X","t0")
        Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
        
        prodp0[kk]<- 1-Py[kk]
        prodp1[kk]<- Py[kk]*cumprod(prodp0[1:kk-1])[kk-1]
      }
      poprisk <- cumsum(prodp1)
    }
    pool<- rbind(pool,cbind(id=jj,t,t0,X, Astar,L,lag.A,lag.L,L.csum,lag.L.csum,lag.cumAstar,Py,prodp0,prodp1,poprisk))
    }
  return(pool)
}
#Function to calculate population risk for each time points.
meanpoprisk_pool<-function(time){mean(pool$poprisk[pool$t0==time])} 

#Bootstrapping Result
for (i in 1:Bsample){
  print("...............")
  print(paste("Simulation:",i,"..............."))
  
  ids<-data.frame(sample(test1$id,no.sample,replace = TRUE))
  ids$newid<-1:no.sample
  colnames(ids)<-c("id","newid")
  test1a <- as.data.table(test1)
  setkey(test1a, "id")
  sample <- test1a[J(ids), allow.cartesian = TRUE]  # create the new data set names "sample"
  sample$id<-NULL
  sample$id<-sample$newid
  
  fitL.1<-pred.func.L(L ~ avglag.L.csum + lag.A.csum + X + ns(t0,2), subset(sample,t0>0))
  fitY.1<-pred.func.Y(Y ~ A + L + lag.A + lag.L + X + t0,subset(sample,t0>0))
  for (g in Regimes)
  if (g == 0){
  pool<-data.frame()
  pool<-gformula_pool_sim(trt=0)
  Result_NT[i,]<-sapply(0:(Ndat-1),meanpoprisk_pool)
  } else if (g ==1){
    pool<-data.frame()
    pool<-gformula_pool_sim(trt=1)
    Result_AT[i,]<-sapply(0:(Ndat-1),meanpoprisk_pool)
  }
 RiskRatio$AT<-Result_AT
 RiskRatio$NT<-Result_NT
}
elapsedtime<-proc.time()-startTime 
elapsedtime 

#RiskRatio for always treat and never treat
RR<-RiskRatio$AT/RiskRatio$NT
RR




