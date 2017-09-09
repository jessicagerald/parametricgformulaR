###R program for g-formula in observational study
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
#test1<-read.table(file=dir,sep=",")
test1<-test2014
test1[is.na(test1)] <- 0
##########################                                                                            
####2. Running Gformula###
##########################

#####Step I: Model fitting
##NOTE: If we read in external dataset, need to find out the number of time points and number of subjects first.
pred.func.L<-function(model,dat){
  fitL<-glm(model, family = binomial, subset(dat,t0>0))
  return(fitL)
}
# fit the parametric model using average cumulative information of L and A, baseline non-vary baseline 
# covariate X and spline function of time points 
pred.func.Y<-function(model,dat){
  fitY<-glm(model, family = binomial, data=dat)
  return(fitY)
}

#####Step II: Monte Carlo simulation for each subject
##Static treatment: always treat vs. always not treat
startTime <-proc.time()
set.seed(9375)
N <- 1000         # number of subjects
Ndat <-6        # number of time points
tp <-Ndat-1
no.sample <- 1000 # Number of subjects, for each bootstrapping. Same as the number of subjects in the real data
Bsample <- 1    # Number of bootstrapping, make this number small first
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
    #X <- rep(subj$X[1],Ndat)
    Astar<-rep(trt,Ndat)
    lag.A<-as.vector(stats::lag(zoo(Astar),-1,na.pad=TRUE))
    
    t0<-seq(0,Ndat-1,by=1)
    t <-seq(1,Ndat,by=1)
    
    for (kk in 1:Ndat){
      if (kk==1){
        L<-lag.L<-Y<-Py<-prodp0<-prodp1<-prodp2<-NULL
        L[kk]<-subj$L[kk]
        lag.A[kk] <- Y[kk]<-0
        d0<-cbind(L[kk],Astar[kk],lag.A[kk])
        colnames(d0)<-c("L","A","lag.A")
        Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
        prodp0[kk]<- 1-Py[kk]
        prodp1[kk]<- Py[kk]
      }
      else if (kk>1){
        d1<-cbind(lag.A[kk])
        colnames(d1)<-c("lag.A")
        L[kk]<-rbinom(1,1,predict(fitL.1,type="response",newdata=data.frame(d1)))
        lag.L[kk]<-L[kk-1]
        d0<-cbind(L[kk],Astar[kk],lag.A[kk])
        colnames(d0)<-c("L","A","lag.A")
        Py[kk]<-predict(fitY.1,type="response",newdata=data.frame(d0))
        prodp0[kk]<- 1-Py[kk]
        prodp1[kk]<- Py[kk]*cumprod(prodp0[1:kk-1])[kk-1]
      }
      poprisk <- cumsum(prodp1)
    }
    pool<- rbind(pool,cbind(id=jj,t,t0, Astar,L,lag.A,Py,prodp0,prodp1,poprisk))
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
  
  fitL.1<-pred.func.L(L ~ lag.A, sample)
  fitY.1<-pred.func.Y(Y ~ L + A + lag.A, sample)
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

apply(RR,2,mean)
apply(RR,2,sd)


