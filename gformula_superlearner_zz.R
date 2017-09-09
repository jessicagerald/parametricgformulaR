##fit non-parametric model using Super Learner R package
#try out glm,CART,Random forest
#SuperLearner is not allowing missing value, need to set them to be 0s first
stest <- test1
stest[is.na(stest)] <- 0
#SL.library<-c("SL.glm","SL.ipredbagg","SL.randomForest","SL.glmnet")   #this will use all the observations and no screening
SL.library<-c("SL.rpart")  
fitL.2<-SuperLearner(stest$L,X=as.data.frame(cbind(stest$avglag.L.csum,stest$lag.A.csum,stest$X,stest$t0)),
                     family="binomial",SL.library=SL.library,verbose=F)

# fit <- rpart(L~avglag.L.csum+lag.A.csum+X+t0, data =stest)
# summary(fit)
# new<-data.frame(cbind(stest$avglag.L.csum,stest$lag.A.csum,stest$X,stest$t0))
# colnames(new)<-c("avglag.L.csum","lag.A.csum","X","t0")
# preds <- predict(fit, newdata=new[1,], type="vector")

set.seed(12345) 
sl1<-data.frame()
for (jj in 1:N){
  subj<-subset(stest,id==jj)
  X <- rep(subj$X[1],Ndat)
  Astar<-rep(1,Ndat)
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
      colnames(d1)<-c("V1","V2","V3","V4")
      L[kk]<-rbinom(1,1,predict(fitL.2,type="vector",newdata=data.frame(d1))$pred)
      #L<-rbinom(1,1,predict(fitL.2,type="response",newdata=as.data.frame(d1)))
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
  sl1<- rbind(sl1,cbind(id=jj,t,t0,X, Astar,L,lag.A,lag.L,L.csum,lag.L.csum,lag.cumAstar,Py,prodp0,prodp1,poprisk))
}
meanpoprisk_sl<-function(time){mean(sl1$poprisk[sl1$t0==time])} 
sapply(0:(Ndat-1),meanpoprisk_sl)
