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