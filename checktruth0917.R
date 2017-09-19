#R function to compute true causal risk ratio by each time t under data generating models of the simulation section of Young and Tchetgen Tchetgen (2014)
#t is for follow-up time (in our example t=1,2...,5). Other parameters are as in 2014 paper.

causalRRt<-function(t,theta0,theta1,theta2,theta3,beta1){
	h0barsum1<-exp(theta0+theta1)/(1+exp(theta0+theta1))
	h0barsum2<-exp(theta0)/(1+exp(theta0))
	h0bar<-.5*(h0barsum1+h0barsum2)
	S0bar<-(1-h0bar)^(t)
	R0bar<-1-S0bar
	
	h1barsum1<-(exp(theta0+theta1+theta2+theta3)/(1+exp(theta0+theta1+theta2+theta3)))*(exp(beta1)/(1+exp(beta1)))
	h1barsum2<-(exp(theta0+theta2+theta3)/(1+exp(theta0+theta2+theta3)))*(1/(1+exp(beta1)))
	h1bar<-h1barsum1+h1barsum2
	S1bar<-(1-h1bar)^(t)
	R1bar<-1-S1bar
	# epsi0numsum1<-exp(theta1+theta2)/(1+exp(theta0+theta1+theta2))
	# epsi0numsum2<-exp(theta2)/(1+exp(theta0+theta2))
	# epsi0num<-epsi0numsum1+epsi0numsum2
	# epsi0densum1<-exp(theta1)/(1+exp(theta0+theta1))
	# epsi0densum2<-1/(1+exp(theta0))
	# epsi0den<-epsi0densum1+epsi0densum2
	# psi0<-log(epsi0num/epsi0den)
	# 
	# 
	# epsi1outnum<-exp(theta3)/(1+exp(beta1))
	# epsi1numsum1<-exp(theta1+beta1)/(1+exp(theta0+theta1+theta3))
	# epsi1numsum2<-1/(1+exp(theta0+theta3))
	# epsi1num<-epsi1outnum*(epsi1numsum1+epsi1numsum2)
	# epsi1densum1<-exp(theta1)/(1+exp(theta0+theta1))
	# epsi1densum2<-1/(1+exp(theta0))
	# epsi1den<-.5*(epsi1densum1+epsi1densum2)
	# psi1<-log(epsi1num/epsi1den)
	# 
	# 
	# epsi2numprod1sum1<-exp(theta1+beta1)/(1+exp(theta0+theta1+theta3))
	# epsi2numprod1sum2<-1/(1+exp(theta0+theta2+theta3))
	# epsi2numprod1<-epsi2numprod1sum1+epsi2numprod1sum2
	# epsi2numprod2sum1<-exp(theta1)/(1+exp(theta0+theta1))
	# epsi2numprod2sum2<-1/(1+exp(theta0))
	# epsi2numprod2<-epsi2numprod2sum1+epsi2numprod2sum2
	# epsi2num<-epsi2numprod1*epsi2numprod2
	# epsi2denprod1sum1<-exp(theta1)/(1+exp(theta0+theta1+theta2))
	# epsi2denprod1sum2<-1/(1+exp(theta0+theta2))
	# epsi2denprod1<-epsi2denprod1sum1+epsi2denprod1sum2
	# epsi2denprod2sum1<-exp(theta1+beta1)/(1+exp(theta0+theta1+theta3))
	# epsi2denprod2sum2<-1/(1+exp(theta0+theta3))	
	# epsi2denprod2<-epsi2denprod2sum1+epsi2denprod2sum2
	# epsi2den<-epsi2denprod1*epsi2denprod2
	# psi2<-log(epsi2num/epsi2den)
	# 
	# h1bar<-exp(psi0+psi1+psi2)*h0bar
	# S1bar<-(1-h1bar)^(t)
	# R1bar<-1-S1bar
	
	RR<-R1bar/R0bar
	
	#return(c(RR,psi0,psi1,psi2))
	return(c(RR,R1bar,R0bar))
}

nulltheta3<-function(theta1,beta1){
testnum<-.5*(exp(theta1)+1)
testden1<-exp(theta1+beta1)/(1+exp(beta1))
testden2<-1/(1+exp(beta1))
testden<-testden1+testden2
choice<-log(testnum/testden)
return(choice)
}

causalRRt(1,-7,-2,-.8,0,-2)
causalRRt(2,-7,-2,-.8,0,-2)
causalRRt(3,-7,-2,-.8,0,-2)
causalRRt(4,-7,-2,-.8,0,-2)
causalRRt(5,-7,-2,-.8,0,-2)

causalRRt(1,-5,-2,-.8,0,-2)
causalRRt(2,-5,-2,-.8,0,-2)
causalRRt(3,-5,-2,-.8,0,-2)
causalRRt(4,-5,-2,-.8,0,-2)
causalRRt(5,-5,-2,-.8,0,-2)
causalRRt(6,-5,-2,-.8,0,-2)
causalRRt(7,-5,-2,-.8,0,-2)

causalRRt(1,-5,-2,0,nulltheta3(-2,-2),-2)
causalRRt(2,-5,-2,0,nulltheta3(-2,-2),-2)
causalRRt(3,-5,-2,0,nulltheta3(-2,-2),-2)
causalRRt(4,-7,-2,0,nulltheta3(-2,-2),-2)
causalRRt(5,-7,-2,0,nulltheta3(-2,-2),-2)