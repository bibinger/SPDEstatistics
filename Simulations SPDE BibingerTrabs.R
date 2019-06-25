#1. Estimators from Section 3 of paper
#set number of Monte Carlo iterations:
N<-1
#Realizations of estimator (14) from paper
volaest<-numeric(N)
#Realizations of estimator (18) from paper
kappaest<-numeric(N)
#Realizations of estimator (16) from paper
quartest<-numeric(N)
#set number of equidistant spatial observations: (from 1/m to 1, note that right boundary here is not bounded away from 1 such that only first (m-1) statistics will be used)
m<-10
#array with realized volatilities:
RVarray<-array(0,dim=c(N,m))
#spatial domain:
xmin<-0
xmax<-1
#set constant curvature parameter: (model formulation as in (3) from paper)
kappa<-.4
#eigenfunctions (4) from paper:
eigenf<-function(j,x){sqrt(2/(xmax-xmin))*sin(j*x*pi/(xmax-xmin))*exp(-x/kappa)}
#eigenvalues (5) from paper:
lambda<-function(j){(2*kappa)^{-1}*(1+j^2*pi^2*kappa^2/(xmax-xmin)^2)}
#set number of equidistant observations over time:
n<-1000
#set constant volatility parameter:
sigma<-.25
#spectral cut-off frequency: (10000 is our standard choice, 100000 more precise but computation slow)
#J<-100000
J<-10000
#model parameters translated to (1) of paper: theta_0=0;theta_1=1;theta_2=kappa/2

for(k in 1:N){
#first simulate Ornstein-Uhlenbeck coordinate processes with exact simulation
OU<-array(0,dim=c(J,n+1))
OU[,1]<-rep(0,J)
for(j in 1:J){
for(i in 1:n){
OU[j,i+1]<-OU[j,i]*exp(-lambda(j)/n)+sigma*sqrt((1-exp(-2*lambda(j)/n))/(2*lambda(j)))*rnorm(1)}
}
#now simulate solution of SPDE via factor model as in (6) from paper
X<-array(0,dim=c(m,n+1))
for(i in 1:(n+1)){
for(l in 1:m){
X[l,i]<-sum(OU[1:J,i]*eigenf(1:J,l/m*(xmax-xmin)))
}}
#optional second-order term correction (factors close to one)
kor<-1-0.5*(2*sqrt(1+2*(1:(n)))-sqrt(2*(1:(n)))-sqrt(2+2*(1:(n))))
#kor<-rep(1,length=n)
RV<-numeric(m)
for(l in 1:m){RV[l]<-1/sqrt(n)*sum(1/kor*(X[l,2:(n+1)]-X[l,1:n])^2)}
RVarray[k,]<-RV
volaest[k]<-sqrt(pi*kappa/2)*mean(RV[1:(m-1)]*exp((1:(m-1))/m*2/kappa))
#this estimates the curvature parameter from Section 4 in the paper, which equals actually 2/kappa with kappa set above
kappaest[k]<-(log(RV[floor(m/4)])-log(RV[floor(m*3/4)]))/(floor(m*3/4)/m-floor(m*1/4)/m)
#quarticity estimator
RQ<-numeric(m-1)
for(l in 1:(m-1)){RQ[l]<-sum((X[l,2:(n+1)]-X[l,1:n])^4)}
quartest[k]<-1/3*pi*kappa/2*mean(RQ*exp((1:(m-1))/m*4/kappa))
}


#plot(colMeans(RVarray),xlab="",ylab="")
#lines(1/sqrt(pi*kappa/2)*sigma^2*exp(-2/kappa*(1:(m))/m),col=2,lwd=2)


#2. Estimators from Section 4 of paper
#set number of Monte Carlo iterations:
N<-1
#simulate model as in 1. above:
kappaest<-numeric(N)
varkappahat<-numeric(N)
IV0hat<-numeric(N)
m<-10
xmin<-0
xmax<-1
kappa<-.4
eigenf<-function(j,x){sqrt(2/(xmax-xmin))*sin(j*x*pi/(xmax-xmin))*exp(-x/kappa)}
lambda<-function(j){(2*kappa)^{-1}*(1+j^2*pi^2*kappa^2/(xmax-xmin)^2)}
n<-1000
sigma<-.25 #then IV0=sigma^2*sqrt(2/kappa)
J<-10000 #100000

for(k in 1:N){
OU<-array(0,dim=c(J,n+1))
OU[,1]<-rep(0,J)
for(j in 1:J){
for(i in 1:n){
OU[j,i+1]<-OU[j,i]*exp(-lambda(j)/n)+sigma*sqrt((1-exp(-2*lambda(j)/n))/(2*lambda(j)))*rnorm(1)}
}
X<-array(0,dim=c(m,n+1))
for(i in 1:(n+1)){
for(l in 1:m){
X[l,i]<-sum(OU[1:J,i]*eigenf(1:J,l/m*(xmax-xmin)))
}}
#use (m-1) spatial observations without y_m=1
RV<-numeric(m-1)
for(l in 1:(m-1)){RV[l]<-1/sqrt(n)*sum((X[l,2:(n+1)]-X[l,1:n])^2)}
#least squares estimator
Data<-data.frame(t(rbind(RV,(1:(m-1))/m)))
ls<-nls(formula=RV~theta1/sqrt(pi)*exp(-theta2*V2),start=list(theta1 = 1, theta2 = 1),data=Data)
IV0hat[k]<-coef(ls)[[1]]
varkappahat[k]<-coef(ls)[[2]]
kappaest[k]<-2/coef(ls)[[2]]
}


#for the computation of the asymptotic variance:
y<-(1:(m-1))/m
U<-matrix(c(mean(exp(-20*y)),-sigma^2*sqrt(2/kappa)*mean(y*exp(-20*y)),-sigma^2*sqrt(2/kappa)*mean(y*exp(-20*y)),sigma^4*2/kappa*mean(y^2*exp(-20*y))), nrow = 2, ncol = 2)
V<-matrix(c(mean(exp(-10*y)),-sigma^2*sqrt(2/kappa)*mean(y*exp(-10*y)),-sigma^2*sqrt(2/kappa)*mean(y*exp(-10*y)),sigma^4*2/kappa*mean(y^2*exp(-10*y))), nrow = 2, ncol = 2)
COV<-solve(V)%*%U%*%solve(V)
COV*0.75*pi*sigma^4*2/kappa
