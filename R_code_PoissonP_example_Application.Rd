install.packages("stats")
install.packages("foreach")
install.packages("iterators")
install.packages("parallel")
install.packages("doParallel")
install.packages("stats")
install.packages("Kendall")
install.packages("doSNOW")
install.packages("mvnfast")
install.packages("Rfast")

library(stats)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(Kendall)
library(doSNOW)
library(mvnfast)
library(Rfast)

### Poisson process example

### Application - Coal mine

N=read.table("data_coal.txt",header=F)
N=t(t(N))


Ny=length(N)

T=112 # final time


g_func = function(theta,s)
{
	a1=theta[1]
	b1=theta[2]
	a2=theta[3]
	b2=theta[4]
	
	((0.3+3.2*(1-pnorm(s,a1,b1))+1.2*dnorm(s,a2,b2)/dnorm(a2,a2,b2))/4.7)
}


Q_func = function(x,MCE)
{
 
	a1=x[1]
	b1=x[2]
	a2=x[3]
	b2=x[4]
	
 	MCM=0
	
	for(i in 1:M)
	{

		MCM2=matrix(unlist(MCE[i]),ncol=1,byrow=T)
		L=MCM2[1]
		s=MCM2[2:(L+1)]
		f=sum(log(1-((0.3+3.2*(1-pnorm(s,a1,b1))+1.2*dnorm(s,a2,b2)/dnorm(a2,a2,b2))/4.7)))

		MCM=MCM+f
	}

mMCM=MCM/M

sgN=sum(log(g_func(x,N)))

(-sgN-mMCM)

}



theta=c(42,12,85,8) # initial values of the parameters (mu1,sigma1,mu2,sigma2)

lambda=4.7

M=5000  # number of Monte Carlo samples for the E-step
eps=10^(-5) # Stopping rule for the algorithm: |theta-theta'| < eps
MM=10 # maximum number of iterations
n=0 # iterations counter


## Run MCEM algorithm

dist=1

print(theta)
print("")

while(dist>=eps & n<MM)
{
	
	n=n+1;

	MCE=numeric(M)
	
	n_cores = detectCores()
	cluster=makeCluster(n_cores-1)
	registerDoParallel(cluster)

	MCE=foreach(i = 1:M) %dopar%
	{
		m=rpois(1,lambda*T)
		s=sort(runif(m,0,T))
		
		fg=g_func(theta[1:4],s)
		ind=rbinom(m,1,1-fg)
		s=s[which(ind==1)]
		Nx=length(s)
		

		MCE[i]=c(Nx,s)
	}
stopCluster(cl=cluster)



theta_hat=optim(theta[1:4],Q_func,gr=NULL,MCE,method="L-BFGS-B",lower=rep(0,4),upper=c(Inf,Inf,Inf,Inf))



print(theta_hat$convergence)

theta_old=theta
theta=theta_hat$par

dist=sqrt(sum((theta_old-theta)^2))

print(theta)
print(dist)
print("")	
}




### MLE=(41.819579, 16.543211, 86.306012 , 7.004204)


curve(lambda*g_func(theta[1:4],x),0,112,lwd=2,ylim=c(0,3.5),xlab="time",ylab="Intensity function",cex.lab=1.5,cex.axis=1.5)
points(N,rep(0,length(N)),pch=1,cex=0.5)


