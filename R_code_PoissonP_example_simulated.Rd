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

sim_data = function(theta,T)
{
N=rpois(1,theta[1]*T)
u=sort(runif(N,0,T))
p=g_func(theta[2:5],u)
ind=rbinom(N,1,p)
N=u[which(ind==1)]
(N)
}



T=100 # final time


g_func = function(theta,s)
{
	a1=theta[1]
	b1=theta[2]
	a2=theta[3]
	b2=theta[4]
	
	((1.1+sin(a1*sin(a2*s+b2)*s+b1))/2.11)
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
		f=sum(log(1-((1.1+sin(a1*sin(a2*s+b2)*s+b1)))/2.11))

		MCM=MCM+f
	}

mMCM=MCM/M

sgN=sum(log(g_func(x,N)))

(-sgN-mMCM)

}



theta=c(4,0.045,pi/2,0.09,0) # initial values of the parameters (lambda,a1,b1,a2,b2)


M=5000  # number of Monte Carlo samples for the E-step
MM=10 # number of iterations

R=30 # number of replications

MLE_dist=matrix(0,R,5)

set.seed(10)

for(j in 1:R)
{

theta=c(4,0.045,pi/2,0.09,0)

## generate data

N=sim_data(theta,T)

## Run MCEM algorithm


for(n in 1:MM)
{
	
	MCE=numeric(M)
	
	n_cores = detectCores()
	cluster=makeCluster(n_cores-1)
	registerDoParallel(cluster)

	MCE=foreach(i = 1:M) %dopar%
	{
		m=rpois(1,theta[1]*T)
		s=sort(runif(m,0,T))
		
		fg=g_func(theta[2:5],s)
		ind=rbinom(m,1,1-fg)
		s=s[which(ind==1)]
		Nx=length(s)
		

		MCE[i]=c(Nx,s)

	}
stopCluster(cl=cluster)


theta_hat=optim(theta[2:5],Q_func,gr=NULL,MCE,method="L-BFGS-B",lower=c(0.001,0,0.001,0),upper=c(Inf,2*pi,Inf,2*pi))


print(theta_hat$convergence)

theta_old=theta
theta=c(theta[1],theta_hat$par)

dist=sqrt(sum((theta_old-theta)^2))

print(theta)
print(dist)
print("")	
}

MLE_dist[j,]=theta

print(j)
print("")

}



curve(theta[1]*g_func(thetaR[2:5],x),0,100,ylim=c(0,4),xlab="time",ylab="Intensity function",cex.lab=1.5,cex.axis=1.5,lwd=2)

for( i in 1:50)
{
curve(MLE_dist[i,1]*g_func(MLE_dist[i,2:5],x),0,100,ylim=c(0,4),xlab="time",ylab="Intensity function",cex.lab=1.5,cex.axis=1.5,add=TRUE,lwd=1,col="grey")
}
curve(theta[1]*g_func(thetaR[2:5],x),0,100,lwd=2,add=TRUE)


