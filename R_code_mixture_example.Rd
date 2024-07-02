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




#### Point-mass mixture example ####


n=1500    # sample size
M=200000   # number of Monte Carlo samples

p=0.6     # probability of zero rain
alpha=0.5 # shape parameter of the Gamma distribution
beta=0.05 # rate parameter of the Gamma distribution

f=function(x)
{
	x1=x[1]
	x2=x[2]
	-n1*x1*log(x2) + n1*lgamma(x1) - (x1-1)*Sly +x2*Sy
}

df=function(x)
{
	x1=x[1]
	x2=x[2]
	c( -n1*log(x2) + n1*digamma(x1) - Sly, -n1*x1/x2 + Sy )
}



### Simulation

Sys.time()

MLE=numeric(M)

	n_cores = detectCores()
	cluster=makeCluster(n_cores-1)
	registerDoParallel(cluster)

	
	MLE=foreach(i = 1:M) %dopar%
	{
		n0=rbinom(1,n,p)
		n1=n-n0

		y=rgamma(n1,shape=alpha,rate=beta)

		Sy=sum(y)
		Sly=sum(log(y))
		
		p_hat=n0/n
		
		beta_0=mean(y)/var(y)
		alpha_0=beta_0*mean(y)
		
		theta_hat=optim(c(alpha_0,beta_0),f,df,method="L-BFGS-B",lower=c(0.001,0.001),upper=c(Inf,Inf))
		theta_hat=c(p_hat,theta_hat$par)

		MLE[i]=theta_hat
		
	}

stopCluster(cl=cluster)

MLE=matrix(unlist(MLE),ncol=3,byrow=T)

Sys.time()

### Plots

par(mfrow=c(1,3))

plot(density(MLE[,1]),lwd=2,main="",xlab="p",ylab="",cex.lab=1.5,cex.axis=1.5)
abline(v=p,lwd=2,col="red")

#quartz()     #for Mac
#win.graph()  #for Windows

plot(density(MLE[,2]),lwd=2,main="",xlab=expression(alpha),ylab="",cex.lab=1.5,cex.axis=1.5)
abline(v=alpha,lwd=2,col="red")

#quartz()     #for Mac
#win.graph()  #for Windows

plot(density(MLE[,3]),lwd=2,main="",xlab=expression(beta),ylab="",cex.lab=1.5,cex.axis=1.5)
abline(v=beta,lwd=2,col="red")



### Application - Rainfall in Rio Grando do Sul State, Brazil


data=read.table("rain_data_00_23.txt",sep="\t",header=F)
data=t(t(data))
data=data[,2]

n=length(data)
n0=length(which(data==0))
n1=n-n0
y=data[which(data!=0)]

Sy=sum(y)
Sly=sum(log(y))
		
p_hat=n0/n
		
beta_0=mean(y)/var(y)
alpha_0=beta_0*mean(y)
		
theta_hat=optim(c(alpha_0,beta_0),f,df,method="L-BFGS-B",lower=c(0.001,0.001),upper=c(Inf,Inf))
theta_hat=c(p_hat,theta_hat$par)

# MLE=(0.58921792 0.51970596 0.05347131)

hist(y,prob=T,breaks=seq(0,150,5),main='',xlab='',ylab='')
curve(dgamma(x,theta_hat[2],theta_hat[3]),0,150,add=TRUE,lwd=2,col="red")

