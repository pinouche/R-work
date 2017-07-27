larva<- read.table(file="C:/Users/pinouche/Downloads/beetlelarva.txt",header=TRUE)
fix(larva)



library(lattice)
library(nlme)
library(MASS)
install.packages("Rcpp")
install.packages("ggplot2")
library(Rcpp)
library(grid)
library(ggplot2)
require(reshape2)
require(ggplot2)
install.packages('sna')
library('sna')
install.packages('coda')
library(coda)

# explanatory analysis with a heatmap of the data

Easting <- larva$easting
Northing <- larva$northing

ggplot(larva, aes(Easting, Northing)) +
  theme(plot.title = element_text(face="bold", 
                                  size=20, hjust=0.5)) + 
  theme(axis.title = element_text(face="bold", size=16,
                                  hjust=0.5)) +
  geom_tile(aes(fill = count)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(text = element_text(size=16)) 

# Prior elicitation

# Prior number 1 

param1 <- rnorm(10000, mean = 0, sd = 1)
param2 <- rnorm(10000, mean = 0, sd = 1)

expvector1 <- 0
for (i in 1:10000)
{
  expvector1[i] <- exp(param1[i]+param2[i])
}

# Priot number 2

param1 <- rnorm(10000, mean = 0, sd = sqrt(5))
param2 <- rnorm(10000, mean = 0, sd = sqrt(5))

expvector <- 0
for (i in 1:10000)
{
  expvector[i] <- exp(param1[i]+param2[i])
}

# Plot histograms of the simulated data

par(mfrow = c(1,2))
hist(expvector1,breaks=200,col='red',xlim=c(0,40),ylim=c(0,8000),xlab = 'Number of larvae',main='Simulated number of larvae for 10000 simulations
      and using N(0,1) as priors')
hist(expvector,breaks=10000,col=2,xlim=c(0,500),ylim=c(0,8000),xlab = 'Number of larvae',main='Simulated number of larvae for 10000 simulations
      and using N(0,5) as priors')

# Create the X matrix 

to_dummy = function(X) {
  out = data.frame(matrix(nrow=nrow(X), ncol=0))
  for (val in unique(X$northing)) {
    out[paste("northing", val, sep="_")] = ifelse(X$northing==val, 1, 0)
  }
  
  for (val in unique(X$easting)) {
    out[paste("easting", val, sep="_")] = ifelse(X$easting==val, 1, 0)
  }
  
  return(data.matrix(out))
}
mat <- to_dummy(larva)

# likelihood 


loglik <- function(larva,beta)
{
  loglik <- 0
  mu <- 0
  
  for (i in 1:48)
  {
    
    mu[i] <- exp(mat[i,]%*%beta)
    loglik <- loglik + log(dpois(larva$count[i],mu[i]))
    
  }
  return(loglik)
}

# First Prior

lpr<-function(beta) {
  sum(log(dnorm(beta)))
} 


# Second Prior

lpr<-function(beta) {
  sum(log(dnorm(beta,sqrt(10))))
} 

#initialise (could use glm fit)

beta0=c(rep(0,times=14))

#MCMC loop - here "beta" is the current state of the Markov chain.
#betap will be the proposed state

MCMC<-function(K=100000,beta=beta0) {
  #returns K samples from posterior using MCMC
  #no subsampling, start state goes in beta
  
  B=matrix(NA,K,14); LP=rep(NA,K); LLK=rep(NA,K)
  #storage, I will write the sampled betas here
  
  lp=loglik(larva,beta)+lpr(beta) 
  #log posterior is log likelihood + log prior + constant
  
  for (k in 1:K) {
    
    #tuned RW MH - I adjusted the step sizes so they were 
    #unequal for beta[1] and beta[2] 
    betap=rnorm(14,beta,0.1)#generate candidate
    
    LLK[k]=loglik(larva,betap)
    lpp=LLK[k]+lpr(betap)        #calculate log post for candidate
    
    MHR=lpp-lp                         #"log(pi(y)/pi(x))" 
    print(MHR)
    if (log(runif(1))<MHR) 
      {           #Metropolis Hastings acceptance step
      beta=betap                       #if we accept update the state
      lp=lpp
    }
    
    B[k,]=beta                         #save the sequence of MCMC states, our samples.
    LP[k]=lp
  }
  return(list(B=B,L=LP,LL=LLK))
} 

# Stored values for 2 runs with prior 1
K <- 100000
Output1 <- MCMC(K,beta=beta0);
beta1<- c(rep(1,times=14))
beta2 <- c(rep(10,times=14))
Output2 <- MCMC(K,beta=beta1);
Output3 <- MCMC(K,beta=beta2);

Run1Prior1B <- Output1$B # Result of the second MCMC run with different starting values
Run1Prior1LP <- Output1$L
Run1Prior1LLLK1 <- Output1$LL
Run2Prior1B <- Output2$B # Result of the second MCMC run with different starting values
Run2Prior1LP <- Output2$L
Run2Prior1LLLK1 <- Output2$LL
Run3Prior1B <- Output3$B # Result of the second MCMC run with different starting values
Run3Prior1LP <- Output3$L
Run3Prior1LLLK1 <- Output3$LL

# Plot the densities for each run on top of each other

par(mfrow=c(1,3))
plot(density(Run1Prior1LP), xlim=c(-158,-130),col='red',ylim=c(0,0.25),xlab='Value of log posterior',
     main = 'Densities of log posteriors for 3 MCMC runs
    with standard normal priors')
lines(density(Run2Prior1LP),col='blue')
lines(density(Run3Prior1LP),col='green')
legend( 'topright', inset=0, 
        legend=c(expression(bold("Densities")),"First run","Second run",'Third run'),
        col=c(NA,'red','blue','green'),
        lty=c(NA,1,1,1), merge=FALSE, cex=0.7)

acf(Run1Prior1B[,1],lag.max=3000,ylab='Autocorrelation',main='Autocorrelation plot for the Beta_Northing6')
acf(Run1Prior1LP,lag.max=3000,ylab='Autocorrelation',main='Autocorrelation plot for the log posterior')


# convergence plot for Beta_Northing6

plot(Run1Prior1B[,1],type='l',xlab='Samples',ylab='Parameter Value',main='Time series of Beta_Northing6 
     for 100000 samples')
abline(h=mean(Run1Prior1B[,1]),col='red')
legend( 'topright', inset=0, 
        legend=c('Posterior mean'),
        col=c('red'),
        lty=c(1), merge=FALSE, cex=0.7)

hist(Run1Prior1B[,1],breaks=50,xlab='Parameter value',ylab='Frequency',main='Histogram of Beta_Northing6 
     for 100000 samples')
abline(v=mean(Run1Prior1B[,1]),col='red')
legend( 'topright', inset=0, 
        legend=c('Posterior mean'),
        col=c('red'),
        lty=c(1), merge=FALSE, cex=0.7)
# Stored values for 1 run with the second prior N(0,10)

Run1Prior2B <- Output1$B # Result of the second MCMC run with different starting values
Run1Prior2LP <- Output1$L
Run1Prior2LLLK1 <- Output1$LL

# Plot the 2 densities for the 2 runs to see if they agree (first prior)

# Computing the harmonic estimates for the two priors and then compute bayes factor
p_hat <- 1/(mean(1/(exp(Run1Prior1LLLK1))))
p_hat2 <- 1/(mean(1/(exp(Run1Prior2LLLK1))))   

# Bayes factor
Bayesfac <- (p_hat/p_hat2)

# Get the HPD 95% confidence itnerval

HPDinterval(as.mcmc(Run1Prior1B))

# Obtain the posterior mean for each parameter

PosteriorMeanVec <- 0
for (i in 1:14)
{
  
 PosteriorMeanVec[i] <- mean(BB[,i])
 
}

PosteriorMeanVec

# Obtain the predicted number of larvas on each of the 48 plots, using the posterior means as point estimates
  
mu2 <- 0

for (i in 1:48)
{
    
  mu2[i] <- exp(mat[i,]%*%PosteriorMeanVec)
    
}






