# Load the libraries

library(lattice)
library(Rcpp)
library(grid)
library(ggplot2)
library(penalized)
library(parallel)
library(latticeExtra)
library(snow)
library(doSNOW)

# Get the data

source("https://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
library("DNAcopy")

# Few plots of the data

data(coriell)
y.all = coriell$Coriell.05296 # y.all is the whole signal
plot(y.all,xlab="Genome Index", ylab="Signal")

y = y.all[900:1500] 
plot(y, xlab="Genome Index", ylab="Signal", main="Osberved signal for the studied genomic region")

# EXERCISE 1 

# a)

# Forward alpha recursion and beta backward recursion

alpha_recursion = function(y, mu, A, m, s)
{
  K = length(mu)
  T = length(y)
  alpha = matrix(0, nrow=T,ncol=K)
  for (j in 1:K) alpha[1,j] = dnorm(y[1],m[j],s[j]) *sum(A[,j]* mu)
  for (t in 2:T) for (j in 1:K) alpha[t,j] = dnorm(y[t],m[j],s[j])*sum(A[,j]* alpha[t-1,])
  return(alpha)
}


beta_recursion = function(y, mu, A, m, s)
{
  K = length(mu)
  T = length(y)
  beta = matrix(0, nrow=T,ncol=K)
  for (j in 1:K) beta[T,j] = 1
  for (t in T:2) for (i in 1:K) beta[t-1,i] = sum(c(dnorm(y[t],m[1],s[1]),dnorm(y[t],m[2],s[2]),dnorm(y[t],m[3],s[3]))*A[i,]* beta[t,])
  return(beta)
}

mu <- c(1/3,1/3,1/3)
A <- matrix(0.005,nrow=3,ncol=3)+diag(3)*(0.99-0.005)
m <- c(0,0.5,-0.5)
s <- c(0.05,0.08,0.08)
y <- y[!is.na(y)]

# b)

alpha <- alpha_recursion(y=y, mu=mu, A=A, m=m, s=s)
beta <- beta_recursion(y=y, mu=mu, A=A, m=m, s=s)
normalisation <- 0
numerator <- matrix(0, nrow = length(y), ncol = length(mu))
marginal_posteriori <- 0


for (t in 1:length(y))
{
  normalisation[t] <- sum(alpha[t,]*beta[t,])
}


# We find that all the normalisation values have the same value for the computer.
# So we do not divide by the normalisation.

for (t in 1:length(y))
{
  for (k in 1:length(mu))
  {
   numerator[t,k] <- (alpha[t,k]*beta[t,k]) 
  }
  marginal_posteriori[t] <- which.max(numerator[t,])
}

# Plot the marginal a posteriori states at all the observations

Nstate_1 <- length(which(marginal_posteriori==1))
Nstate_2 <- length(which(marginal_posteriori==2))
Nstate_3 <- length(which(marginal_posteriori==3))

plot(marginal_posteriori, ylab='Most likely hidden state {1,2,3}', main='Marginal a posteriori position with the Normal
    emission distributions')

# EXERCISE 2 

# PART 1

# Define the new vector of emission densities  with the t distribution for state 1

pdf <- function(tt,j,y)
{
  if(j==1) {return((gamma(5/2)/(sqrt(pi)*0.1))*(1+(y[tt]^2)/(0.01))^(-5/2))}
  if(j==2) {return(dnorm(y[tt],0.5,0.08))}
  if(j==3) {return(dnorm(y[tt],-0.5,0.08))}
}

# alpha recursion

alpha_recursion = function(y, mu, A)
{
  K = length(mu)
  T = length(y)
  alpha = matrix(0, nrow=T,ncol=K)
  for (j in 1:K) alpha[1,j] = pdf(tt=1,j=j,y=y) *sum(A[,j]* mu)
  for (t in 2:T) for (j in 1:K) alpha[t,j] = pdf(tt=t,j=j,y=y)*sum(A[,j]* alpha[t-1,])
  return(alpha)
}

# beta recursion
beta_recursion = function(y, mu, A)
{
  K = length(mu)
  T = length(y)
  beta = matrix(0, nrow=T,ncol=K)
  for (j in 1:K) beta[T,j] = 1
  for (t in T:2) for (i in 1:K) beta[t-1,i] = sum(c((gamma(5/2)/(sqrt(pi)*0.1))*(1+(y[t]^2)/(0.01))^(-5/2),dnorm(y[t],m[2],s[2])
                                                    ,dnorm(y[t],m[3],s[3]))*A[i,]* beta[t,])
  return(beta)
}


mu <- c(1/3,1/3,1/3)
A <- matrix(0.005,nrow=3,ncol=3)+diag(3)*(0.99-0.005)
m <- c(0,0.5,-0.5)
s <- c(0.05,0.08,0.08)
y <- y[!is.na(y)]

alpha <- alpha_recursion(y=y, mu=mu, A=A)
beta <- beta_recursion(y=y, mu=mu, A=A)
normalisation <- 0
numerator <- matrix(0, nrow = length(y), ncol = length(mu))
marginal_posteriori <- 0

for (t in 1:length(y))
{
  normalisation[t] <- sum(alpha[t,]*beta[t,])
}

# We find that all the normalisation values have the same value for the computer.
# So we do not divide by the normalisation.

for (t in 1:length(y))
{
  for (k in 1:length(mu))
  {
    numerator[t,k] <- (alpha[t,k]*beta[t,k]) 
  }
  marginal_posteriori[t] <- which.max(numerator[t,])
}

# Plot the marginal maximum a posteriori states

plot(marginal_posteriori, ylab='Most likely hidden state {1,2,3}', main='Marginal a posteriori position with 
     t emission pdf for state 1')

# Get the counts of the number of times the hidden states occur

NNstate_1 <- length(which(marginal_posteriori==1))
NNstate_2 <- length(which(marginal_posteriori==2))
NNstate_3 <- length(which(marginal_posteriori==3))

# Compare the counts for the 2 models

c(Nstate_1,Nstate_2,Nstate_3)
c(NNstate_1,NNstate_2,NNstate_3)


# PART 2

# Compute the probabilities to be used in the EM algorithm

proba_xt_xt_1 <- function(y, alpha, beta, mu, A){
  K = length(mu)
  T = length(y)
  proba_xt_xt_1 = matrix(0, nrow = 3*K, ncol = T)
  for (t in 2:T){
    if(is.na(pdf(tt=t,j=1,y=y))){
      for (j in 1:K)  proba_xt_xt_1[j, t] = alpha[t-1,1]*A[1,j]*beta[t,j]
      for (j in (K+1):(2*K)) proba_xt_xt_1[j, t] = alpha[t-1,2]*A[2,j-K]*beta[t,j-K]
      for (j in (2*K+1):(3*K)) proba_xt_xt_1[j, t] = alpha[t-1,3]*A[3,j-2*K]*beta[t,j-2*K]
    }
    else{
      for (j in 1:K)  proba_xt_xt_1[j, t] = alpha[t-1,1]*pdf(tt=t,j=j,y=y)*A[1,j]*beta[t,j]
      for (j in (K+1):(2*K)) proba_xt_xt_1[j, t] = alpha[t-1,2]*pdf(tt=t,j=(j-K),y=y)*A[2,j-K]*beta[t,j-K]
      for (j in (2*K+1):(3*K)) proba_xt_xt_1[j, t] = alpha[t-1,3]*pdf(tt=t,j=(j-2*K),y=y)*A[3,j-2*K]*beta[t,j-2*K]
    }
  }
  return(proba_xt_xt_1)
}

# EM algorithm

EM_algorithm <- function(y, mu, epsilon, A_start){
  K = length(mu)
  T = length(y)
  A_hat = A_start
  alpha = alpha_recursion(y, mu, A_hat)
  beta = beta_recursion(y, mu, A_hat)
  proba_xt_xt_1 = proba_xt_xt_1(y, alpha, beta, mu, A_hat)
  N = rowSums(proba_xt_xt_1)
  A_hat = matrix(0, nrow = 3, ncol = 3)
  for (j in 1:K) for (i in 1:K) A_hat[i, j] = N[i+3*(j-1)]/(N[i] + N[i+3] + N[i+6])
  likelihood = sum(alpha[T,])
  print(likelihood)
  print(A_hat)
  n=1
  
  repeat{
    alpha = alpha_recursion(y, mu, A_hat)
    beta = beta_recursion(y, mu, A_hat)
    proba_xt_xt_1 = proba_xt_xt_1(y, alpha, beta, mu, A_hat)
    N = rowSums(proba_xt_xt_1)
    A_hat = matrix(0, nrow = 3, ncol = 3)
    for (j in 1:K) for (i in 1:K) A_hat[i, j] = N[i+3*(j-1)]/(N[i] + N[i+3] + N[i+6])
    new.likelihood = sum(alpha[T,])
    print(new.likelihood)
    print(A_hat)
    
    if(abs(new.likelihood - likelihood) < likelihood*epsilon){
      return(A_hat)
      break
    }
    else{
      likelihood = new.likelihood
    }
  }
}

# pick a starting matrix and use the EM to estimate the transition matrix (make sure that it is a valid
# transition matrix)

AA <- matrix(0.5,nrow=3,ncol=3)

EM_algorithm(y=y, mu=mu, epsilon=0.001, A_start=AA)
