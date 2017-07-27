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
library(doParallel)
install.packages("doSNOW")
install.packages("doParallel")


# Get the data to test the algorithm

data("nki70")
nki70
dim(nki70)
markers.df <- as.data.frame(nki70[,-c(2:7)])
dim(markers.df)

# First step. Function no parallelised to obtain the initial predictive correlation. 

xval = function(data, k, lambda1, lambda2) {
  
  n = nrow(data)
  
  folds = split(sample(n), seq_len(k))
  
  xval.fold = function(fold) {
    
    dtrain = data[-fold, ]
    dtest = data[fold, ]
    
    ElasticNet.fold = penalized(time, penalized = ~ . - time, data = dtrain, lambda1 = lambda1, lambda2 = lambda2, trace = FALSE)
    
    pred = predict(ElasticNet.fold, data = dtest)
    
    return(data.frame(PRED = pred[, "mu"], OBS = dtest$time))
    
  } 
  
  pairs = lapply(folds, xval.fold)
  return(do.call("rbind", pairs))
} 

pairs = xval(markers.df,k=10, lambda1 = 10, lambda2 = 10)
cor(pairs)[1,2]

# Make 10 runs of the 10 fold cross validation function and take the average

mean(replicate(10,cor(xval(markers.df, k=10, lambda1 = 10, lambda2 = 10))[1,2]))

# Parallel function to implement the algorithm

parallel.xval = function(data, cluster, k, lambda1, lambda2) 
{
  
  n = nrow(data)
  
  folds = split(sample(n), seq_len(k))
  
  xval.fold = function(fold, lambda1, lambda2) {
    
    dtrain = data[-fold, ]
    dtest = data[fold, ]
    
    ElasticNet.fold = penalized(time, penalized = ~ . - time, data = dtrain, lambda1 = lambda1, lambda2 = lambda2, trace = FALSE)
    
    pred = predict(ElasticNet.fold, data = dtest)
    
    return(data.frame(PRED = pred[, "mu"], OBS = dtest$time))
  } 
  
  # Added code to implement parallelisation
  
  clusterExport(cluster, list("data"))
  clusterEvalQ(cluster, library(penalized))
  
  pairs = parLapply(cluster, folds, xval.fold, lambda1 = lambda1, lambda2 = lambda2)
  return(do.call("rbind", pairs))
  
}
cl = makeCluster(2)
pairs = parallel.xval(markers.df, cl, k=10, lambda1 = 10, lambda2 = 10)
stopCluster(cl)
cor(pairs)

mean(replicate(10,cor(parallel.xval(markers.df,cl, k=10, lambda1 = 10, lambda2 = 10))[1,2]))


# Find the right parameters by now implementing the grid-search algorithm

tune.elasticNet = function(data, cluster, k, initial.lambda1, initial.lambda2, stepping1, stepping2, epsilon) 
  {
  
  cur.lambda1 = initial.lambda1
  cur.lambda2 = initial.lambda2
  
  predcorVec <- 0
  
  #registerDoParallel(cluster,list('parralel.xval'))
  #clusterExport(cluster, list("pararlel.xval"))
  #foreach(i=1:10) %dopar% {predcorVec[i] <- cor(parallel.xval(data, cluster = cluster, k = k, lambda1 = cur.lambda1, lambda2 = cur.lambda2))[1,2]}
  
  #cur.predcor <- mean(predcorVec)
  
  cur.predcor <- mean(replicate(10,cor(parallel.xval(data, cluster = cluster, k = k, lambda1 = cur.lambda1, lambda2 = cur.lambda2))[1,2]))
  
  # repeat the body until the condition in if is not met for the 9 combinations of lambdas at each iteration
  new.lambda1 <- 0
  new.lambda2 <- 0

  repeat {
    
    new.lambda1 <- cur.lambda1 + seq(-stepping1, stepping1, stepping1)
    new.lambda2 <- cur.lambda2 + seq(-stepping2, stepping2, stepping2)
    correl <- matrix(0, nrow=3, ncol=3)
    
    for (j in 1:3)
    {
      for (m in 1:3)
      {
      correl[j,m] <- mean(replicate(10,cor(parallel.xval(data, cluster=cluster, k=k, lambda1 = new.lambda1[j], lambda2 = new.lambda2[m]))[1,2]))
      }
    }
  
      if (max(correl) - cur.predcor >= epsilon) 
    {
      
      lambda1_ind <- which(correl == max(correl), arr.ind = TRUE)[1]
      lambda2_ind <- which(correl == max(correl), arr.ind = TRUE)[2]
      cur.lambda1 <- new.lambda1[lambda1_ind]
      cur.lambda2 <- new.lambda2[lambda2_ind]
      cur.predcor <- max(correl) 
      
    } #THEN
    
    
    else {
    
      break
    } #ELSE
    
    if(cur.lambda1 < stepping1) {stepping1 <- cur.lambda1/2}
    
    if(cur.lambda2 < stepping2) {stepping2 <- cur.lambda2/2}  
    
  } #REPEAT
  
  return(c(cur.predcor = cur.predcor, cur.lambda1 = cur.lambda1, cur.lambda2 = cur.lambda2))
  
} #TUNE.ELASTICNET

cl = makeCluster(2)
tune.elasticNet(markers.df, cluster = cl, k = 10, initial.lambda1 = 5, initial.lambda2 = 5, stepping1=6, stepping2=6, epsilon=1e-4)
stopCluster(cl)

system.time(tune.elasticNet(markers.df, cluster = cl, k = 10, initial.lambda1 = 25, initial.lambda2 = 25, stepping1=6, stepping2=6, epsilon=1e-4))


# Plot the regression with the initial parameters and the new ones

# PLot the above regression 

xyplot(nki70[,1] ~ fitted(ElasticNet2), main = "Elastic-Net Regression",
       xlab = "Fitted values", ylab = "Observed values",
       panel = function(...) {
         panel.xyplot(...)
         panel.lmline(..., col = "green", lty = 2)
       })

# Initial 

ElasticNet = penalized(time, penalized = markers.df[,2:71] , data = nki70, lambda1 = 20, lambda2 = 20, trace = FALSE)

# After cross validation 

ElasticNet2 = penalized(time, penalized = markers.df[,2:71] , data = nki70, lambda1 = 3.96, lambda2 = 2.4, trace = FALSE)


cl <- makeCluster(2)
registerDoParallel(cl)
foreach(i=1:3) %dopar% sqrt(i)

