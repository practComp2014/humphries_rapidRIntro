#==========================================================================
#
#  Practical Computing 2014-04-25
#
#  R implementation of a factor Model with MCMC. 
#  Author: John Eric Humphries
#
#===========================================================================



#================================
# Section 0: Setup
#================================

# set up
#setwd("C:/Users/Acer/Dropbox/JohnEricHumphries_EconGradSchool/CourseMaterial/BayesianMetrics")
#setwd("/mnt/ide0/home/johneric/sbox/courses/BayesianMetrics/finalProject")
setwd("~/sbox/projects/factorXcorr")
rm(list=ls())           # Clear the workspace
set.seed(907)           # set random seed

# Libraries
library(mnormt)
#library(compiler)
#library(doMC)
#enableJIT(3) # precompiling R Code


#================================
# Section 1: Generating Data
#================================

# Model Structure 
n            <- 500          # Number of observations
K            <- 1            # number of factors
J            <- 12           # number of measures

#Building Xs and fs
Var         <- matrix(c(1.2, 0,
                         0,1.1),2,2) #
Mean        <- c(1.5,3)
dta         <- rmnorm(n, mean=Mean,varcov =Var )
x           <- dta[,1]   # First covariate
f.true      <- dta[,2] 
alpha.true  <- matrix(c(   1   , 2.14,   .22,  .93,  .253,   .17, 1.5,   .85,3,2,1,4),J,1)                           # The factor loadings
beta.true   <- matrix(c(   .234, 1.11, 1.514,  4.83,2.596,  3.66, 3.542, 2.32,1,2,3,1.5),J,1)
X           <- x
sigma2.true <-        c(.5,   3, .8, .7,  1, .3, 0.8, 1,.5,1.2,.8,1.8)
#building Y
Y <- matrix(NA,n,J)
for (j in 1:J) {
  Y[,j] <- X*beta.true[j,] +   f.true*alpha.true[j,]  + rnorm(n, mean=0, sd=sqrt(sigma2.true[j]))   
}

#================================
# Section 2: The Gibbs Sampler
#================================

#--------------------------------
# Section 2.0: Defining Values
#-------------------------------
nfactors <- 1
measures <- J



#--------------------------------
# Section 2.1: Starting Values
#-------------------------------
alpha  <-  alpha.true * 1.2 * c(1,.1,3,.1,.5,2) # starting at true values to begin.  
f      <-  f.true * c(.5,1.4)
sigma2 <-  sigma2.true * 5
F1     <-  Var[1,1] 
beta   <-  beta.true * .7 *  c(.6,1.2,.2,1,3,.2)


#--------------------------------
# Section 2.3: The gibbs sampler (R VERSION)
#--------------------------------

factorRGibbs <- function(Niter,burnin,thin,alpha,f,sigma2,F1,beta,Y,X) { 
  
  #Alpha Update#
  alphaUpdate = function(Ya,fa,sigma2)   {
    for (miter in (1+nfactors):measures) {
      A1                 <- solve(  1/sigma2[miter] * (t(fa)%*%fa))
      alpha1             <- A1 %*% (  1/sigma2[miter] * t(fa) %*% (Ya[,miter]) )
      alpha[miter,] <- rmnorm(1,alpha1,A1)    
    }
    alpha
  }
  
  #Beta Update#
  betaUpdate = function(Yb,Xb,sigma2) {
    for (miter in (1):measures) {
      B1                 <- solve(  1/sigma2[miter] * (t(Xb)%*%Xb))  
      beta1             <- B1 %*% ( 1/sigma2[miter] * t(Xb) %*% (Yb[,miter]) )
      beta[miter,] <- rmnorm(1,beta1,B1)    
    }
    beta
  }
  
  #Factor Update
  factorUpdate <- function(Y,alpha,sigma2) { 
    F1  <- solve( t(alpha) %*%diag(1/sigma2)%*% (alpha))
    #     f <- factorStepRcpp(n,F1,alphaf,F0invf0,Yf,sigma2)
    #f       <- apply(Ya,1,factorUpdate)
    for (m in 1:n) { # n = number of agents
      #F1      <- solve(F0inv + t(alpha) %*%diag(1/sigma2)%*% (alpha))
      f1  <- F1 %*% ( t(alpha) %*%diag(1/sigma2)%*% (Yf[m,])  )
      f[m]       <- rnorm(1,f1,F1)
    }
    f
  }
  
  
    
  # The Loop ----------------------------------------
  draws  = matrix(NA,Niter,length(alpha) + length(sigma2) + length(f) +length(beta) )  
  for (iter in 1:(burnin + Niter*thin)){ 
  
    #Draw from Sigma (informal "flat" prior)
    a1      <- ( apply(Y,2,length))/2
    b1      <- rep(NA,measures)
    for (i in 1:measures) 
      b1[i] <- (  sum((Y[,i] - f*alpha[i,] - X*beta[i,])^2))/2
    sigma2  <-  1/rgamma(measures,a1,b1)
    
   # Draw from alpha
    Ya = Y - X%*%t(beta)
    alpha <- alphaUpdate(Ya,f,sigma2)

    # Draw from beta
    Yb    <- Y - f%*%t(alpha)
    beta  <- betaUpdate(Yb,X,sigma2)
    
    
    # Updating Factor
    Yf = Y - X%*%t(beta)
    f  <- factorUpdate(Yf,alpha,sigma2)
     

    # Saving output
    if( (iter - burnin)>0 & iter%%thin==0) draws[(iter-burnin)/thin,] <- c(alpha,sigma2,beta,f)
    if (0 == (iter %% 100)) print(paste("Iteration ", iter))
  }
  return(draws)
}

#Running the function
results <- factorRGibbs(1000,1000,1,alpha,f,sigma2,F1,beta,Y,X)


#===========================
# Section 3: Checking Results
#===========================

alpha.mean = matrix(apply(results,2,mean)[1:J],J,1)
sigma2.mean = apply(results,2,mean)[13:24]
beta.mean = matrix(apply(results,2,mean)[25:36],J,1)
f.mean = (apply(results,2,mean)[37:dim(results)[2]])
cbind(alpha.true,alpha.mean)
cbind(beta.true,beta.mean)
cor(c(f.mean),c(f.true))
cbind(sigma2.mean,sigma2.true)
c( cor(x,f.mean), cor(x,f.true))




#----
# In parallel
#-----
#registerDoMC(10)
#results <- foreach(i = 1:10, .combine=rbind) %dopar% {
#  factorRGibbs(200,2000,20,alpha,f,sigma2,F1,beta,Y,X)
#}
#
#results2 <- factorRGibbs(5000,10000,10,alpha,f,sigma2,F1,beta,Y,X)


#library(doParallel)
#cl <- makeCluster(4)
#registerDoParallel(cl)
#registerDoMC(10)
#results <- foreach(i = 1:10, .combine=rbind) %dopar% {
#  library(mnormt)
#  factorRGibbs(100,100,1,alpha,f,sigma2,F1,beta,Y,X)
#}
#

#stopCluster(cl)
