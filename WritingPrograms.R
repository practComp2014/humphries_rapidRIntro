#===========================================
## @knitr code_part0
#  TITLE:  computational economics: assignment 1 solutions
#  AUTHOR: John Eric Humphries
#  (this code draws from a pset I assigned in econ 21410 as well as oliver Browne's solutions)
#  abstract: an applied introductino to R.
#
#  Date: April-8-14
#================================


#========================
# Section 0: setup
#========================

#setwd("/mnt/ide0/home/johneric/sbox/projects/neighborhoodVision/")
rm(list=ls())           # Clear the workspace
set.seed(907) 
library(ggplot2)
library(mnormt)
library(sandwich)
library(car)
library(xtable)
library(aod)
library(AER)
library(MASS)

#====================
# Section 1: generating data
#====================
GenerateDataSimple = function(n){
  # Generates data with from a simple linear model
  # args: n - number of observations to generate
  # output: A list containing a Y vector and an X matrix
  x1     <- rbinom(n,10,.5)
  x2     <- rnorm(n,20,10)
  X      <- cbind(1,x1,x2)
  theta2 <- 5
  eps    <- rnorm(n,0,sqrt(4))
  beta   <- matrix(cbind(2,3,1.4),3,1)
  Y   <-  X %*% beta + eps
  colnames(X) <- c("const", "x1","x2")
  colnames(Y) <- c("Y")
  return(list(Y,X))
}

# Simulating simple data
data <- data.frame(GenerateDataSimple(100))
attach(data)
X = cbind(const,x1,x2)

#=========================
# Section 2: Defining sub functions
#=========================

SimpleOLS <- function(y=y,x=x) {
  beta_hat1 = solve(t(x) %*% x) %*% (t(x)%*% y)
  se1       = sqrt(diag((1/(length(y) - length(beta_hat1))) * (t(y - x%*%beta_hat1) %*% (y-x%*%beta_hat1))[1] * solve(t(x) %*% x)))
  out = list(t(beta_hat1),se1)
  return(out)
}



SumOfSquares = function(b,x=X,y=Y)
{
  #Sum of squares function for numerical optimizer
  sum((y - x%*%b)^2)
}

Llik = function(pars,x=X,y=Y)
{
  #Log Likelihood function for normal OLS with numerical optimizer
  b = pars[1:3]
  t2 = pars[4]
  n  = length(y)
  if (t2 <=.000001) value = 1000000000000000000000
  else  value = -1 * (-n / 2 * log(2*pi) - n/2* log(t2) - 1/(2* t2)* (t(y - x%*%b) %*% (y - x%*%b)))
  value
}

OLS.gradient <- function(b,x=X,y=Y){
  #Returns the analytic gradient of OLS
  return(-2*t(X)%*%(Y - X%*%b))
}

# Bayesian Regression
bayesRegress <- function(x=X,y=Y,ndraws=5000){
  # A function that calculates bayesian linear regression with "flat" priors built in.
  #posteriors
  pn <- solve(t(x)%*%x)
  mn <- (pn)%*%(t(x)%*%y)
  an <- 1/2 * length(y) 
  bn <- 1/2 * ( t(y)%*%(diag(rep(1,length(y)))  -  x%*%pn%*%t(x))%*%y )
  draws = matrix(NA,ndraws,dim(x)[2]+1)
  draws[,4]   <- 1/rgamma(ndraws,an,bn)
  draws[,1:3] <- t(sapply(c(1:ndraws), function(i) rmnorm(1, mn ,draws[i,4]*pn)))
  beta <- apply(draws,2,mean)
  se   <- apply(draws,2,sd)
  return(list(beta,se))
}

lik_bootstrap = function(pars = c(1,1,1,1000), x=X,y=Y,bsnum=bsnum) {   
  bs_estimates = matrix(NA,bsnum, (dim(x)[2] +1) )
  n = length(y)
  for (i in 1:bsnum) {
    samp  = sample(n,n,replace=T)
    Xboot = x[samp,]
    Yboot = y[samp] 
    bs_estimates[i,] <- optim(par=pars, Llik,  x=Xboot, y=Yboot )$par
  }   
  standard_errors = apply(bs_estimates,2,sd)
  return (standard_errors)
}


#=================================
#
#=================================

OLS = function(x=X,y=Y,method=1,optim.method ="Nelder-Mead",gradient=NULL)
{
  #Calculates OLS using one of four methods:
  # Method == 0 uses the standard lm function
  # Method == 1 Calculates OLS algebraically
  # Method == 2 Uses an optimizer to minimize the sum of squares
  # Method == 3 Uses an optimizer to maximize a likelihood
  # Method == 4 Run Bayesian Regression.
  if(method==0){
    result = lm(y ~ x -1 )
    beta_hat0 = result$coefficients
    se0       = summary(result)$coefficients[,2]
    out = list(beta_hat0,se0)
  }else if(method==1){
    beta_hat1 = solve(t(x) %*% x) %*% (t(x)%*% y)
    se1       = sqrt(diag((1/(length(y) - length(beta_hat1))) * (t(y - x%*%beta_hat1) %*% (y-x%*%beta_hat1))[1] * solve(t(x) %*% x)))
    out = list(t(beta_hat1),se1)
  }else if(method==2){
    beta_hat2 <- optim( c(0,0,0), SumOfSquares, method = "BFGS", x=x, y=y, hessian=T, gr=OLS.gradient)
    out = list(beta_hat2$par,sqrt(diag(solve(beta_hat2$hessian))))
  }else if(method==3){
    beta_hat3 = optim( c(1,1,1,1000), Llik, method = optim.method,  x=x, y=y)$par
    se3 = lik_bootstrap(pars = c(1,1,1,1000), x=X,y=Y,bsnum=200)
    out = list(beta_hat3[1:3],se3[1:3])
  }else if(method==4){
    beta_hat4 = bayesRegress(x=x, y=y)
    out = list(beta_hat4[[1]][1:3],beta_hat4[[2]][1:3])    
  }else{
    stop("Error! Method not found")
  }
  class(out) <- "olsoutput"
  return(out)
}

#============================
# Formatting output of a function (classes)
#============================

#
a = c(1:1000)
b = rep(c("a", "b", "c"),20)
class(a)
class(b)
summary(a)
summary(b)

#
b = as.factor(b)
class(b)
class(a)

# Formatting output for summary command on my new object class. 
summary.olsoutput <- function(object) {
  beta  <- t(as.matrix(object[[1]]))
  se    <- as.matrix(object[[2]])
  out   <- cbind(beta,se)
  colnames(out) <- c("beta","se")
  out   <- round(out, 3)
  print("Regression Results")
  print(paste("",round(beta,3),""))
  print(paste("(",round(se,3),")", sep= ""))
  #return(out)
}

