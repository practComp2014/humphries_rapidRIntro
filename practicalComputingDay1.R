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
library(systemfit)
library(MASS)


#=======================
## @knitr code_part1
# Section 1: intro to R.
#=======================

# Getting started
draws <- rnorm(25)
summary(draws)
class(draws)

#Making and manipulating matricies
draws.m <- matrix(draws,5,5)
class(draws.m)
draws.m[1,1]
draws.m[7]
draws.m[4,]
diag(draws.m)
a = c(2,3)
b = c(4,5)
d = c(7)
c(a,b)
cbind(a,b)
rbind(a,b)

c(a,d)
# broadcasting!
cbind(a,d)

#logical statements on matricies
draws[draws>0]
draws>0
which(draws>0)
any(draws>0)
all(draws>0)
all.equal(c(2,2),c(2,2))
as.character(draws)

# Getting to know lists
colors <- c("red", "blue", "yello")
mylist <- list(draws,draws.m,colors)
mylist
mylist[[1]]
mylist[[2]]
mylist[[3]]
mylist[[3]][1]

#============================
## @knitr code_part2
# Section 2: Create data generating functions and generate data
#=============================

#generate data and learn about control-flow.

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
data <- GenerateDataSimple(100)
X <- data[[2]]
Y <- data[[1]]

GenerateData = function(n=100,J=1,nX=5,het=T) {
    #This functino generates data from multiple linear equations
    #
    # args: n = sample size; J = number of equations , X = number of Xs
    # 
    # Notes: generates multiple linear equatinos, need to add exclusions.
    # Error handling:
    if(length(n)>1 | length(J)>1 | length(nX)>1 | class(n)!="numeric" | class(J)!="numeric" | class(nX)!="numeric")
        stop("wrong length or class of inputs!")
    covar        <- matrix(runif(nX^2,min=0, max =.9),nX,nX)
    diag(covar)  <- runif(nX,min=1,max=8) # this function is super overloaded, which makes it a pain, always check output.
    mean         <- runif(nX,min=-5,max=5) + rnorm(nX)
    x            <- cbind(1,matrix(rmnorm(n, mean=mean,varcov = covar ),n,nX))
    beta.true    <- cbind(1,matrix(runif(J*nX,-5,10),J,nX,byrow=T))
    if (het ==F)
        sigma2.true  <- matrix(runif(J,0,12),1,J)
    if (het ==T)
        sigma2.true  <- matrix(runif(n*J,6,18),n,J)
    # Generating outcomes
    y <-matrix(NA,n,J)
    for (iter in 1:J) {
        y[,iter]     <- x%*% as.matrix(beta.true[iter,]) +   rnorm(n, mean=0, sd=sqrt(sigma2.true[,iter])) 
    } # end of for loop
    return(list(x,y,beta.true,sigma2.true))
} # end of function call

# Simulating heterogeneous data
data.big  <- GenerateData(n=100,J=1,nX=2)
x <- data.big[[1]]
y <- data.big[[2]]
beta.true = data.big[[3]]



#=======================
## @knitr code_part3
# Section 3: OLS with pre-built commands
#=======================

#-------------------------
# standard linear model
#--------------------------
reg1 = lm( y ~ x -1) # weights and na.action options can be important
# playing around
reg1
names(reg1)
regout1 <- summary(reg1)
regout1 
names(regout1)
regout1$coefficients
xtable(regout1)
plot(reg1)
abline(reg1)

#could also run:
summary( lm( y ~ x + I(X^2) -1  )) # I() means take input "as is" , ie, don't overlaod to equatiol notaion meaning of "^"

#---------------------------------
# heteroskedastic "robust" standard errors:
#---------------------------------
reg1.robust <- rlm(y ~ x -1)
regout1.robust <-summary(reg1.robust)
regout1.robust
plot(reg1.robust)
abline(reg1.robust)

#----------------------------
# Testing linear hypotheses
#-----------------------------
linearHypothesis(reg1, c("x1=0", "x2=0"))

#-----------------------
# fitted values
#-----------------------
reg1.hat        <- fitted(reg1)
reg1.robust.hat <- fitted(reg1.robust)
reg1.resids     <- residuals(reg1)

#------------------------
# Dummy regression with "factor" variable types (i.e., categorical variables)
#------------------------
types <- c("blue", "green","red")
type  <- sample(types,length(y),replace=T)
reg1.dummies <- lm(y ~ x[,2:dim(x)[2]] + type -1)
summary(reg1.dummies)
reg1.dummies2 <- lm(y ~ x[,2:dim(x)[2]] + type +  x[,2]*type -1)
summary(reg1.dummies2)

#-------------------
# A very short summary table introducing tapply
#-------------------

# equivalent to "bys ___ : in STATA:
tapply(X[,2],type,mean) # assuming row defines individual!
tapply(X[,2],type,summary) # assuming row defines individual!


#Diagnostics for regressions
reg0 <- lm(Y ~ X[,2] + X[,3])
summary(reg0) 

avPlots(reg0, id.n=2,id.cex=.7) # partial regression / influence plots
qqPlot(reg0, id.n=2)
influenceIndexPlot(reg0, id.n=2)
influencePlot(reg0, id.n=2)
ncvTest(reg0) # test for heteroskedasticity (homosk data)
ncvTest(reg1) # test for heteroskedasticity (heterok data)

# A probit and a logit
Ybin <- Y>44
probit1 <- glm(Ybin ~ X - 1, family = binomial(link="probit"))
summary(probit1)
confint(probit1)

wald.test(b = coef(probit1), Sigma = vcov(probit1), Terms = 2:3)

logit1 <- glm(Ybin ~ X - 1, family = binomial(link="logit"))
summary(logit1)

# generate new data to predict on
data2 <- GenerateDataSimple(100)
X2 <- data2[[2]]
colnames
Y2 <- data2[[1]]
Y2bin <- Y2>44
logLik(probit1)
# predicted values:
data2.predvalues = predict(probit1, newdata=data2, type = "response") #returns probability of success, can add standard error with se.fit = T
cbind(Y2bin,data2.predvalues)


# partial linear regression
GenerateDataNonlinear = function(n) {   
    x1     <- rnorm(n,6,2)
    X      <- cbind(1,x1)
    theta2 <- 5
    eps    <- rnorm(n,0,sqrt(140))
    Y   <-  as.matrix(5 + 3*x1  + 1.1 * x1^2  + 30* sin(x1) - .15 * x1^3  + eps )
    colnames(X) <- c("const", "x1")
    colnames(Y) <- c("Y")
    return(list(Y,X))
}

data3 <- data.frame(GenerateDataNonlinear(1000))
reg3.lm   <- lm(Y ~ x1 , data=data3)
reg3.poly <- lm(Y ~ poly(x1, 4)  , data=data3)
reg3.bs   <- lm(Y ~  bs(x1, df = 5) , data=data3)
reg3.true <- lm( Y ~ x1  + I(x1^3) +  I( sin(x1))  + I(x1^3) , data=data3)


# cross validate the two models (plus data frames!:
data4 <- data.frame(GenerateDataNonlinear(5000))
data4$yhat.lm <- predict(reg3.lm, newdata= data4)
data4$yhat.poly <- predict(reg3.poly, newdata= data4)
data4$yhat.bs <- predict(reg3.bs, newdata= data4)
data4$yhat.true <- predict(reg3.true, newdata= data4)

# Plot showing fit. 
plot(Y ~ jitter(x1, factor = 3), pch=18, col=rgb(0.5, 0.5, 0.5, alpha=0.2), data=data4, ylim=c(-80,80), xlim = c(0,12))
vals <- cbind(data4$yhat.lm,data4$x1)
v2   <- vals[order(vals[,2]),]
lines(v2[,1] ~ v2[,2], col="red")
vals <- cbind(data4$yhat.poly,data4$x1)
v2   <- vals[order(vals[,2]),]
lines(v2[,1] ~ v2[,2], col="blue")
vals <- cbind(data4$yhat.bs,data4$x1)
v2   <- vals[order(vals[,2]),]
lines(v2[,1] ~ v2[,2], col="green")
vals <- cbind(data4$yhat.true,data4$x1)
v2   <- vals[order(vals[,2]),]
lines(v2[,1] ~ v2[,2], col="black")
# y ~ (a + b + c)^2 would model all two-way interactions

# Systems of equations (SUR):
data5 <- GenerateData(n=1000,J=6,nX=2,het=F)
Y     <- data5[[2]]
X     <- data5[[1]]
m     <- list()
m[[1]]    <- Y[,1] ~ X -1
m[[2]]    <- Y[,2] ~ X - 1
m[[3]]    <- Y[,3] ~ X - 1

lm.sur<- systemfit(m, method="SUR")


#=======================
## @knitr code_part4
# Section 4: Building OLS Commands
#=======================

# Thank you Oliver Browne for cleaning up and extending my messy code, much of what is below he helped write


# generating some fresh data:
data5 <- data.frame(GenerateDataSimple(1000))
rm(Y)
attach(data5)
X = cbind(const,x1,x2)
x = X
y = Y

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

# Formatting output for summary command on my new object class. 
summary.olsoutput <- function(object) {
    beta  <- as.matrix(object[[1]])
    se    <- as.matrix(object[[2]])
    out   <- cbind(beta,se)
    colnames(out) <- c("beta","se")
    out   <- round(out, 3)
    print("Regression Results")
    print(paste("",round(beta,3),""))
    print(paste("(",round(se,3),")", sep= ""))
    return(out)
}
    


#============================
## @knitr code_part3
# OLS Estimates
#============================

method   <- c(0,1,2,3,4)

#Generate OLS estimates for each method
sapply(method,function(m,x,y) OLS(x,y,m),x=X,y=Y,simplify="array")
beta_hat <- sapply(method,function(m,x,y) OLS(x,y,m),x=X,y=Y,simplify="array")
beta_hat



