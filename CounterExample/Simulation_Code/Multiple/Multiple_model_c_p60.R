# path setup

func.path <- ""  # function definition path
save.path <- ""  # result saving path
model.path <- ""  # model saving path

################## Part 0: Preparation #########################################

### Reads in arguments passed in from command line
### args is a vector of strings 
args <- (commandArgs(TRUE))

for(i in 1:length(args)){
  # run each element of the vector as if passed directly to the console
  # have run.ind
  eval(parse(text = args[[i]]))
}

timestart <- proc.time()

# Load the library and functions
source(paste(func.path,"ProxAlg_DFGM.R", sep="/"))
source(paste(func.path,"KFoldCVSplit2.R", sep="/"))

library(Matrix)
library(MASS)
library(matrixcalc)
library(fda)
library(quadprog)
library(foreach)
library(doParallel)

# Set the constants

m <- 5  ## number of function basis
p <- 60  ## dimension
n <- 100  ## sample size
obo <- 15 ## original observation number on each function

lambdamax <- 15
lambdamin <- 0
steplength <- 0.2
lambda.v <- seq(lambdamax, lambdamin, -steplength)

M.paral <- length(lambda.v)

################### Part 1: Generate Random Function #########################

# Decide Precision Matrix
load(paste(model.path, "/model_c_p", p, ".RData", sep=""))
Theta.L <- g

Theta.X <- Theta.L$X
Theta.Y <- Theta.L$Y

Delta.true <- Theta.X - Theta.Y

# generate the covariance matrix
Covariance.X <- solve(Theta.X)
Covariance.Y <- solve(Theta.Y)

# generate Gaussian random vector with zero mean and Covariance matrix
set.seed(run.ind)
X <- mvrnorm(n, mu=rep(0, (m * p)), Sigma=Covariance.X)
Y <- mvrnorm(n, mu=rep(0, (m * p)), Sigma=Covariance.Y)

# generate the observations
u.X <- seq(1/obo, 1, 1/obo)  ## a vector recording the obervation time points of X
u.Y <- seq(1/obo, 1, 1/obo)  ## a vector recording the obervation time points of Y

fbasis <- create.fourier.basis(nbasis = m)  ## set the fourier basis

observ.X <- array(rep(0, (n * p * obo)), c(n, p, obo))  ## array recording X observations
observ.Y <- array(rep(0, (n * p * obo)), c(n, p, obo))  ## array recording Y observations

for (i in c(1:n)){
  for (j in c(1:p)){
    coeff.X <- X[i, (1 + (j - 1) * m):(j * m)]
    fdatabasis.X <- fd(coeff.X, fbasis)
    observ.X[i, j, ] <- eval.fd(u.X, fdatabasis.X) + rnorm(obo, 0, 0.5)
    
    coeff.Y <- Y[i, (1 + (j - 1) * m):(j * m)]
    fdatabasis.Y <- fd(coeff.Y, fbasis)
    observ.Y[i, j, ] <- eval.fd(u.Y, fdatabasis.Y) + rnorm(obo, 0, 0.5)
  }
}

## As a conclude, what we have for analysis in the following steps are 
## u.x, u.y, observ.X, observ.Y

# setting parallel computing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

################### Part 2: Estimate the differential graphs #########################
Paral.result <- foreach(lambda.ind=1:M.paral, .combine='rbind') %do% {
  lambda.choose <- lambda.v[lambda.ind]
  
  result.L <- list()
  for (t in 1:obo){
    X <- scale(observ.X[, , t], center=T, scale=F)
    SX <- (t(X) %*% X) / (n-1)
    
    Y <- scale(observ.Y[, , t], center=T, scale=F)
    SY <- (t(Y) %*% Y) / (n-1)
    
    g <- ProxAlg_DFGM(SX, SY, p, M=1, lambda=lambda.choose)
    g$lambda <- lambda.choose
    result.L[[t]] <- g
  }
  result <- list(Deltas=result.L)
  result
}

# stop cluster
stopCluster(cl)

timeend <- proc.time()
runningtime <- timeend-timestart
print(runningtime)

# save result
save(Paral.result,
     file=paste(save.path, "/Sep_Time_Est_p", p, "_runind", run.ind, ".RData", sep=""))
