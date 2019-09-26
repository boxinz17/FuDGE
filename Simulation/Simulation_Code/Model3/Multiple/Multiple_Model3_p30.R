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
source(paste(func.path,"local_sup_basis_func.R", sep="/"))
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
p <- 30  ## dimension
n <- 100  ## sample size
obo <- 15 ## original observation number on each function

lambdamax <- 15
lambdamin <- 0
steplength <- 0.1
lambda.v <- seq(lambdamax, lambdamin, -steplength)

M.paral <- length(lambda.v)

################### Part 1: Generate Random Function #########################

# Decide Precision Matrix
load(paste(model.path, "/model3_p", p, ".RData", sep=""))
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

basis.val.M.X <- local_sup_basis_func(u.X, m)
basis.val.M.Y <- local_sup_basis_func(u.Y, m)

fbasis <- create.fourier.basis(nbasis = m)  ## set the fourier basis

observ.X <- array(rep(0, (n * p * obo)), c(n, p, obo))  ## array recording X observations
observ.Y <- array(rep(0, (n * p * obo)), c(n, p, obo))  ## array recording Y observations

for (i in c(1:n)){
  for (j in c(1:p)){
    coeff.X <- X[i, (1 + (j - 1) * m):(j * m)]
    coeff.X <- matrix(coeff.X, ncol=1)
    val.x <- basis.val.M.X %*% coeff.X
    val.x <- as.vector(val.x)
    observ.X[i, j, ] <- val.x + rnorm(obo, 0, 0.5)
    
    coeff.Y <- Y[i, (1 + (j - 1) * m):(j * m)]
    coeff.Y <- matrix(coeff.Y, ncol=1)
    val.Y <- basis.val.M.Y %*% coeff.Y
    val.Y <- as.vector(val.Y)
    observ.Y[i, j, ] <- val.Y + rnorm(obo, 0, 0.5)
  }
}

## As a conclude, what we have for analysis in the following steps are 
## u.x, u.y, observ.X, observ.Y

################### Part 2: Basis Expansion and FPCA #########################

# choose hyperparameters of Graph X & Y
# Select the number of basis L and number of principle components M via Cross Validation

# First for Graph X
potential.L.X <- 5:10  ## candidates for L of X
potential.M.X <- 4:6  ## candidates for M of X
k.cv <- 5  ## number of folds for Cross Validation

RSS.record.X <- matrix(rep(Inf, (length(potential.L.X ) * length(potential.M.X))), 
                       nrow=length(potential.L.X ))  ## record the RSS via matrix

L.num <- 0
M.num <- 0
for (L in potential.L.X){
  L.num <- L.num + 1
  M.num <- 0
  for (M in potential.M.X){
    if (M > L){
      M.num <- M.num + 1
    }
    else {
      M.num <- M.num + 1
      
      RSS.record.X[L.num, M.num] <- 0
      
      bbasis <- create.bspline.basis(rangeval=c(0, 1), nbasis=L)  ## Use bspline basis 
      ## to fit data
      for (j in 1:p){
        ## For every dimension j, we make the observations as a matrix
        obs.val.matrix <- matrix(rep(0, (n * obo)), nrow=obo)
        for (i in c(1:n)){
          obs.val.vec <- as.vector(observ.X[i, j, ])
          obs.val.matrix[, i] <- obs.val.vec
        }
        
        ## CV Split the observation matrix
        g.temp <- KFoldCVSplit2(obs.val.matrix, u.X, k.cv)
        obs.val.matrix.L <- g.temp$fval
        time.pts.L <- g.temp$t
        
        ## Fit the model and do fpca seperately, then calculate the RSS in validation set
        for (l in 1:k.cv){
          
          valid.matrix <- obs.val.matrix.L[[l]]
          valid.time <- time.pts.L[[l]]
          
          nrow.valid <- length(valid.time)
          
          if (l == 1){
            train.matrix <- obs.val.matrix.L[[2]]
            train.time <- time.pts.L[[2]]
            for (temp.count in 3:k.cv){
              train.matrix <- rbind(train.matrix, obs.val.matrix.L[[temp.count]])
              train.time <- c(train.time, time.pts.L[[temp.count]])
            }
          } else if (l == k.cv){
            train.matrix <- obs.val.matrix.L[[1]]
            train.time <- time.pts.L[[1]]
            for (temp.count in 2:(k.cv-1)){
              train.matrix <- rbind(train.matrix, obs.val.matrix.L[[temp.count]])
              train.time <- c(train.time, time.pts.L[[temp.count]])
            }
          } else {
            train.matrix <- obs.val.matrix.L[[1]]
            train.time <- time.pts.L[[1]]
            for (temp.count in c(1:(l-1), (l+1):k.cv)){
              train.matrix <- rbind(train.matrix, obs.val.matrix.L[[temp.count]])
              train.time <- c(train.time, time.pts.L[[temp.count]])
            }
          }
          
          ## As a result, we have train.matrix for trainning, valid.matrix for validation
          ## train.time records the time points in train set, 
          ## valid.time records validation set
          
          ## Train the model using train.matrix
          ### Fit with basis expansion
          fd.object.array <- Data2fd(argvals=train.time, y=train.matrix, basisobj=bbasis)
          ### FPCA
          fd.pca.obj <- pca.fd(fd.object.array, nharm=M)
          
          ## Compute RSS with Validation matrix
          predict.valid <- apply(matrix(eval.fd(fd.pca.obj$harmonics, valid.time), ncol=M), 
                                 1, sum)
          
          
          predict.valid <- matrix(predict.valid, nrow=nrow.valid) %*% matrix(rep(1, n), 
                                                                             ncol=n)
          RSS.temp <- valid.matrix - predict.valid
          RSS.record.X[L.num, M.num] <- RSS.record.X[L.num, M.num] + sum(RSS.temp ^ 2) / 
            (dim(RSS.temp)[1] * dim(RSS.temp)[2])
        }
      }
    }
  }
}

L.selected.X <- potential.L.X[which(RSS.record.X == min(RSS.record.X), 
                                    arr.ind = TRUE)[1, 1]]
M.selected.X <- potential.M.X[which(RSS.record.X == min(RSS.record.X), 
                                    arr.ind = TRUE)[1, 2]]

# Then for Graph Y
potential.L.Y <- 5:10  ## candidates for L of Y
potential.M.Y <- 4:6  ## candidates for M of Y
k.cv <- 5  ## number of folds for Cross Validation

RSS.record.Y <- matrix(rep(Inf, (length(potential.L.Y ) * length(potential.M.Y))), 
                       nrow=length(potential.L.Y ))  ## record the RSS via matrix

L.num <- 0
M.num <- 0
for (L in potential.L.Y){
  L.num <- L.num + 1
  M.num <- 0
  for (M in potential.M.Y){
    if (M > L){
      M.num <- M.num + 1
    }
    else {
      M.num <- M.num + 1
      
      RSS.record.Y[L.num, M.num] <- 0
      
      bbasis <- create.bspline.basis(rangeval=c(0, 1), nbasis=L)  ## Use bspline basis 
      ## to fit data
      for (j in 1:p){
        ## For every dimension j, we make the observations as a matrix
        obs.val.matrix <- matrix(rep(0, (n * obo)), nrow=obo)
        for (i in c(1:n)){
          obs.val.vec <- as.vector(observ.Y[i, j, ])
          obs.val.matrix[, i] <- obs.val.vec
        }
        
        ## CV Split the observation matrix
        g.temp <- KFoldCVSplit2(obs.val.matrix, u.Y, k.cv)
        obs.val.matrix.L <- g.temp$fval
        time.pts.L <- g.temp$t
        
        ## Fit the model and do fpca seperately, then calculate the RSS in validation set
        for (l in 1:k.cv){
          
          valid.matrix <- obs.val.matrix.L[[l]]
          valid.time <- time.pts.L[[l]]
          
          nrow.valid <- length(valid.time)
          
          if (l == 1){
            train.matrix <- obs.val.matrix.L[[2]]
            train.time <- time.pts.L[[2]]
            for (temp.count in 3:k.cv){
              train.matrix <- rbind(train.matrix, obs.val.matrix.L[[temp.count]])
              train.time <- c(train.time, time.pts.L[[temp.count]])
            }
          } else if (l == k.cv){
            train.matrix <- obs.val.matrix.L[[1]]
            train.time <- time.pts.L[[1]]
            for (temp.count in 2:(k.cv-1)){
              train.matrix <- rbind(train.matrix, obs.val.matrix.L[[temp.count]])
              train.time <- c(train.time, time.pts.L[[temp.count]])
            }
          } else {
            train.matrix <- obs.val.matrix.L[[1]]
            train.time <- time.pts.L[[1]]
            for (temp.count in c(1:(l-1), (l+1):k.cv)){
              train.matrix <- rbind(train.matrix, obs.val.matrix.L[[temp.count]])
              train.time <- c(train.time, time.pts.L[[temp.count]])
            }
          }
          
          ## As a result, we have train.matrix for trainning, valid.matrix for validation
          ## train.time records the time points in train set, 
          ## valid.time records validation set
          
          ## Train the model using train.matrix
          ### Fit with basis expansion
          fd.object.array <- Data2fd(argvals=train.time, y=train.matrix, basisobj=bbasis)
          ### FPCA
          fd.pca.obj <- pca.fd(fd.object.array, nharm=M)
          
          ## Compute RSS with Validation matrix
          predict.valid <- apply(matrix(eval.fd(fd.pca.obj$harmonics, valid.time), ncol=M), 
                                 1, sum)
          
          
          predict.valid <- matrix(predict.valid, nrow=nrow.valid) %*% matrix(rep(1, n), 
                                                                             ncol=n)
          RSS.temp <- valid.matrix - predict.valid
          RSS.record.Y[L.num, M.num] <- RSS.record.Y[L.num, M.num] + sum(RSS.temp ^ 2) / 
            (dim(RSS.temp)[1] * dim(RSS.temp)[2])
        }
      }
    }
  }
}

L.selected.Y <- potential.L.Y[which(RSS.record.Y == min(RSS.record.Y), 
                                    arr.ind = TRUE)[1, 1]]
M.selected.Y <- potential.M.Y[which(RSS.record.Y == min(RSS.record.Y), 
                                    arr.ind = TRUE)[1, 2]]

# set the L and M for X and Y be the same
L.selected <- ceiling((L.selected.X + L.selected.Y) / 2)
M.selected <- ceiling((M.selected.X + M.selected.Y) / 2)

# Use selected L and M to estimate the principal component scores matrix

## First for X
principle.score.X <- matrix(rep(0, n), nrow=n)
for (j in 1:p){
  obs.val.matrix <- matrix(rep(0, (n * obo)), nrow=obo)
  for (i in c(1:n)){
    obs.val.vec <- as.vector(observ.X[i, j, ])
    obs.val.matrix[, i] <- obs.val.vec
  }
  bbasis <- create.bspline.basis(rangeval=c(0, 1), nbasis=L.selected)
  fd.object.array <- Data2fd(argvals=u.X, y=obs.val.matrix, basisobj=bbasis)
  fd.pca.obj <- pca.fd(fd.object.array, nharm=M.selected)
  principle.score.X <- cbind(principle.score.X, fd.pca.obj$scores)
}
principle.score.X <- principle.score.X[, -1]

## Then for Y
principle.score.Y <- matrix(rep(0, n), nrow=n)
for (j in 1:p){
  obs.val.matrix <- matrix(rep(0, (n * obo)), nrow=obo)
  for (i in c(1:n)){
    obs.val.vec <- as.vector(observ.Y[i, j, ])
    obs.val.matrix[, i] <- obs.val.vec
  }
  bbasis <- create.bspline.basis(rangeval=c(0, 1), nbasis=L.selected)
  fd.object.array <- Data2fd(argvals=u.Y, y=obs.val.matrix, basisobj=bbasis)
  fd.pca.obj <- pca.fd(fd.object.array, nharm=M.selected)
  principle.score.Y <- cbind(principle.score.Y, fd.pca.obj$scores)
}
principle.score.Y <- principle.score.Y[, -1]

# Use estimated principle components scores matrix to estimate covariance matrix of
# principle components, both for X and Y
principle.score.X.cen <- scale(principle.score.X, center=T, scale=F)
principle.score.Y.cen <- scale(principle.score.Y, center=T, scale=F)

estimated.cov.pc.X <- (t(principle.score.X.cen) %*% principle.score.X.cen) / (n-1)
estimated.cov.pc.Y <- (t(principle.score.Y.cen) %*% principle.score.Y.cen) / (n-1)

# In the following steps, we use estimated.cov.pc.X and estimated.cov.pc.Y to
# estimated the differential matrix

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
