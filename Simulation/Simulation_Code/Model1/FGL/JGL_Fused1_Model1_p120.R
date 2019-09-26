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
source(paste(func.path,"KFoldCVSplit2.R", sep="/"))
source(paste(func.path,"ROC_JGL.R", sep="/"))
source(paste(func.path,"AUC_Func.R", sep="/"))

library(Matrix)
library(MASS)
library(matrixcalc)
library(fda)
library(quadprog)
library(foreach)
library(doParallel)
library(JGL)

# Set the constants

m <- 5  ## number of function basis
p <- 120  ## dimension
n <- 100  ## sample size
obo <- 200 ## original observation number on each function

#lambda1.v <- seq(from=0, to=1, by=0.1)  # The tuning parameter for the graphical lasso penalty
lambda2.v <- seq(from=1, to=0, by=-0.05)  # The tuning parameter for the fused or group lasso penalty

################### Part 1: Generate Random Function #########################

# Decide Precision Matrix
load(paste(model.path, "/model1_p", p, ".RData", sep=""))
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

# In the following steps, we use principle.score.X and principle.score.Y to
# estimated the differential matrix

################### Part 3: Estimate the Differential Graph #########################

prin.score <- list(X=principle.score.X, Y=principle.score.Y)

# lambda1.v <- seq(from=0, to=1, by=0.01)
# n.lam1 <- length(lambda1.v)
# suppErr.v <- rep(Inf, n.lam1)
# 
# for (a in 1:n.lam1){
#   lambda1 <- lambda1.v[a]
#   
#   JGL.result <- JGL(prin.score, penalty="fused", lambda1=lambda1, lambda2=0, 
#                     penalize.diagonal=FALSE, return.whole.theta=TRUE)
#   Supp.list <- Supp_JGL(JGL.result, p, M.selected)
#   
#   errX <- abs((sum(Supp.list$X)-sum(diag(Supp.list$X)))/2 - 
#                 (sum(Theta.L$SupportX)-sum(diag(Theta.L$SupportX)))/2)
#   errY <- abs((sum(Supp.list$Y)-sum(diag(Supp.list$Y)))/2 - 
#                 (sum(Theta.L$SupportY)-sum(diag(Theta.L$SupportY)))/2)
#   if (a>1){
#     if (errX+errY > suppErr.v[a-1]){
#       lambda1.choose <- lambda1.v[a-1]
#       suppErr.v[a] <- errX+errY
#       break
#     }
#   }
#   
#   suppErr.v[a] <- errX+errY
# }

lambda1.choose <- 0

# setting parallel computing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

M.paral <- length(lambda2.v)

Paral.result <- foreach(lambda2.ind=1:M.paral, .combine='rbind') %do% {
  lambda2 <- lambda2.v[lambda2.ind]
  JGL.result <- JGL(prin.score, penalty="fused", lambda1=lambda1.choose, lambda2=lambda2, 
                    penalize.diagonal=FALSE, return.whole.theta=TRUE)
  ROC.point <- ROC_JGL(JGL.result, Theta.L$SupportDelta, p, M.selected)
  
  temp.result <- matrix(ROC.point, nrow=1)
  colnames(temp.result) <- c("FPR", "TPR")
  temp.result
}

Paral.result <- rbind(Paral.result, matrix(c(1,1), nrow=1))

stopCluster(cl)

save(Paral.result, file=paste(save.path, "/Para_result_p", p, 
                              "_runind", run.ind, "_lam1_", lambda1.choose, ".RData", sep=""))

# # test
# load(paste(save.path, "/Para_result_p", p, 
#            "_runind", run.ind, "_lam1_", lambda1.choose, ".RData", sep=""))
# 
# plot(Paral.result[, 1], Paral.result[, 2], type="b")
# AUC_Func(Paral.result[, 1], Paral.result[, 2])


lambda1.choose <- 0.01

# setting parallel computing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

M.paral <- length(lambda2.v)

Paral.result <- foreach(lambda2.ind=1:M.paral, .combine='rbind') %do% {
  lambda2 <- lambda2.v[lambda2.ind]
  JGL.result <- JGL(prin.score, penalty="fused", lambda1=lambda1.choose, lambda2=lambda2, 
                    penalize.diagonal=FALSE, return.whole.theta=TRUE)
  ROC.point <- ROC_JGL(JGL.result, Theta.L$SupportDelta, p, M.selected)
  
  temp.result <- matrix(ROC.point, nrow=1)
  colnames(temp.result) <- c("FPR", "TPR")
  temp.result
}

Paral.result <- rbind(Paral.result, matrix(c(1,1), nrow=1))

stopCluster(cl)

save(Paral.result, file=paste(save.path, "/Para_result_p", p, 
                              "_runind", run.ind, "_lam1_", lambda1.choose, ".RData", sep=""))

# # test
# load(paste(save.path, "/Para_result_p", p, 
#            "_runind", run.ind, "_lam1_", lambda1.choose, ".RData", sep=""))
# 
# plot(Paral.result[, 1], Paral.result[, 2], type="b")
# AUC_Func(Paral.result[, 1], Paral.result[, 2])

timeend <- proc.time()
runningtime <- timeend-timestart
print(runningtime)