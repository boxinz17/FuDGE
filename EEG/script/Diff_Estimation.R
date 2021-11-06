####################### Part 1 Preparation #############################################
file.path <- paste("E:/Research/Projects/FuDGE/Code/EEG", sep="")
func.path <- paste("E:/Research/Projects/FuDGE/Code/EEG/Function_Definitions", sep="")
# Load the data, library and self-defined functions
library(Matrix)
library(MASS)
library(matrixcalc)
library(fda)
library(quadprog)
library(poweRlaw)

load(paste(file.path, "/alco_filtered_array.Rdata", sep=""))  # load alco.filtered.array
load(paste(file.path, "/contrl_filtered_array.Rdata", sep=""))  # load contrl.filtered.array

source(paste(func.path,"ProxAlg_DFGM.R", sep="/"))
source(paste(func.path,"KFoldCVSplit2.R", sep="/"))
source(paste(func.path,"Graph_Genertation.R", sep="/"))
source(paste(func.path,"ROC_plot2.R", sep="/"))

N.alco <- 77
N.contrl <- 45
p <- 64
N.sample <- 256

u.alco <- seq(1/N.sample, 1, 1/N.sample)
u.contrl <- seq(1/N.sample, 1, 1/N.sample)

# set paramters for old code
observ.X <- alco.filtered.array
observ.Y <- contrl.filtered.array
u.X <- u.alco
u.Y <- u.contrl
obo <- N.sample

## As a conclude, what we have for analysis in the following steps are 
## u.x, u.y, observ.X, observ.Y

################### Part 2: Basis Expansion and FPCA ########################

# Select the number of basis functions (dimension of B-Splines) L

# First for Alcoholic Group
n <- N.alco
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

# Then for Control Group
n <- N.contrl
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

principle.score.X <- matrix(0, nrow=n, ncol=1)
principle.score.Y <- matrix(0, nrow=n, ncol=1)
for (j in 1:p){
  obs.val.matrix.X <- matrix(rep(0, (n * obo)), nrow=obo)
  obs.val.matrix.Y <- matrix(rep(0, (n * obo)), nrow=obo)
  for (i in c(1:n)){
    obs.val.matrix.X[, i] <- as.vector(observ.X[i, j, ])
    obs.val.matrix.Y[, i] <- as.vector(observ.Y[i, j, ])
  }
  bbasis <- create.bspline.basis(rangeval=c(0, 1), nbasis=L.selected)
  
  fd.object.array.X <- Data2fd(argvals=u.X, y=obs.val.matrix.X, basisobj=bbasis)
  fd.object.array.Y <- Data2fd(argvals=u.X, y=obs.val.matrix.Y, basisobj=bbasis)
  
  fd.pca.obj.X <- pca.fd(fd.object.array.X, nharm=M.selected)
  fd.pca.obj.Y <- pca.fd(fd.object.array.Y, nharm=M.selected)
  fd.pca.obj <- pca.fd(fd.object.array.X+fd.object.array.Y, nharm=M.selected)
  
  proj.score.X <- inprod(fd.object.array.X, fd.pca.obj$harmonics, rng=c(0,1))
  proj.score.Y <- inprod(fd.object.array.Y, fd.pca.obj$harmonics, rng=c(0,1))
  
  principle.score.X <- cbind(principle.score.X, proj.score.X)
  principle.score.Y <- cbind(principle.score.Y, proj.score.Y)
}
principle.score.X <- principle.score.X[,-1]
principle.score.Y <- principle.score.Y[,-1]

# Use estimated principle components scores matrix to estimate covariance matrix of
# principle components, both for X and Y
principle.score.X.cen <- scale(principle.score.X, center=T, scale=F)
principle.score.Y.cen <- scale(principle.score.Y, center=T, scale=F)

estimated.cov.pc.X <- (t(principle.score.X.cen) %*% principle.score.X.cen) / (N.alco-1)
estimated.cov.pc.Y <- (t(principle.score.Y.cen) %*% principle.score.Y.cen) / (N.contrl-1)

# In the following steps, we use estimated.cov.pc.X and estimated.cov.pc.Y to
# estimated the differential matrix

################### Part 3: Estimate the Differential Graph #########################
timestart <- proc.time()
s.rate <- 0.01
num.change.exp <- floor(s.rate * (p*(p+1)/2))  # expected number of changed edges

lambda.min <- 0
lambda.max <- max(max(abs(estimated.cov.pc.X)), max(abs(estimated.cov.pc.Y)))

count <- 0
for (lambda in seq(lambda.min, lambda.max, by=1e-4)){
  count <- count + 1
  g <- ProxAlg_DFGM(cov.X=estimated.cov.pc.X, cov.Y=estimated.cov.pc.Y, 
               p=p, M=M.selected, lambda=lambda)
  print(paste("Convergence:", g$converge, "Number of iterations: ", g$num.iter))
  FrobMat <- g$blockFrob
  AdjMat <- matrix(0, nrow=p, ncol=p)
  AdjMat[which(FrobMat > 0)] <- 1
  diag(AdjMat) <- 0
  est.num.ch <- (sum(AdjMat) - sum(diag(AdjMat))) / 2
  print(est.num.ch)
  if (est.num.ch <= num.change.exp){
    save(AdjMat, file=paste(file.path, "/AdjMat.Rdata", sep=""))
    break
  }
}

timeend <- proc.time()
runningtime <- timeend-timestart
round(as.numeric(runningtime[3]), 2)
