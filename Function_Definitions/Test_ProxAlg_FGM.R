################## Part 0: Preparation #########################################

# Load the library and functions
source("E:/学习工作/硕士/graduate research/code/Differential Functional Graphical Model/Function_Definitions/Graph_Genertation_PowLaw.R")
source("E:/学习工作/硕士/graduate research/code/Differential Functional Graphical Model/Function_Definitions/ProxAlg_FGM.R")


library(Matrix)
library(MASS)
library(matrixcalc)
library(fda)
library(quadprog)

# Set the constants

m <- 5  ## number of function basis
p <- 30  ## dimension
n <- 100  ## sample size
obo <- 200 ## original observation number on each function

gammamax <- 10
gammamin <- 0
steplength <- 0.1
gamma.v <- seq(gammamax, gammamin, -steplength)

################### Part 1: Generate Random Function #########################

# Decide Precision Matrix
Theta.true <- diag(1, nrow=(p*m), ncol=(p*m))
Support.True <- diag(1, nrow=p, ncol=p)
for (i in 1:(p-1)){
  Support.True[i, i+1] <- 1
  Support.True[i+1, i] <- 1
  Theta.true[((i-1)*m+1):(i*m), (i*m+1):((i+1)*m)] <- diag(0.4, nrow=m, ncol=m)
  Theta.true[(i*m+1):((i+1)*m), ((i-1)*m+1):(i*m)] <- diag(0.4, nrow=m, ncol=m)
}

# generate the covariance matrix
Covariance.X <- solve(Theta.true)

# generate Gaussian random vector with zero mean and Covariance matrix
set.seed(666)
X <- mvrnorm(n, mu=rep(0, (m * p)), Sigma=Covariance.X)

# generate the observations
u.X <- seq(1/obo, 1, 1/obo)  ## a vector recording the obervation time points of X

fbasis <- create.fourier.basis(nbasis = m)  ## set the fourier basis

observ.X <- array(rep(0, (n * p * obo)), c(n, p, obo))  ## array recording X observations

for (i in c(1:n)){
  for (j in c(1:p)){
    coeff.X <- X[i, (1 + (j - 1) * m):(j * m)]
    fdatabasis.X <- fd(coeff.X, fbasis)
    observ.X[i, j, ] <- eval.fd(u.X, fdatabasis.X) + rnorm(obo, 0, 0.5)
  }
}

## As a conclude, what we have for analysis in the following steps are 
## u.x,observ.X

################### Part 2: Basis Expansion and FPCA #########################

# set the L and M for X and Y be the same
L.selected <- 6
M.selected <- 5

# Use selected L and M to estimate the principal component scores matrix

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

# Use estimated principle components scores matrix to estimate covariance matrix of
# principle components
principle.score.X.cen <- scale(principle.score.X, center=T, scale=F)

estimated.cov.pc.X <- (t(principle.score.X.cen) %*% principle.score.X.cen) / (n-1)

# In the following steps, we use estimated.cov.pc.X

################### Part 3: Estimate the Differential Graph #########################

temp.result.X <- matrix(NA, nrow=length(gamma.v), ncol=2)
colnames(temp.result.X) <- c("TPR", "FPR")
count <- 0
for (gam in gamma.v){
  count <- count + 1
  ## Set the matrix representing true differential edges
  Supp.Mat.X <- Support.True
  gamma.choose <- gam
  
  print(paste("gamma is:", gamma.choose))
  
  gX <- ProxAlg_FGM(S=estimated.cov.pc.X, p=p, M=M.selected, 
                    gamma=gamma.choose)
  print(gX$converge)
  print(gX$num.iter)
  # Compute TPR and FPR
  FrobMat <- gX$blockFrob
  
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat[which(FrobMat > 0)] <- 1
  
  v.estimate <- c()
  v.true <- c()
  for (i in 1:(p - 1)){
    v.estimate <- c(v.estimate, Edgeshat[i, (i + 1):p])
    v.true <- c(v.true, Supp.Mat.X[i, (i + 1):p])
  }
  
  n <- length(v.estimate)
  TP <- rep(0, n)
  TN <- rep(0, n)
  FP <- rep(0, n)
  FN <- rep(0, n)
  
  TP <- as.double(sum(v.estimate & v.true))
  TN <- as.double(sum(!v.estimate & !v.true))
  FP <- as.double(sum(v.estimate & !v.true))
  FN <- as.double(sum(!v.estimate & v.true))
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  
  temp.result.X[count, 1] <- TPR
  temp.result.X[count, 2] <- FPR
}

plot(temp.result.X[, 2], temp.result.X[, 1], xlab="FPR", ylab="TPR", main="ROC Curve", type="b")
