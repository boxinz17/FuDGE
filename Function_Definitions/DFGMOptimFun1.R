source("E:/学习工作/硕士/graduate research/code/Differential Functional Graphical Model/experiment_1/Function_Definitions/ConstrainedL2.R")
source("E:/学习工作/硕士/graduate research/code/Differential Functional Graphical Model/experiment_1/Function_Definitions/MatOrg.R")
source("E:/学习工作/硕士/graduate research/code/Differential Functional Graphical Model/experiment_1/Function_Definitions/IRGL2.R")


DFGMOptimFun1 <- function(SigmahatX, SigmahatY, lambda, gamma, p, M, n.iteration=100, stop.criterion=0.05){
  # Optimzation function to solve the optimization problem in DFGM, using L2 group loss with constraint
  # 
  # Args:
  #   SigmahatX: The estimated covariance matrix of Graph X
  #   SigmahatY: The estimated covariance matrix of Graph Y
  #   lambda: The hyperparamter used in constraint
  #   gamma: The hyperparameter used in identifying edges
  #   p: The number of vertices
  #   M: The number of principle components
  #   n.iteration: The maximum number of iterations for solving optimization problem
  #   stop.criterion: The ratio of lambda to stop the iteration
  # Returns:
  #   A list with:
  #     DelatMathat: The estimated differential matrix between ThetaX and ThetaY
  #     blockFrob: A matrix of Frobenius norm of each block matrix
  #     Edges: A matrix indicates the estimated edges
  g1 <- MatOrg(SigmahatX, SigmahatY, p, M)
  g2 <- IRGL2(g1$matrixList, g1$vector, lambda, p, M, n.iteration, stop.criterion)
  
  Delta <- matrix(0, nrow=(p * M), ncol=(p * M))
  for (i in 1:p){
    for (j in 1:p){
      Delta[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 
        matrix(g2$betahat[i, j, ], nrow=M, byrow=FALSE)
    }
  }
  ## Symmertize the Delta
  for (i in 1:(p *M)){
    for(j in 1:(p * M)){
      if(abs(Delta[i, j]) < abs(Delta[j, i])){
        Delta[j, i] <- Delta[i, j]
      } else {
        Delta[i, j] <- Delta[j, i]
      }
    }
  }
  
  Delta.frob <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      Delta.frob[i, j] <- sqrt(sum(g2$betahat[i, j, ] ^ 2))
    }
  }
  ## Symmetrize the Delta.frob
  for (i in 1:p){
    for(j in 1:p){
      if((Delta.frob[i, j]) < Delta.frob[j, i]){
        Delta.frob[j, i] <- Delta.frob[i, j]
      } else {
        Delta.frob[i, j] <- Delta.frob[j, i]
      }
    }
  }
  
  EdgeMat <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      if (sqrt(sum(g2$betahat[i, j, ] ^ 2)) > gamma){
        EdgeMat[i, j] <- 1
      } else {
        EdgeMat[i, j] <- 0
      }
    }
  }
  ## Symmetrize the EdgeMat
  for (i in 1:p){
    for (j in 1:p){
      if((EdgeMat[i, j] == 0) | (EdgeMat[j, i] == 0)){
        EdgeMat[i, j] <- 0
        EdgeMat[j, i] <- 0
      }
    }
  }
  
  return(list(DelatMathat=Delta, blockFrob=Delta.frob, Edges=EdgeMat))
}

# test
#library(Matrix)
#library(MASS)
#library(matrixcalc)
#library(fda)
#library(quadprog)

## Set the constants

#M <- 5  ## number of function basis
#p <- 3  ## dimension

## Decide Precision Matrix
#Theta.X <- diag(1, M * p)
#Theta.X[(M + 1):(2 * M), 1:M] <- diag(0.9, M)
#Theta.X[1:M, (M + 1):(2 * M)] <- diag(0.9, M)

#Theta.Y <- diag(1, M * p)

#Delta.true <- Theta.X - Theta.Y

#CovMatX <- solve(Theta.X)
#CovMatY <- solve(Theta.Y)

#g <- DFGMOptimFun1(CovMatX, CovMatY, 0.1, 0.01, p, M)