#' @param p number of nodes
#' @param m number of function basis to generate random function
#' @param n number of observations
#' @param s number of changed edges
#' @param seed random seed
#' 
#' @return List of precision matrix for graph X and Y. 
#' 

library(expm)
library(matrixcalc)

Graph_Generation_Model_Examp <- function(m, alpha, beta, seed=996){
  set.seed(seed=seed)
  
  scale.param <- 1.0
  
  Omega.X.eps <- matrix(0, nrow=(3*m), ncol=(3*m))
  for (i in 1:3){
    for (k in 1:m){
      Omega.X.eps[(i-1)*m+k, (i-1)*m+k] = k**alpha / scale.param
    }
  }
  
  Omega.Y.eps <- Omega.X.eps
  for (i in 1:3){
    Omega.Y.eps[i*m, i*m] <- 1 / ((1/m)**alpha + (1/m)**beta) / scale.param
  }
  
  A = diag(1, nrow=(3*m), ncol=(3*m))
  A[1:m, (m+1):(2*m)] <- - diag(1, nrow=m, ncol=m)
  A[(2*m+1):(3*m), (m+1):(2*m)] <- - diag(1, nrow=m, ncol=m)
  
  Pre.mat.X = t(A) %*% Omega.X.eps %*% A
  Pre.mat.Y = t(A) %*% Omega.Y.eps %*% A
  
  support.X <- matrix(0, ncol=3, nrow=3)
  support.Y <- matrix(0, ncol=3, nrow=3)
  support.delta <- matrix(0, nrow=3, ncol=3)
  for (i in 1:3){
    for (j in 1:3){
      FnormX <- frobenius.norm(Pre.mat.X[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)])
      FnormY <- frobenius.norm(Pre.mat.Y[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)])
      Fnormdelta <- frobenius.norm(Pre.mat.X[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)]-frobenius.norm(Pre.mat.Y[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)]))
      
      if (FnormX > 0){
        support.X[i, j] <- 1
        support.X[j, i] <- 1
      }
      
      if (FnormY > 0){
        support.Y[i, j] <- 1
        support.Y[j, i] <- 1
      }
      
      if (Fnormdelta > 0){
        support.delta[i, j] <- 1
        support.delta[j, i] <- 1
      }
    }
  }
  
  
  return(list(X=Pre.mat.X, Y=Pre.mat.Y, SupportX=support.X, 
              SupportY=support.Y, SupportDelta=support.delta))
}
