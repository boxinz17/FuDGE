# This is the model1 from part 5 of Tony Cai paper: 
# http://www-stat.wharton.upenn.edu/~tcai/paper/Testing-Differential-Network.pdf
# modifed to functional graph case

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

Graph_Generation_Model2 <- function(p, m, s=4, seed=233){
  set.seed(seed=seed)
  
  Omega.star <- diag(1, nrow=(p*m), ncol=(p*m))
  for (i in 1:(p-1)){
    Omega.star[((i-1)*m+1):(i*m), (i*m+1):((i+1)*m)] <- diag(0.4, nrow=m, ncol=m)
    Omega.star[(i*m+1):((i+1)*m), ((i-1)*m+1):(i*m)] <- diag(0.4, nrow=m, ncol=m)
  }
  
  for (i in 1:(p-2)){
    Omega.star[((i-1)*m+1):(i*m), ((i+1)*m+1):((i+2)*m)] <- diag(0.2, nrow=m, ncol=m)
    Omega.star[((i+1)*m+1):((i+2)*m), ((i-1)*m+1):(i*m)] <- diag(0.2, nrow=m, ncol=m)
  }
  
  support.delta <- matrix(0, ncol=p, nrow=p)
  U <- matrix(0, nrow=(p*m), ncol=(p*m))
  for (i in 1:s){
    support.delta[i, i+3] <- 1
    support.delta[i+3, i] <- 1
    
    c <- runif(1, min=0.6, max=1)
    
    if (p <= 30){
      c <- c / 2
    } else if (p <= 60){
      c <- c / 3
    } else if (p <= 90){
      c <- c / 4
    } else {
      c <- c / 5
    }
    
    U[((i-1)*m+1):(i*m), ((i+2)*m+1):((i+3)*m)] <- diag(c, nrow=m, ncol=m)
    U[((i+2)*m+1):((i+3)*m), ((i-1)*m+1):(i*m)] <- diag(c, nrow=m, ncol=m)
  }
  
  Pre.mat.X <- Omega.star
  Pre.mat.Y <- Pre.mat.X + U
  
  support.X <- matrix(0, ncol=p, nrow=p)
  support.Y <- matrix(0, ncol=p, nrow=p)
  for (i in 1:p){
    for (j in 1:p){
      FnormX <- frobenius.norm(Pre.mat.X[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)])
      FnormY <- frobenius.norm(Pre.mat.Y[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)])
      
      if (FnormX > 0){
        support.X[i, j] <- 1
        support.X[j, i] <- 1
      }
      
      if (FnormY > 0){
        support.Y[i, j] <- 1
        support.Y[j, i] <- 1
      }
    }
  }
  
  return(list(X=Pre.mat.X, Y=Pre.mat.Y, SupportX=support.X, 
              SupportY=support.Y, SupportDelta=support.delta))
}




