# This is the model3 from part 5 of Tony Cai paper: 
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

rescale_func <- function(c, p){
  if (p <= 30){
    c <- c / 2
  } else if (p <= 60){
    c <- c / 3
  } else if (p <= 90){
    c <- c / 4
  } else {
    c <- c / 5
  }
  return(c)
}

Graph_Generation_Model3 <- function(p, m, s=4, seed=667){
  set.seed(seed=seed)
  
  Omega.star <- diag(1, nrow=(p*m), ncol=(p*m))
  support.X <- matrix(0, nrow=p, ncol=p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      c <- 0.1 * rbinom(1, 1, 0.3)
      if (c > 0){
        support.X[i, j] <- 1
        support.X[j, i] <- 1

        Omega.star[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] <- diag(c, nrow=m, ncol=m)
        Omega.star[((j-1)*m+1):(j*m), ((i-1)*m+1):(i*m)] <- diag(c, nrow=m, ncol=m)
      }
    }
  }
  
  Omega <- Omega.star
  
  # for (i in 1:p){
  #   for (j in 1:p){
  #     if (i == j){
  #       next
  #     }
  #     
  #     z1 <- rbinom(1, 1, 0.6)
  #     if (z1 > 0){
  #       
  #       z2 <- rbinom(1, 1, 0.5)
  #       if (z2 > 0){
  #         c <- runif(1, 0.2, 0.5)
  #       } else {
  #         c <- runif(1, -0.5, -0.2)
  #       }
  #       
  #       c <- rescale_func(c, p)
  #       
  #       if (support.X[i, j] == 0){
  #         support.X[i, j] <- 1
  #         support.X[j, i] <- 1
  #       }
  #       
  #       Omega.star[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] <- diag(c, nrow=m, ncol=m)
  #     }
  #   }
  # }
  # 
  # Omega <- 0.5 * (Omega.star + t(Omega.star))
  
  support.delta <- matrix(0, ncol=p, nrow=p)
  U <- matrix(0, nrow=(p*m), ncol=(p*m))
  num.change <- 0
  while(num.change < s){
    i <- sample(1:(p-1), 1)
    if (i == (p-1)){
      j <- p
    } else {
      j <- sample((i+1):p, 1)
    }
    
    if(support.X[i, j]==1){
      next
    }
    
    support.delta[i, j] <- 1
    support.delta[j, i] <- 1
    num.change <- sum(support.delta)/2
    
    c <- runif(1, 0.6, 1)
    
    c <- rescale_func(c, p)
    
    U[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] <- diag(c, nrow=m, ncol=m)
    U[((j-1)*m+1):(j*m), ((i-1)*m+1):(i*m)] <- diag(c, nrow=m, ncol=m)
    
  }
  
  Pre.mat.X <- Omega
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
