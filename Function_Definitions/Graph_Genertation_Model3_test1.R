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


rescale_func <- function(x, p){
  if (p <= 30){
    x <- x
  } else if (p <= 60){
    x <- x / 2
  } else if (p <= 90){
    x <- x / 4
  } else if (p <= 120){
    x <- x / 8
  }
  return(x)
}

Graph_Generation_Model3 <- function(p, m, n, s=4, seed=2333){
  set.seed(seed=seed)
  
  D <- matrix(0, nrow=(p*m), ncol=(p*m))
  D.diag.v <- runif(n=p, min=0.5, max=2.5)
  for (i in 1:p){
    D[((i-1)*m+1):(i*m), ((i-1)*m+1):(i*m)] <- diag(D.diag.v[i], nrow=m, ncol=m)
  }
  
  Omega.star <- diag(1, nrow=(p*m), ncol=(p*m))
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      c <- 0.8 * rbinom(1, 1, 0.2)
      
      c <- rescale_func(c, p)
      
      if (c > 0){
        Omega.star[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] <- diag(c, nrow=m, ncol=m)
        Omega.star[((j-1)*m+1):(j*m), ((i-1)*m+1):(i*m)] <- diag(c, nrow=m, ncol=m)
      }
    }
  }
  
  eg0 <- eigen(Omega.star)$values
  delta <- abs(eg0[length(eg0)])+0.05
  Omega <- round(sqrtm(D) %*% ((Omega.star+diag(delta, nrow=(p*m), ncol=(p*m)))/(1+delta)) %*%
                   sqrtm(D), 10)  # round to correct asymmetricity due to numerical problem

  #Omega <- round(sqrtm(D) %*% Omega.star %*% sqrtm(D), 10)
  
  #Omega <- Omega.star
  
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
    
    support.delta[i, j] <- 1
    support.delta[j, i] <- 1
    num.change <- sum(support.delta)/2
    
    # z <- sample(c(-1, 1), 1)
    # if (z < 0){
    #   c <- runif(1, min=-2*5.5*sqrt((1.5/n)), max=-5.5*sqrt((1.5/n)))
    # } else {
    #   c <- runif(1, min=5.5*sqrt((1.5/n)), max=2*5.5*sqrt((1.5/n)))
    # }
    c <- 0.8
    c <- rescale_func(c, p)
    
    U[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] <- diag(c, nrow=m, ncol=m)
    U[((j-1)*m+1):(j*m), ((i-1)*m+1):(i*m)] <- diag(c, nrow=m, ncol=m)
  }
  
  eg1 <- eigen(Omega)$values
  eg2 <- eigen(Omega+U)$values

  delta <- abs(min(eg1[length(eg1)], eg2[length(eg2)])) + 0.05

  Pre.mat.X <- Omega + diag(delta, ncol=(p*m), nrow=(p*m))
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
