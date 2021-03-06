#' Genertae Precision Matrix for graph X and Y
#' 
#' The basic setting is following TT Cai paper, but modified to a functional setting.
#' 
#' @param p number of nodes
#' @param m number of function basis to generate random function
#' @param alpha power parameter for power-law degree distribution, the larger the sparser
#' @param seed random seed
#' 
#' @return List of precision matrix for graph X and Y. 
#' 
#' @examples
#' Preci.Mat <- Graph_Generation(40, 5, 2, 111)
#' Preci.Mat$X
#' Preci.Mat$Y

library(poweRlaw)
library(matrixcalc)

Graph_Generation <- function(p, m, alpha, seed){
  set.seed(seed)
  
  # Generate Precision matrix for graph X
  num.edg <- 0
  ## randomly choose the support for Graph X. See TT Cai paper.
  support.X <- matrix(0, nrow=p, ncol=p)  # support matrix for graph X
  for (i in 1:p){
    support.X[i, i] <- 1
  }
  for (i in 1:(p-1)){
    d <- rpldis(1, xmin=1, alpha=alpha)
    if (d > p-i){
      d <- p - i
    }
    neigh <- sample((i+1):p, size=d, replace=F)
    
    for (k in 1:d){
      support.X[i, neigh[k]] <- 1
      support.X[neigh[k], i] <- 1
    }
    
    num.edg <- num.edg + d
    if (num.edg > p*(p-1)/10){
      break
    }
  }
  edg.num.v <- apply(support.X, 1, sum) - 1
  
  neighbor.list <- list()
  for (i in 1:p){
    if (edg.num.v[i] == 0){
      neighbor.list[[i]] <- "No Edge"
      next
    }
    neighbors <- c()
    for (j in 1:p){
      if (j != i){
        if (support.X[i, j] > 0){
          neighbors <- c(neighbors, j)
        }
      }
    }
    neighbor.list[[i]] <- neighbors
  }
  
  ## Generate precision matrix
  Pre.mat.X <- matrix(0, nrow=(m*p), ncol=(m*p))
  
  for (i in 1:p){
    if (length(neighbor.list[[i]])==1 & neighbor.list[[i]][1] == "No Edge"){
      next
    }
    for (j in neighbor.list[[i]]){
      a <- sample(c(-1, 1), size=1)
      if (a > 0){
        b <- runif(1, min = 0.2, max = 0.5)
      } else {
        b <- runif(1, min = -0.5, max = -0.2)
      }
      Pre.mat.X[((i - 1) * m + 1):(i * m), ((j - 1) * m + 1):(j * m)] <- diag(b, m)
    }
  }
  Pre.mat.X <- 0.5 * (Pre.mat.X + t(Pre.mat.X))  # average to get symmetric matrix
  
  ### Adjust diagonal element to make sure of postive definite
  for (g in 1:10){
    for (i in 1:(m*p)){
      Pre.mat.X[i, i] <- g
    }
    if (is.positive.definite(Pre.mat.X)){
      break
    }
  }
  
  # Generate Precision matrix for graph Y. Change top 2 node hub.
  Pre.mat.Y <- Pre.mat.X
  support.Y <- support.X # support matrix for graph Y, which should be identical to X
  support.delta <- matrix(0, nrow=p, ncol=p)  # support matrix for Delta
  
  ## Node 1
  i <- order(edg.num.v, decreasing=T)[1]
  
  magnitude.v <- numeric(edg.num.v[[i]]) # select top 20%, in magnitude, edges to change
  
  count <- 0
  for (j in neighbor.list[[i]]){
    count <- count + 1
    magnitude.v[count] <- sum(abs(Pre.mat.X[((i - 1) * m + 1):(i * m), 
                                            ((j - 1) * m + 1):(j * m)]))
  }
  
  edg.change <- neighbor.list[[i]][order(magnitude.v, decreasing=T)[1:floor(0.2*edg.num.v[i])]]
  
  for (j in edg.change){
    Pre.mat.Y[((i - 1) * m + 1):(i * m), ((j - 1) * m + 1):(j * m)] <- 
      -Pre.mat.X[((i - 1) * m + 1):(i * m), ((j - 1) * m + 1):(j * m)]
    Pre.mat.Y[((j - 1) * m + 1):(j * m), ((i - 1) * m + 1):(i * m)] <- 
      -Pre.mat.X[((j - 1) * m + 1):(j * m), ((i - 1) * m + 1):(i * m)]
    support.delta[i, j] <- 1
    support.delta[j, i] <- 1
  }
  
  ## Node 2
  i <- order(edg.num.v, decreasing=T)[2]
  
  magnitude.v <- numeric(edg.num.v[[i]]) # select top 20%, in magnitude, edges to change
  
  count <- 0
  for (j in neighbor.list[[i]]){
    count <- count + 1
    magnitude.v[count] <- sum(abs(Pre.mat.X[((i - 1) * m + 1):(i * m), 
                                            ((j - 1) * m + 1):(j * m)]))
  }
  
  edg.change <- neighbor.list[[i]][order(magnitude.v, decreasing=T)
                                   [1:floor(0.2*edg.num.v[i])]]
  
  for (j in edg.change){
    Pre.mat.Y[((i - 1) * m + 1):(i * m), ((j - 1) * m + 1):(j * m)] <- 
      -Pre.mat.X[((i - 1) * m + 1):(i * m), ((j - 1) * m + 1):(j * m)]
    Pre.mat.Y[((j - 1) * m + 1):(j * m), ((i - 1) * m + 1):(i * m)] <- 
      -Pre.mat.X[((j - 1) * m + 1):(j * m), ((i - 1) * m + 1):(i * m)]
    support.delta[i, j] <- 1
    support.delta[j, i] <- 1
  }
  
  return(list(X=Pre.mat.X, Y=Pre.mat.Y, SupportX=support.X, 
              SupportY=support.Y, SupportDelta=support.delta))
}


