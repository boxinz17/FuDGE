library(matrixcalc)

blockwise_Frob <- function(ma, M){
  # ma: matrix
  # block size
  p <- dim(ma)[1]
  if (p %% M != 0){
    stop("dimension of matrix cannot be divided by block size")
  }
  n.b <- p / M
  result <- matrix(0, nrow=n.b, ncol=n.b)
  for (i in 1:n.b){
    for (j in 1:n.b){
      result[i, j] <- frobenius.norm(ma[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)])
    }
  }
  return(result)
}

object_func <- function(cov.X, cov.Y, Delta, lambda, M){
  part1 <- sum(diag(0.5 * cov.Y %*% t(Delta) %*% cov.X %*% Delta))
  part2 <- sum(diag(t(Delta) %*% (cov.Y - cov.X)))
  part3 <- lambda * sum(blockwise_Frob(Delta, M))
  return(part1-part2+part3)
}

ProxAlg_DFGM <- function(cov.X, cov.Y, p, M, lambda, Eta="Auto", n.iteration=2000, Intialization=matrix(0, nrow=(p*M), ncol=(p*M)), e_tol=1e-3){
  # Optimzation function to solve the optimization problem in DFGM
  # Use Proximal Algorithm
  # 
  # Args:
  #   cov.X: the estimated covariance matrix of Graph X
  #   cov.Y: the estimated covariance matrix of Graph Y
  #   p: the number of vertices
  #   M: the number of principle components
  #   lambda: hyperparameter.
  #   Eta: hyperparamter. If Eta="Auto", then Eta would be set as inverse of the 
  #        product of the maximum singular value of cov.X and cov.Y
  #   n.iteration: maximum times of iteration
  #
  # Returns:
  #   A list:
  #     DelatMathat: The estimated differential matrix between ThetaX and ThetaY
  #     blockFrob: A matrix of Frobenius norm of each block matrix
  if(Eta == "Auto"){
    Eta <- 1 / (svd(cov.X)$d[1] * svd(cov.Y)$d[1])
  }
  
  converge.indicator <- FALSE
  
  Delta.old <- Intialization
  for (t in 1:n.iteration){
    obj.old <- object_func(cov.X, cov.Y, Delta.old, lambda, M)
    Delta <- Delta.old
    
    # calculate gradient
    grad.Mat <- cov.X %*% Delta %*% cov.Y - (cov.Y - cov.X)
    A <- Delta - Eta * grad.Mat
    
    # update Delta
    for (i in 1:p){
      for (j in 1:p){
        Asub <- A[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        Fnorm <- frobenius.norm(Asub)
        if(Fnorm <= lambda*Eta){
          Delta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, nrow=M, ncol=M)
        } else {
          Delta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- ((Fnorm - lambda*Eta) / Fnorm) * Asub
        }
      }
    }
    
    Delta.new <- Delta
    obj.new <- object_func(cov.X, cov.Y, Delta.new, lambda, M)
    
    if (abs((obj.new-obj.old)/(obj.old+1e-15)) < e_tol){
      converge.indicator <- TRUE
      Delta.hat <- 0.5 * (Delta.new + t(Delta.new))
      Delta.frob <- blockwise_Frob(Delta.hat, M)
      
      return (list(DelatMathat=Delta.hat, blockFrob=Delta.frob, converge=converge.indicator, 
                   num.iter=t))
      break
    }
    Delta.old <- Delta.new
  }
  
  Delta.hat <- 0.5 * (Delta.new + t(Delta.new))
  Delta.frob <- blockwise_Frob(Delta.hat, M)
  
  return (list(DelatMathat=Delta.hat, blockFrob=Delta.frob, converge=converge.indicator, 
               num.iter=t))
}