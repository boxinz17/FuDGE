DFGMOptimFun2 <- function(cov.X, cov.Y, lambda, p, M, alpha=T){
  # Optimzation function to solve the optimization problem in DFGM, method 2
  # 
  # Args:
  #   cov.X: the estimated covariance matrix of Graph X
  #   cov.Y: the estimated covariance matrix of Graph Y
  #   lambda: hyperparamter
  #   p: the number of vertices
  #   M: the number of principle components
  #   alpha: the step length of gradient descent, either a number of TRUE.
  #          If True, then we will set automatic step length.
  #
  # Returns:
  #   A list:
  #     DelatMathat: The estimated differential matrix between ThetaX and ThetaY
  #     blockFrob: A matrix of Frobenius norm of each block matrix
  n.iteration <- 1000
  delta.record <- array(0, c(n.iteration, (M * p), (M * p)))
  delta.record[1, , ] <- matrix(0, nrow=(M * p), ncol=(M * p))
  grad.L.record <- array(0, c(n.iteration, (M * p), (M * p)))
  grad.P.record <- array(0, c(n.iteration, (M * p), (M * p)))
  
  if (identical(alpha, T)){
    alpha <- min(1e-3 /(lambda + 1e-10), 1e-2)
  }
  
  subject.record <- numeric(n.iteration)
  
  for (t in 1:(n.iteration - 1)){
    delta.old <- as.matrix(delta.record[t, , ])
    
    grad.L <- t(0.5 * (cov.Y %*% delta.old %*% cov.X) + 
                  0.5 * (cov.X %*% delta.old %*% cov.Y) - (cov.Y - cov.X))
    grad.L.record[t, , ] <- grad.L
    
    grad.P <- matrix(0, nrow=(M * p), ncol=(M * p))
    for (i in 1:p){
      for (j in 1:p){
        submat <- delta.old[((i - 1) * M + 1):(i * M), 
                            ((j - 1) * M + 1):(j * M)] 
      }
      grad.P[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)]  <- 
        submat / (frobenius.norm(submat) + 1e-8)  # prevent denominator to be zero
    }
    grad.P.record[t, , ] <- grad.P
    
    part2.temp <- 0
    for (i in 1:p){
      for (j in 1:p){
        submat <- delta.old[((i - 1) * M + 1):(i * M), 
                            ((j - 1) * M + 1):(j * M)]
        part2.temp <- part2.temp + frobenius.norm(submat)
      }
    }
    subject <- sum(diag(0.5 * (cov.X %*% delta.old %*% cov.Y %*% delta.old) - 
                          delta.old %*% (cov.Y - cov.X))) + lambda * part2.temp
    subject.record[t] <- subject
    
    if (abs(subject.record[t]) > 1e8){
      break
    } else if ((t > 1)){
      if ((subject.record[t] - subject.record[t - 1]) < 1e-8){
        break
      }
    }
    
    delta.record[t + 1, , ] <- delta.old - alpha * (grad.L + lambda * grad.P)
  }
  Delta <- delta.record[t, , ]
  
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
      Delta.frob[i, j] <- frobenius.norm(Delta[(((i - 1) * M) + 1):(i * M), 
                                               (((j - 1) * M) + 1):(j * M)])
    }
  }
  
  return (list(DelatMathat=Delta, blockFrob=Delta.frob))
}