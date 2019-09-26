ProxAlg <- function(cov.X, cov.Y, p, M, lambda, Eta=0.1, n.iteration=1000){
  # Optimzation function to solve the optimization problem in DFGM
  # Use Proximal Algorithm
  # 
  # Args:
  #   cov.X: the estimated covariance matrix of Graph X
  #   cov.Y: the estimated covariance matrix of Graph Y
  #   p: the number of vertices
  #   M: the number of principle components
  #   lambda: hyperparameter
  #   Eta: hyperparamter
  #   n.iteration: maximum times of iteration
  #
  # Returns:
  #   A list:
  #     DelatMathat: The estimated differential matrix between ThetaX and ThetaY
  #     blockFrob: A matrix of Frobenius norm of each block matrix
  delta.record <- array(0, c(n.iteration, (M * p), (M * p)))
  delta.record[1, , ] <- matrix(0, nrow=(M * p), ncol=(M * p))
  grad.L.record <- array(0, c(n.iteration, (M * p), (M * p)))
  
  for (t in 1:(n.iteration - 1)){
    delta.old <- as.matrix(delta.record[t, , ])
    
    grad.L <- t(0.5 * (cov.Y %*% delta.old %*% cov.X) + 
                  0.5 * (cov.X %*% delta.old %*% cov.Y) - (cov.Y - cov.X))
    grad.L.record[t, , ] <- grad.L
    A <- delta.old - Eta * grad.L
    
    delta.new <- delta.old
    
    for (i in 1:p){
      for (j in 1:p){
        Asub <- A[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)]
        
        Fnorm <- frobenius.norm(Asub)
        if(Fnorm <= lambda){
          delta.new[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 
            matrix(0, nrow=M, ncol=M)
        } else{
          delta.new[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 
            ((Fnorm - lambda) / Fnorm) * Asub
        }
      }
    }
    
    ## Symmertize the delta.new
    for (i in 1:(p *M)){
      for(j in 1:(p * M)){
        if(abs(delta.new[i, j]) < abs(delta.new[j, i])){
          delta.new[j, i] <- delta.new[i, j]
        } else {
          delta.new[i, j] <- delta.new[j, i]
        }
      }
    }
    
    delta.record[t + 1, , ] <- delta.new
    
    ## check overflow
    if (frobenius.norm(delta.new) / (M * p) > 1e+12){
      break
    }
    
    ## check convergence
    diff <- delta.record[t+1, , ] - delta.record[t, , ]
    if(frobenius.norm(diff) / (M * p) < 
       1e-3 * (frobenius.norm(delta.record[t, , ]) / (M * p))){
      break
    }
  }
  Delta <- delta.record[t + 1, , ]
  delta.record <- delta.record[1:(t+1), , ]
  
  Delta.frob <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      Delta.frob[i, j] <- frobenius.norm(Delta[(((i - 1) * M) + 1):(i * M), 
                                               (((j - 1) * M) + 1):(j * M)])
    }
  }
  
  return (list(DelatMathat=Delta, blockFrob=Delta.frob, DeltaRecord=delta.record))
}
