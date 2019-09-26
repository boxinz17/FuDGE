# This is an updated version of old ProxAlg.R
ProxAlg_DFGM <- function(cov.X, cov.Y, p, M, lambda, Eta="Auto", n.iteration=2000){
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
        if(Fnorm <= lambda*Eta){
          delta.new[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 
            matrix(0, nrow=M, ncol=M)
        } else{
          delta.new[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 
            ((Fnorm - lambda*Eta) / Fnorm) * Asub
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
    converge.indicator <- FALSE
    diff <- delta.record[t+1, , ] - delta.record[t, , ]
    if(frobenius.norm(diff) / (M * p) < 
       1e-3 * (frobenius.norm(delta.record[t, , ]) / (M * p))){
      converge.indicator <- TRUE
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
  
  return (list(DelatMathat=Delta, blockFrob=Delta.frob, DeltaRecord=delta.record, 
               converge=converge.indicator, num.iter=t))
}
