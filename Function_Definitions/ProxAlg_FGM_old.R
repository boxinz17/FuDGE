ProxAlg_FGM <- function(S, p, M, gam, Eta="Auto", n.iteration=2000){
  # Optimzation function to solve the optimization problem in FGM
  # Use Proximal Algorithm
  # 
  # Args:
  #   S: the estimated covariance matrix of Graph
  #   p: the number of vertices
  #   M: the number of principle components
  #   gam: hyperparameter.
  #   Eta: hyperparamter. If Eta="Auto", then Eta would be set as inverse of the 
  #        product of the maximum singular value of cov.X and cov.Y
  #   n.iteration: maximum times of iteration
  #
  # Returns:
  #   A list:
  #     ThetaMathat: The estimated differential matrix between ThetaX and ThetaY
  #     blockFrob: A matrix of Frobenius norm of each block matrix
  #     Support: Support Matrix
  if(Eta == "Auto"){
    Eta <- 0.1 # this needs more careful design
  }
  
  Theta.record <- array(0, c(n.iteration, (M * p), (M * p)))
  Theta.record[1, , ] <- diag(1, nrow=(M * p), ncol=(M * p))
  grad.L.record <- array(0, c(n.iteration, (M * p), (M * p)))
  
  for (t in 1:(n.iteration - 1)){
    Theta.old <- as.matrix(Theta.record[t, , ])
    
    grad.L <- solve(Theta.old) - S
    grad.L.record[t, , ] <- grad.L
    A <- Theta.old - Eta * grad.L
    
    Theta.new <- Theta.old
    
    for (i in 1:p){
      for (j in 1:p){
        Asub <- A[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)]
        
        if (i == j){
          Theta.new[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- Asub
          next
        }
        
        Fnorm <- frobenius.norm(Asub)
        if(Fnorm <= gam*Eta){
          Theta.new[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 
            matrix(0, nrow=M, ncol=M)
        } else{
          Theta.new[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 
            ((Fnorm - gam*Eta) / Fnorm) * Asub
        }
      }
    }
    
    ## Symmertize the Theta.new
    for (i in 1:(p *M)){
      for(j in 1:(p * M)){
        if(abs(Theta.new[i, j]) < abs(Theta.new[j, i])){
          Theta.new[j, i] <- Theta.new[i, j]
        } else {
          Theta.new[i, j] <- Theta.new[j, i]
        }
      }
    }
    
    Theta.record[t + 1, , ] <- Theta.new
    
    ## check overflow
    if (frobenius.norm(Theta.new) / (M * p) > 1e+12){
      break
    }
    
    converge.indicator <- FALSE
    ## check convergence
    diff <- Theta.record[t+1, , ] - Theta.record[t, , ]
    if(frobenius.norm(diff) / (M * p) < 
       1e-3 * (frobenius.norm(Theta.record[t, , ]) / (M * p))){
      converge.indicator <- TRUE
      break
    }
  }
  Theta <- Theta.record[t + 1, , ]
  Theta.record <- Theta.record[1:(t+1), , ]
  
  Theta.frob <- matrix(0, nrow=p, ncol=p)
  SupMat <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      Theta.frob[i, j] <- frobenius.norm(Theta[(((i - 1) * M) + 1):(i * M), 
                                               (((j - 1) * M) + 1):(j * M)])
      if (Theta.frob[i, j] > 0){
        SupMat[i, j] <- 1
      }
    }
  }
  
  return (list(ThetaMathat=Theta, blockFrob=Theta.frob, Support=SupMat, 
               converge=converge.indicator, num.iter=t))
}
