ConstrainedL2 <- function(Amat, bvec, lambda){
  # Function to solve the constrained L2 minimization problem
  #
  # Args:
  #   Amat: coefficient matrix before beta
  #   bvec: Intercept vector
  #   lambda: threshold of constraint
  #
  # Returns:
  #   A solution vector
  if (is.vector(Amat)){
    Amat <- as.matrix(Amat)
  }
  constrained.dim <- dim(Amat)[1]
  x.dim <- dim(Amat)[2]
  Dmat <- diag(x.dim)
  dvec <- rep(0, x.dim)
  Amat <- t(rbind(Amat, -Amat))
  cvec <- rep(lambda, constrained.dim)
  bvec <- c(bvec - cvec, - (cvec +  bvec))
  sol  <- solve.QP(Dmat, dvec, Amat, bvec, meq=0)
  
  return(sol$solution)
}

# test
#Amat <- diag(1, 4)
#bvec <- rep(1, 4)
#lambda <- 0.5
#ConstrainedL2(Amat, bvec, lambda)