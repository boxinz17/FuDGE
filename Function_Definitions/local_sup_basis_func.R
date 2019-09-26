MyLocalFourier <- function(x, m, k){
  n <- length(x)
  y <- numeric(n)
  for (i in 1:n){
    if (x[i]>(k-1)/m & x[i]<k/m){
      y[i] <- cos((x[i]-(2*k-1)/(2*m))*2*m*pi)+1
    } else {
      y[i] <- 0
    }
  }
  return(y)
}

local_sup_basis_func <- function(u, m){
  # Param:
  #  u: the observation grid of length n
  #  m: number of basis
  # Return:
  #  A nxm matrix of function values
  n <- length(u)
  ObservMat <- matrix(0, nrow=n, ncol=m)
  for (i in 1:n){
    for (j in 1:m){
      ObservMat[i, j] <- MyLocalFourier(u[i], m, j)
    }
  }
  return(ObservMat)
}
