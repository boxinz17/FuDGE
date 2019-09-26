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