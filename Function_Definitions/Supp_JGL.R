library(matrixcalc)
Supp_JGL <- function(JGL.result, p, M){
  # Para:
  #   JGL.result: The object returned by JGL function
  #   p: Number of nodes
  #   M: Number of selected principle components
  # Returns:
  #   A list of support for graph X & Y
  Supp.X <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      block.norm <- frobenius.norm(JGL.result$theta[[1]][((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)])
      if (block.norm>0){
        Supp.X[i, j] <- 1
      }
    }
  }
  
  Supp.Y <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      block.norm <- frobenius.norm(JGL.result$theta[[2]][((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)])
      if (block.norm>0){
        Supp.Y[i, j] <- 1
      }
    }
  }
  
  return(list(X=Supp.X, Y=Supp.Y))
}