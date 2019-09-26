## Matrix organization function definition
MatOrg <- function(A, B, p, M){
  # Function to convert two matrices to the format we need
  #
  # Args:
  #   A: The Theta hat X matrix
  #   B: The Theta hat Y matrix
  #   p: The number of vertices
  #   M: The number of principle components
  #
  # Returns:
  #   A list containing the matrix list and vector we need
  Aarr <- array(0, c(p, p, ((p ^ 2) * (M ^ 2)), (M ^ 2)))
  for (i in 1:p){
    for (j in 1:p){
      for (k in 1:p){
        for (l in 1:p){
          submatrix1 <- A[((k - 1) * M + 1):(k * M), ((i - 1) * M + 1):(i * M)]
          submatrix2 <- B[((j - 1) * M + 1):(j * M), ((l - 1) * M + 1):(l * M)]
          submatrix <- kronecker(submatrix2, submatrix1)
          Aarr[i, j, (((k - 1) * p + l - 1) * (M ^ 2) + 1):
                 (((k - 1) * p + l) * (M ^ 2)), ] <- submatrix
        }
      }
    }
  }
  
  b <- numeric((p ^ 2) * (M ^ 2))
  for (k in 1:p){
    for (l in 1:p){
      submatrix1 <- A[((k - 1) * M + 1):(k * M), ((l - 1) * M + 1):(l * M)]
      submatrix2 <- B[((k - 1) * M + 1):(k * M), ((l - 1) * M + 1):(l * M)]
      subvector <- as.vector(submatrix1 - submatrix2)
      b[(((k - 1) * p + l - 1) * (M ^ 2) + 1):(((k - 1) * p + l) * (M ^ 2))] <- subvector
    }
  }
  
  return (list(matrixList=Aarr, vector=b))
}

### test
#A <-matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1), nrow=4, byrow=F)
#B <-matrix(c(2, 2, 0, 0, 2, 2, 0, 0, 0, 0, 2, 2, 0, 0, 2, 2), nrow=4, byrow=F)
#p <- 2
#M <- 2
#g <- MatOrg(A, B, p, M)
#g$matrixList[1, 1, , ]
#g$matrixList[1, 2, , ]
#g$matrixList[2, 1, , ]
#g$matrixList[2, 2, , ]
#g$vector

