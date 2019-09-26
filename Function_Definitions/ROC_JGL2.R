library(matrixcalc)
ROC_JGL2 <- function(JGL.result, SuppDelta, p, M, thresh=0.02){
  # Para:
  #   JGL.result: The object returned by JGL function
  #   SuppDelta: Support Matrix for Delta
  #   p: Number of nodes
  #   M: Number of selected principle components
  # Retunrs:
  #   A vector of c(FPR, TPR)
  ThetaX <- JGL.result$est[[1]]
  ThetaY <- JGL.result$est[[2]]
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  block.norm <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- frobenius.norm(ThetaX[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] - 
                               ThetaY[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)])/(M^2)
      block.norm[i, j] <- diff
      if (diff>thresh){
        Edgeshat[i,j] <- 1
      }
    }
  }
  
  v.estimate <- c()
  v.true <- c()
  for (i in 1:(p - 1)){
    v.estimate <- c(v.estimate, Edgeshat[i, (i + 1):p])
    v.true <- c(v.true, SuppDelta[i, (i + 1):p])
  }
  
  n <- length(v.estimate)
  TP <- rep(0, n)
  TN <- rep(0, n)
  FP <- rep(0, n)
  FN <- rep(0, n)
  
  TP <- as.double(sum(v.estimate & v.true))
  TN <- as.double(sum(!v.estimate & !v.true))
  FP <- as.double(sum(v.estimate & !v.true))
  FN <- as.double(sum(!v.estimate & v.true))
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  
  return(list(ROC.point=c(FPR, TPR), block.norm=block.norm))
}