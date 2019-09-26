F1_score <- function(SuppMat.hat, SuppMat.true){
  # Param:
  #  SuppMat.hat: Estimated Support Matrix
  #  SuppMat.true: True Support Matrix
  p <- dim(SuppMat.hat)[1]
  
  v.estimate <- c()
  v.true <- c()
  for (i in 1:(p - 1)){
    v.estimate <- c(v.estimate, SuppMat.hat[i, (i + 1):p])
    v.true <- c(v.true, SuppMat.true[i, (i + 1):p])
  }
  
  n <- length(v.estimate)
  TP <- rep(0, n)
  TN <- rep(0, n)
  FP <- rep(0, n)
  FN <- rep(0, n)
  
  TP <- sum(v.estimate & v.true)
  TN <- sum(!v.estimate & !v.true)
  FP <- sum(v.estimate & !v.true)
  FN <- sum(!v.estimate & v.true)
  
  recal <- TP / (TP + FN)
  Preci <- TP / (TP + FP)
  
  f1 <- 2*Preci*recal/(Preci+recal)
  return(f1)
}

