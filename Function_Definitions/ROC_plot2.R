ROC_plot2 <- function(cov.X, cov.Y, TrueMat, p, M, lambdamax=2, lambdamin=0, steplength=0.1){
  # Plot the ROC curve for estimated differential edges against lambda
  #
  # Args:
  #   cov.X: estimated covariance of Graph X
  #   cov.Y: estimated covariance of Graph Y
  #   TrueMat: matrix representing true differential edges
  #   p: the number of edges
  #   M: the number of principle components
  #   lambdamax: the maximum number of lambda
  #
  # Returns:
  #   A list containing the FPR and TPR record
  
  lambda.candidates <- seq(lambdamax, lambdamin, -steplength)
  TPR.record <- numeric(length(lambda.candidates))
  FPR.record <- numeric(length(lambda.candidates))
  
  k <- 0
  for (lambda in lambda.candidates){
    k <- k + 1
    g <- ProxAlg(cov.X, cov.Y, p, M, lambda)
    FrobMat <- g$blockFrob
    
    Edgeshat <- matrix(0, nrow=p, ncol=p)
    Edgeshat[which(FrobMat > 0)] <- 1
    
    v.estimate <- c()
    v.true <- c()
    for (i in 1:(p - 1)){
      v.estimate <- c(v.estimate, Edgeshat[i, (i + 1):p])
      v.true <- c(v.true, TrueMat[i, (i + 1):p])
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
    
    TPR.record[k] <- TPR
    FPR.record[k] <- FPR
  }
  
  return(list(TPR=TPR.record, FPR=FPR.record))
}
