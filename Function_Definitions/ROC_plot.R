ROC_plot <- function(FrobMat, TrueMat){
  # Plot the ROC curve for estimated differential edges
  # 
  # Args:
  #   FrobMat: The blockwise Frobenius norm matrix of estimated delat matrix
  #   TrueMat: The matrix showing the true differential graph edges
  #   
  # Returns:
  #   A list:
  #     TPR: The vector of True Postive Ratios
  #     FPR: The vector of False Postive Ratios
  p <- dim(FrobMat)[1]
  low_bound <- min(FrobMat)
  up_bound <- max(FrobMat)
  
  TPR.record <- numeric(length(seq(from=(up_bound + 0.01), to=(low_bound - 0.01), length.out=300)))
  FPR.record <- numeric(length(seq(from=(up_bound + 0.01), to=(low_bound - 0.01), length.out=300)))
  i <- 0
  for (t in seq(from=(up_bound + 0.01), to=(low_bound - 0.01), length.out=300)){
    i <- i + 1
    Edgeshat <- matrix(0, nrow=p, ncol=p)
    Edgeshat[which(FrobMat > t)] <- 1
    
    v.estimate <- as.vector(Edgeshat)
    v.true <- as.vector(TrueMat)
    
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
    
    TPR.record[i] <- TPR
    FPR.record[i] <- FPR
  }
  return(list(TPR=TPR.record, FPR=FPR.record))
}

# test
#FrobMat <- g$blockFrob
#TrueMat <- matrix(0, nrow=p, ncol=p)
#for (i in 1:(p - 1)){
#  TrueMat[i, i + 1] <- 1
#  TrueMat[i + 1, i] <- 1
#}

#g3 <- ROC_plot(FrobMat, TrueMat)
#plot(g3$FPR, g3$TPR, xlab="False Postive Rate", ylab="True Postive Rate", type="l", 
#     main=paste("ROC curve when p=", p))
