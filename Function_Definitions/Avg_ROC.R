Avg_ROC <- function(ROC_list, v.FPR){
  # Param:
  #  ROC_list: A list of ROC records
  #  v.FPR: the given FPR vector
  # Return:
  #  A vector of averaged TPR
  M <- length(ROC_list)
  n <- length(v.FPR)
  v.TPR <- numeric(n)
  
  for (i in 1:n){
    brk.pt <- v.FPR[i]
    if (brk.pt==0){
      v.TPR[i] <- 0
    } else if (brk.pt==1){
      v.TPR[i] <- 1
    } else {
      TPR.record <- numeric(M)
      for (j in 1:M){
        ROC_FPR <- ROC_list[[j]][, "FPR"]
        ROC_TPR <- ROC_list[[j]][, "TPR"]
        
        if (ROC_FPR[length(ROC_FPR)] < 1){
          ROC_FPR <- c(ROC_FPR, 1)
          ROC_TPR <- c(ROC_TPR, 1)
        } else if (ROC_TPR[length(ROC_TPR)] < 1){
          ROC_FPR <- c(ROC_FPR, 1)
          ROC_TPR <- c(ROC_TPR, 1)
        }
        
        index <- sum(brk.pt>ROC_FPR)
        if (index == 0){
          next
        }
        
        if (ROC_FPR[index+1]-ROC_FPR[index] > 0){
          slope <- (ROC_TPR[index+1]-ROC_TPR[index])/(ROC_FPR[index+1]-ROC_FPR[index])
          TPR.record[j] <- slope * (brk.pt-ROC_FPR[index]) + ROC_TPR[index]
        } else {
          TPR.record[j] <- (ROC_TPR[index]+ROC_TPR[index+1])/2
        }
        
        
      }
      v.TPR[i] <- mean(TPR.record)
    }
  }
  
  return(v.TPR)
}
