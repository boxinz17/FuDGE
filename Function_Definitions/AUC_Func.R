AUC_Func <- function(v.FPR, v.TPR){
  v.FPR <- v.FPR[!is.na(v.FPR)]
  v.TPR <- v.TPR[!is.na(v.TPR)]
  
  if (v.FPR[1]!=0 | v.TPR[1]!=0){
    if (v.TPR[1]==1){
      v.FPR <- c(0, v.FPR)
      v.TPR <- c(1, v.TPR)
    } else {
      v.FPR <- c(0, v.FPR)
      v.TPR <- c(0, v.TPR)
    }
  }
  n <- length(v.FPR)
  S <- 0
  for (k in 1:(n-1)){
    # if (v.TPR[k] > v.TPR[k+1]){
    #   stop("Error: TPR decreases")
    # } else if (v.FPR[k] > v.FPR[k+1]){
    #   stop("Error: FPR decreases")
    # }
    
    # if (v.FPR[k] > v.FPR[k+1]){
    #   stop("Error: FPR decreases")
    # }
    
    h1 <- v.TPR[k]
    h2 <- v.TPR[k+1]
    w <- v.FPR[k+1] - v.FPR[k]
    S <- S + 0.5 * (h1+h2) * w
  }
  return(S)
}