AUC_Func2 <- function(v.FPR, v.TPR){
  n <- length(v.FPR)
  S <- 0
  for (k in 1:(n-1)){
    h1 <- v.TPR[k]
    h2 <- v.TPR[k+1]
    w <- v.FPR[k+1] - v.FPR[k]
    S <- S + max(h1, h2) * w
  }
  return(S)
}