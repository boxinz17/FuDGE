AIC_JGL <- function(DataList, JGL.result){
  # Para:
  #   DataList: A list of nxp data matrices
  #   JGL.result: The object returned by JGL function
  # Returns:
  #   An AIC value
  K <- length(DataList)
  Theta.List <- JGL.result$theta
  AIC <- 0
  for (k in 1:K){
    Xk <- DataList[[k]]
    nk <- dim(Xk)[1]
    Xk.c <- scale(Xk, center=T, scale=F)
    Sk <- (t(Xk.c) %*% Xk.c) / (nk-1)
    Thetak <- Theta.List[[k]]
    AIC <- AIC + (nk*sum(diag(Sk%*%Thetak)) - nk*log(det(Thetak)) + 2*sum(Thetak>0))
  }
  return(AIC)
}