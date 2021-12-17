library(matrixcalc)

sub_obj <- function(W, B, lambda1, lambda2, rho){
  # Para:
  #   W: A list of length 2 of W matrices
  #   B: A list of length 2 of B matrices
  #   lambda1, lambda2, rho: parameters passed by main function
  # Return:
  #   A scalar function value
  part1 <- 0.5 * (sum((W[[1]]-B[[1]])^2) + sum((W[[2]]-B[[2]])^2))
  part2 <- (lambda1/rho)*(frobenius.norm(W[[1]])+frobenius.norm(W[[2]]))
  part3 <- (lambda2/rho)*(sum(abs(W[[1]]-W[[2]])))
  return(part1+part2+part3)
}

# # test
# W1 <- matrix(c(1,0,0,1), nrow=2)
# W2 <- matrix(c(1,0,0,1), nrow=2)
# W <- list(W1, W2)
# 
# B1 <- matrix(c(0.5,0,0,0.5), nrow=2)
# B2 <- matrix(c(0.5,0,0,0.5), nrow=2)
# B <- list(B1, B2)
# 
# lambda1 <- 1
# lambda2 <- 1
# rho <- 1
# 
# sub_obj(W, B, lambda1, lambda2, rho)

sub_ADMM <- function(B, lambda1, lambda2, rho, tol.abs=1e-3, tol.rel=1e-3){
  # Para:
  #   B: A list of length 2 of B matrices
  #   lambda1, lambda2, rho: parameters passed by main function
  #   tol: the tolerance used to check convergence
  # Return:
  #   A list of length 2 of minimizer matrices
  B1 <- B[[1]]
  B2 <- B[[2]]
  
  M <- dim(B1)[1]
  
  W1_old <- diag(M)
  W2_old <- diag(M)
  R1_old <- matrix(0, nrow=M, ncol=M)
  R2_old <- matrix(0, nrow=M, ncol=M)
  V1_old <- matrix(0, nrow=M, ncol=M)
  V2_old <- matrix(0, nrow=M, ncol=M)
  
  rho.prime <- 1
  
  ind.primal <- FALSE
  ind.dual <- FALSE
  
  obj.value.record <- c()
  num.iter <- 0
  
  while((!ind.primal)|(!ind.dual)){
    num.iter <- num.iter+1
    
    # (i) update W
    C1 <- (1/(1+rho.prime))*(B1+rho.prime*(R1_old-V1_old))
    fnorm <- frobenius.norm(C1)
    if (fnorm<=lambda1/(rho*(1+rho.prime))){
      W1_new <- matrix(0, nrow=M, ncol=M)
    } else {
      s <- (fnorm-lambda1/(rho*(1+rho.prime)))/fnorm
      W1_new <- s*C1
    }
    
    C2 <- (1/(1+rho.prime))*(B2+rho.prime*(R2_old-V2_old))
    fnorm <- frobenius.norm(C2)
    if (fnorm<=lambda1/(rho*(1+rho.prime))){
      W2_new <- matrix(0, nrow=M, ncol=M)
    } else {
      s <- (fnorm-lambda1/(rho*(1+rho.prime)))/fnorm
      W2_new <- s*C2
    }
    
    # (ii) update R
    D1 <- W1_new + V1_old
    D2 <- W2_new + V2_old
    
    R1_new <- matrix(0, nrow=M, ncol=M)
    R2_new <- matrix(0, nrow=M, ncol=M)
    for (a in 1:M){
      for (b in 1:M){
        if (D1[a,b]>D2[a,b]+(2*lambda2)/(rho*rho.prime)){
          R1_new[a,b] <- D1[a,b]-lambda2/(rho*rho.prime)
          R2_new[a,b] <- D2[a,b]+lambda2/(rho*rho.prime)
        } else if (D1[a,b]<D2[a,b]-(2*lambda2)/(rho*rho.prime)){
          R1_new[a,b] <- D1[a,b]+lambda2/(rho*rho.prime)
          R2_new[a,b] <- D2[a,b]-lambda2/(rho*rho.prime)
        } else {
          R1_new[a,b] <- 0.5*(D1[a,b]+D2[a,b])
          R2_new[a,b] <- 0.5*(D1[a,b]+D2[a,b])
        }
      }
    }
    
    # (iii) update V
    V1_new <- V1_old + (W1_new - R1_new)
    V2_new <- V2_old + (W2_new - R2_new)
    
    # check convergence
    #obj.val.old <- sub_obj(list(W1_old, W2_old), B, lambda1, lambda2, rho)
    obj.val.new <- sub_obj(list(W1_new, W2_new), B, lambda1, lambda2, rho)
    
    obj.value.record <- c(obj.value.record, obj.val.new)

    #diff <- abs(obj.val.new-obj.val.old)/(abs(obj.val.old)+1e-18)
    
    res.primal <- sqrt(frobenius.norm(W1_new-R1_new)^2+frobenius.norm(W2_new-R2_new)^2)
    
    res.dual <- rho.prime*(sqrt(frobenius.norm(R1_new-R1_old)^2+frobenius.norm(R2_new-R2_old)^2))
    
    tol.primal <- M*tol.abs + tol.rel*max(sqrt(frobenius.norm(W1_new)^2+frobenius.norm(W2_new)^2), 
                                          sqrt(frobenius.norm(R1_new)^2+frobenius.norm(R2_new)^2))
    
    tol.dual <- M*tol.abs+tol.rel*sqrt(frobenius.norm(V1_new)^2+frobenius.norm(V2_new)^2)
    
    if (res.primal<=tol.primal){
      ind.primal <- TRUE
    }
    
    if (res.dual<=tol.dual){
      ind.dual <- TRUE
    }
    
    # Finally, update olds
    W1_old <- W1_new
    W2_old <- W2_new
    R1_old <- R1_new
    R2_old <- R2_new
    V1_old <- V1_new
    V2_old <- V2_new
  }
  
  result <- list(est=list(W1_new, W2_new), est2=list(R1_old, R2_old), est3=list(V1_old, V2_old),
                 value=obj.value.record, num.iter=num.iter)
  return(result)
}

# # test
# B1 <- matrix(c(0.5,0,0,0.5), nrow=2)
# B2 <- matrix(c(1.5,0,0,1), nrow=2)
# B <- list(B1, B2)
# 
# lambda1 <- 0.1
# lambda2 <- 100
# rho <- 1
# 
# res <- sub_ADMM(B, lambda1, lambda2, rho)


main_obj <- function(Theta, S, n.sam, p, M, lambda1, lambda2){
  # para:
  #   Theta: A list of length 2 of two precision matrices
  #   S: A list of length 2 of two sample matrices
  #   n.sam: A list of length 2 of sample sizes
  #   p: the number of vertices
  #   M: the number of principle components used
  #   lambda1: the parameter controlling individual sparsity
  #   lambda2: the parameter controlling similarity
  # Return:
  #   A scalar function value
  Theta1 <- Theta[[1]]
  Theta2 <- Theta[[2]]
  S1 <- S[[1]]
  S2 <- S[[2]]
  n1 <- n.sam[[1]]
  n2 <- n.sam[[2]]
  
  # The loss part
  Los <- - n1*(log(det(Theta1))-sum(diag(S1%*%Theta1))) - 
    n2*(log(det(Theta2))-sum(diag(S2%*%Theta2)))
  
  # The penalty part
  Pen <- 0
  for (i in 1:p){
    for (j in 1:p){
      if (i != j){
        sub.mat1 <- Theta1[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        sub.mat2 <- Theta2[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        Pen <- Pen + lambda1*(frobenius.norm(sub.mat1)+frobenius.norm(sub.mat2))
        Pen <- Pen + lambda2*sum(abs(sub.mat1-sub.mat2))
      } else {
        sub.mat1 <- Theta1[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        sub.mat2 <- Theta2[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        Pen <- Pen + lambda2*sum(abs(sub.mat1-sub.mat2))
      }
    }
  }
  
  # Return result
  return(Los+Pen)
}

# # test
# Theta1 <- diag(4)
# Theta2 <- diag(4)
# Theta <- list(Theta1, Theta2)
# 
# S1 <- diag(rep(2, 4))
# S2 <- diag(rep(2, 4))
# S <- list(S1, S2)
# 
# n.sam <- c(100, 100)
# p <- 2
# M <- 2
# lambda1 <- 0.1
# lambda2 <- 0.1
# 
# main_obj(Theta, S, n.sam, p, M, lambda1, lambda2)

FFGL2_ADMM <- function(S, n.sam, p, M, lambda1, lambda2, tol.abs=1e-2, tol.rel=1e-3, mu=10, 
                      tau.incr=1.1, tau.dec=1.1, pre.res=NULL){
  # para:
  #   S: A list of length 2 of two sample matrices
  #   n.sam: A list of length 2 of sample sizes
  #   p: the number of vertices
  #   M: the number of principle components used
  #   lambda1: the parameter controlling individual sparsity
  #   lambda2: the parameter controlling similarity
  #   tol: the tolerance used to check convergence
  # Return:
  #   A list of results
  S1 <- S[[1]]
  S2 <- S[[2]]
  n1 <- n.sam[[1]]
  n2 <- n.sam[[2]]
  
  if (is.null(pre.res)){
    Theta1_old <- diag(p*M)
    Theta2_old <- diag(p*M)
    Z1_old <- matrix(0, nrow=(p*M), ncol=(p*M))
    Z2_old <- matrix(0, nrow=(p*M), ncol=(p*M))
    U1_old <- matrix(0, nrow=(p*M), ncol=(p*M))
    U2_old <- matrix(0, nrow=(p*M), ncol=(p*M))
    rho <- 1
  } else {
    Theta1_old <- pre.res$est[[1]]
    Theta2_old <- pre.res$est[[2]]
    Z1_old <- pre.res$est2[[1]]
    Z2_old <- pre.res$est2[[2]]
    U1_old <- pre.res$est3[[1]]
    U2_old <- pre.res$est3[[2]]
    rho <- pre.res$rho
  }
  
  ind.primal <- FALSE
  ind.dual <- FALSE
  
  obj.value.record <- c()
  num.iter <- 0
  
  res.pri.v <- c()
  res.dual.v <- c()
  
  convg.indc <- TRUE
  
  stop.update.count <- 0
  
  while((!ind.primal)|(!ind.dual)){
    num.iter <- num.iter+1
    
    # (i) update Theta
    mat.temp1 <- S1 - (rho/n1)*Z1_old + (rho/n1)*U1_old
    mat.temp2 <- S2 - (rho/n2)*Z2_old + (rho/n2)*U2_old
    
    eigen.not.real <- TRUE
    count.temp <- 0
    while(eigen.not.real){
      count.temp <- count.temp+1
      
      eigen.res1 <- eigen(mat.temp1)
      eigen.res2 <- eigen(mat.temp2)
      
      if (is.numeric(sum(eigen.res1$vectors)) & is.numeric(sum(eigen.res2$vectors))){
        eigen.not.real <- FALSE
      } else {
        rho <- rho/2
        U1_old <- U1_old*2
        U2_old <- U2_old*2
        
        mat.temp1 <- S1 - (rho/n1)*Z1_old + (rho/n1)*U1_old
        mat.temp2 <- S2 - (rho/n2)*Z2_old + (rho/n2)*U2_old
      }
      
      if (count.temp>10){
        break
      }
    }
    
    D1 <- diag(eigen.res1$values)
    V1 <- eigen.res1$vectors
    
    D.tilde1 <- D1
    for (j in 1:dim(D1)[1]){
      D.tilde1[j,j] <- (n1/(2*rho))*(-D1[j,j]+sqrt(D1[j,j]^2+4*rho/n1))
    }
    
    Theta1_new <- V1 %*% D.tilde1 %*% t(V1)
    
    D2 <- diag(eigen.res2$values)
    V2 <- eigen.res2$vectors
    
    D.tilde2 <- D2
    for (j in 1:dim(D2)[1]){
      D.tilde2[j,j] <- (n2/(2*rho))*(-D2[j,j]+sqrt(D2[j,j]^2+4*rho/n2))
    }
    
    Theta2_new <- V2 %*% D.tilde2 %*% t(V2)
    
    # (ii) update Z
    A1 <- Theta1_new + U1_old
    A2 <- Theta2_new + U2_old
    
    Z1_new <- matrix(0, nrow=(p*M), ncol=(p*M))
    Z2_new <- matrix(0, nrow=(p*M), ncol=(p*M))
    
    for (i in 1:p){
      for (j in 1:p){
        B1 <- A1[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        B2 <- A2[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        
        if (i == j){
          for (a in 1:M){
            for (b in 1:M){
              if (B1[a,b]>B2[a,b]+(2*lambda2)/rho) {
                Z1_new[((i-1)*M+a), ((j-1)*M+b)] <- B1[a,b]-lambda2/rho
                Z2_new[((i-1)*M+a), ((j-1)*M+b)] <- B1[a,b]+lambda2/rho
              } else if (B1[a,b]<B2[a,b]-(2*lambda2)/rho) {
                Z1_new[((i-1)*M+a), ((j-1)*M+b)] <- B1[a,b]+lambda2/rho
                Z2_new[((i-1)*M+a), ((j-1)*M+b)] <- B1[a,b]-lambda2/rho
              } else {
                Z1_new[((i-1)*M+a), ((j-1)*M+b)] <- 0.5*(B1[a,b]+B2[a,b])
                Z2_new[((i-1)*M+a), ((j-1)*M+b)] <- 0.5*(B1[a,b]+B2[a,b])
              }
            }
          }
        } else {
          B <- list(B1, B2)
          res.temp <- sub_ADMM(B, lambda1, lambda2, rho)
          
          Z1_new[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- res.temp$est[[1]]
          Z2_new[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- res.temp$est[[2]]
        }
      }
    }
    
    # (iii) update U
    U1_new <- U1_old + (Theta1_new - Z1_new)
    U2_new <- U2_old + (Theta2_new - Z2_new)
    
    # check convergence
    obj.val.old <- main_obj(list(Theta1_old, Theta2_old), S, n.sam, p, M, lambda1, lambda2)
    obj.val.new <- main_obj(list(Theta1_new, Theta2_new), S, n.sam, p, M, lambda1, lambda2)
    
    res.primal <- sqrt(frobenius.norm(Theta1_new-Z1_new)^2+frobenius.norm(Theta2_new-Z2_new)^2)
    res.dual <- rho*(sqrt(frobenius.norm(Z1_new-Z1_old)^2+frobenius.norm(Z2_new-Z2_old)^2))
    tol.primal <- (M*p)*tol.abs + tol.rel*max(sqrt(frobenius.norm(Theta1_new)^2+frobenius.norm(Theta2_new)^2), 
                                              sqrt(frobenius.norm(Z1_new)^2+frobenius.norm(Z2_new)^2))
    tol.dual <- (M*p)*tol.abs+tol.rel*sqrt(frobenius.norm(U1_new)^2+frobenius.norm(U2_new)^2)
    
    if (res.primal<=tol.primal){
      ind.primal <- TRUE
    }
    
    if (res.dual<=tol.dual){
      ind.dual <- TRUE
    }
    
    # adjust rho
    if (res.primal>mu*res.dual){
      rho <- tau.incr*rho
      U1_new <- U1_new/tau.incr
      U2_new <- U2_new/tau.incr
    } else if (res.dual>mu*res.primal){
      rho <- rho/tau.dec
      U1_new <- U1_new*tau.dec
      U2_new <- U2_new*tau.dec
    } else {
      rho <- rho
    }
    
    # Finally, update olds
    update.ind <- TRUE
    if (num.iter>1){
      if (5*res.primal.old<res.primal | 5*res.dual.old<res.dual){
        update.ind <- FALSE
      }
    }
    
    if (update.ind){
      Theta1_old <- Theta1_new
      Theta2_old <- Theta2_new
      Z1_old <- Z1_new
      Z2_old <- Z2_new
      U1_old <- U1_new
      U2_old <- U2_new
      
      res.primal.old <- res.primal
      res.dual.old <- res.dual
      
      obj.value.record <- c(obj.value.record, obj.val.new)
      res.pri.v <- c(res.pri.v, res.primal)
      res.dual.v <- c(res.dual.v, res.dual)
      
      stop.update.count <- 0
      
      print(num.iter)
      print(paste("tol.primal=", tol.primal, "res.primal=", res.primal))
      print(paste("tol.dual=", tol.dual, "res.dual=", res.dual))
    } else {
      rho <- rho/tau.dec
      U1_old <- U1_old*tau.dec
      U2_old <- U2_old*tau.dec
      
      obj.value.record <- c(obj.value.record, obj.val.old)
      res.pri.v <- c(res.pri.v, res.primal.old)
      res.dual.v <- c(res.dual.v, res.dual.old)
      
      stop.update.count <- stop.update.count + 1
      
      print(num.iter)
      print(paste("tol.primal=", tol.primal, "res.primal=", res.primal.old))
      print(paste("tol.dual=", tol.dual, "res.dual=", res.dual.old))
    }
    
    print(paste("if update:", update.ind))
    
    # set up maximum iteration numbers
    if (num.iter>200){
      convg.indc <- FALSE
      break
    }
    
    if (stop.update.count>5){
      break
    }
  }
  
  result <- list(est=list(Theta1_new, Theta2_new), est2=list(Z1_new, Z2_new), 
                 est3=list(U1_new, U2_new), res.pri=res.pri.v, res.dual=res.dual.v, 
                 value=obj.value.record, num.iter=num.iter, rho=rho, convergence=convg.indc)
  return(result)
}

# # test
# n.sam <- c(100, 100)
# p <- 2
# M <- 2
# lambda1 <- 0.1
# lambda2 <- 1000
# 
# # S1 <- diag(rep(1, 4))
# # S2 <- diag(rep(1, 4))
# # S <- list(S1, S2)
# 
# set.seed(111)
# S1 <- matrix(rnorm((p*M)^2), nrow=(p*M), ncol=(p*M))
# S1 <- 0.5*(S1%*%t(S1))
# S2 <- matrix(rnorm((p*M)^2), nrow=(p*M), ncol=(p*M))
# S2 <- 0.5*(S2%*%t(S2))
# S <- list(S1, S2)
# 
# res <- FFGL2_ADMM(S, n.sam, p, M, lambda1, lambda2)

#################################################################
# S <- cov.list
# n.sam<-c(n,n)
# p<-p
# M<-M.selected
# lambda1<-lambda1.choose
# lambda2<-30
# tol.abs=1e-3
# tol.rel=1e-3
# mu=10 
# tau.incr=2
# tau.dec=2
