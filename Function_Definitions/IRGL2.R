source("E:/学习工作/硕士/graduate research/code/Differential Functional Graphical Model/experiment_1/Function_Definitions/ConstrainedL2.R")

IRGL2 <- function(Aarr, bvec, lambda, p, M, n.iteration=100, stop.criterion=0.05){
  # Use Iterative approach to solve group L2 norm with max element constraint
  # Args:
  #   Aarr: An array containing the matrices we need
  #   bvec: The vector we need
  #   lambda: threshold for max element constraint
  #   p: The number of vertices
  #   M: The number of principle components
  #   n.iteration: The number of maxmium iteration times
  #   stop.criterion: The critertion to stop the iteration, should be a postive number
  #                   indicating the ratio of lambda
  #
  # Returns:
  #   A list with: 
  #     betahat: the estimation of beta
  #     betarecord: the history of estimations of beta
  #     lossrecord: the history of loss function values
  #     numiter: the number of iterations completed
  Beta <- array(0, c(n.iteration, p, p, M ^ 2))
  beta.record <- array(0, c(n.iteration, (p^2) * (M ^ 2)))
  loss.record <- rep(Inf, n.iteration)
  
  # Intialize Beta 1
  for (i in 1:p){
    for (j in 1:p){
      if ((i == 1) & (j == 1)){
        A.temp <- Aarr[i, j, , ]
      } else {
        A.temp <- cbind(A.temp, Aarr[i, j, , ])
      }
    }
  }
  if (is.singular.matrix(A.temp)){
    beta.temp <- bvec %*% t(ginv(A.temp))
  } else {
    beta.temp <- solve(A.temp, bvec)
  }
  
  beta.record[1, ] <- beta.temp
  
  for(i in 1:p){
    for (j in 1:p){
      Beta[1, i, j, ] <- beta.temp[(((i - 1) * p + j - 1) * (M ^ 2) + 1):(((i - 1) * p + j) * (M ^ 2))]
    }
  }
  
  for(i in 1:p){
    for(j in 1:p){
      if ((i == 1) & (j == 1)){
        loss.record[1] <- sqrt(sum(Beta[1, i, j, ] ^ 2))
      }
    }
  }
  
  # Solve beta hat by iterative approach
  for (t in 1:(n.iteration - 1)){
    beta.temp <- Beta[t, , , ]
    
    for (i in 1:p){
      for (j in 1:p){
        
        ## compute dij(t+1), denoted by d.temp
        for (i1 in 1:p){
          for (j1 in 1:p){
            if ((i1 != i) | (j1 != j)){
              d.temp <- d.temp <- beta.temp[i1, j1, ] %*% t(Aarr[i1, j1, , ])
            }
          }
        }
        d.temp <- bvec - d.temp
        
        ## compute betaij(t+1)
        A.temp <- Aarr[i, j, , ]
        beta.renew <-tryCatch(ConstrainedL2(A.temp, d.temp, lambda), 
                              error=function(e){return (beta.temp[i, j, ])})
        beta.temp[i, j, ] <- beta.renew
      }
    }
    
    Beta[t + 1, , , ] <- beta.temp
    
    # renew beta record
    for(i in 1:p){
      for (j in 1:p){
        if ((i == 1) & (j == 1)){
          beta.record[t + 1, (((i - 1) * p + j - 1) * (M ^ 2) + 1):(((i - 1) * p + j) * (M ^ 2))] <- Beta[t + 1, i, j, ]
        }
      }
    }
    
    # renew loss record
    for(i in 1:p){
      for(j in 1:p){
        if ((i == 1) & (j == 1)){
          loss.record[t+1] <- sqrt(sum(Beta[t+1, i, j, ] ^ 2))
        }
      }
    }
    
    # check if convergence
    if (t > 3){
      if ((loss.record[t + 1] - loss.record[t - 2]) < lambda * stop.criterion){
        return(list(betahat=Beta[t + 1, , , ], betarecord=beta.record[1:(t + 1), ], lossrecord=loss.record[1:(t + 1)], numiter=(t + 1)))
      }
    }
    
  }
  return(list(betahat=Beta[t + 1, , , ], betarecord=beta.record, lossrecord=loss.record, numiter=(t + 1)))
}

# test
#p <- 2
#M <- 2
#Aarr <- array(0, c(p, p, ((p ^ 2) * (M ^ 2)), (M ^ 2)))
#Aarr[1, 1, , ] <- rbind(diag(1, 4), matrix(0, nrow=12, ncol=4))
#Aarr[1, 2, , ] <- rbind(matrix(0, nrow=4, ncol=4), diag(1, 4), matrix(0, nrow=8, ncol=4))
#Aarr[2, 1, , ] <- rbind(matrix(0, nrow=8, ncol=4), diag(1, 4), matrix(0, nrow=4, ncol=4))
#Aarr[2, 2, , ] <- rbind(matrix(0, nrow=12, ncol=4), diag(1, 4))
#bvec <- c(rep(1, 4), rep(0, 12))

#g <- IRGL2(Aarr, bvec, 0.5, p, M)
#g$numiter
#g$betarecord
#g$lossrecord
#for(i in 1:p){
#  for(j in 1:p){
#    print(g$betahat[i, j, ])
#  }
#}
