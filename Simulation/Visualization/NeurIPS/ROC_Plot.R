save.path <- ""  # ROC curve saving path

func.path <- ""  # function definitions path

##################################### Model 1 #################################################
DFGM.result.path <- ""  # FuDGE result for model1 path

Multiple.result.path <- ""  # Multipe result for model1 path

Single.result.path <- ""  # Single result for model1 path

model.path <- ""  # model 1 saving path

source(paste(func.path, "AUC_Func.R", sep="/"))
source(paste(func.path, "Avg_ROC.R", sep="/"))
source(paste(func.path,"blockwise_Frob.R", sep="/"))

n.sample <- 100

M <- 30  # number of simulations

obo <- 15

p.v <- c(30, 60, 90, 120)

v.FPR <- seq(0, 1, 0.01)

thresh <- 0.01

##################### First for DFGM mehod #################################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(DFGM.result.path, "/Para_result_p", p, "_runind", i, ".RData", sep=""))
    ROC_list[[i]] <- Paral.result
    AUC.v[i] <- AUC_Func(Paral.result[, 2], Paral.result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- median(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.my <- data.frame(p=p.v, Mean.my=AUC.mean.v, Standard_Error.my=AUC.sd.v)

L.v.TPR.my.mod1 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.my.mod1[[i]] <- v.TPR
}

################# Second for AIC selected Single Method ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(Single.result.path, "/ThetaEst_Model1_p", p, "_runind", i, ".RData", sep=""))
    load(paste(Single.result.path, "/Otherinfo_Model1_p", p, "_runind", i, ".RData", sep=""))
    
    SX <- Otherinfo$SX
    SY <- Otherinfo$SY
    Theta.L <- Otherinfo$Theta.L
    gam.v <- Otherinfo$gam.v
    M.selected <- Otherinfo$M.selected
    
    # Plot ROC by using AIC
    N <- length(Paral.result)
    v.AIC.X <- numeric(N)
    v.AIC.Y <- numeric(N)
    
    for (j in 1:N){
      result <- Paral.result[[j]]
      gX <- result$gX
      gY <- result$gY
      
      v.AIC.X[j] <- n.sample*(-log(det(gX$ThetaMathat)) + 
                                sum(diag(SX %*% gX$ThetaMathat))) + sum(gX$Support)*(M.selected^2)
      v.AIC.Y[j] <- n.sample*(-log(det(gY$ThetaMathat)) + 
                                sum(diag(SY %*% gY$ThetaMathat))) + sum(gY$Support)*(M.selected^2)
    }
    
    best.index.X <- which.min(v.AIC.X)
    best.index.Y <- which.min(v.AIC.Y)
    
    Theta.hat.final.X <- Paral.result[[best.index.X]]$gX$ThetaMathat
    Theta.hat.final.Y <- Paral.result[[best.index.Y]]$gY$ThetaMathat
    
    Supp.Mat <- Theta.L$SupportDelta
    
    diff.Mat <- abs(Theta.hat.final.X - Theta.hat.final.Y)
    FrobMat <- blockwise_Frob(diff.Mat, M.selected)
    
    lambda.v <- seq(2*max(diff.Mat), 0, length.out=200)
    ROC_result <- matrix(NA, nrow=length(lambda.v), ncol=2)
    colnames(ROC_result) <- c("TPR", "FPR")
    count2 <- 0
    for (lambda in lambda.v){
      count2 <- count2 + 1
      
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[which(FrobMat >= lambda)] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.Mat[k, (k + 1):p])
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
      
      ROC_result[count2, 1] <- TPR
      ROC_result[count2, 2] <- FPR
    }
    
    ROC_list[[i]] <- ROC_result
    AUC.v[i] <- AUC_Func(ROC_result[, 2], ROC_result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.AIC <- data.frame(Mean.AIC=AUC.mean.v, Standard_Error.AIC=AUC.sd.v)

L.v.TPR.AIC.mod1 <- list()
len.p.v <- length(p.v)
for (i in 1:len.p.v){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.AIC.mod1[[i]] <- v.TPR
}

################# Third for BIC selected Single Method ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(Single.result.path, "/ThetaEst_Model1_p", p, "_runind", i, ".RData", sep=""))
    load(paste(Single.result.path, "/Otherinfo_Model1_p", p, "_runind", i, ".RData", sep=""))
    
    SX <- Otherinfo$SX
    SY <- Otherinfo$SY
    Theta.L <- Otherinfo$Theta.L
    gam.v <- Otherinfo$gam.v
    M.selected <- Otherinfo$M.selected
    
    # Plot ROC by using BIC
    N <- length(Paral.result)
    v.BIC.X <- numeric(N)
    v.BIC.Y <- numeric(N)
    
    for (j in 1:N){
      result <- Paral.result[[j]]
      gX <- result$gX
      gY <- result$gY
      
      v.BIC.X[j] <- n.sample*(-log(det(gX$ThetaMathat)) + 
                                sum(diag(SX %*% gX$ThetaMathat))) + sum(gX$Support)*log(n.sample)*(M.selected^2)
      v.BIC.Y[j] <- n.sample*(-log(det(gY$ThetaMathat)) + 
                                sum(diag(SY %*% gY$ThetaMathat))) + sum(gY$Support)*log(n.sample)*(M.selected^2)
    }
    
    best.index.X <- which.min(v.BIC.X)
    best.index.Y <- which.min(v.BIC.Y)
    
    Theta.hat.final.X <- Paral.result[[best.index.X]]$gX$ThetaMathat
    Theta.hat.final.Y <- Paral.result[[best.index.Y]]$gY$ThetaMathat
    
    Supp.Mat <- Theta.L$SupportDelta
    
    diff.Mat <- abs(Theta.hat.final.X - Theta.hat.final.Y)
    FrobMat <- blockwise_Frob(diff.Mat, M.selected)
    
    lambda.v <- seq(2*max(diff.Mat), 0, length.out=200)
    ROC_result <- matrix(NA, nrow=length(lambda.v), ncol=2)
    colnames(ROC_result) <- c("TPR", "FPR")
    count2 <- 0
    for (lambda in lambda.v){
      count2 <- count2 + 1
      
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[which(FrobMat >= lambda)] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.Mat[k, (k + 1):p])
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
      
      ROC_result[count2, 1] <- TPR
      ROC_result[count2, 2] <- FPR
    }
    
    ROC_list[[i]] <- ROC_result
    AUC.v[i] <- AUC_Func(ROC_result[, 2], ROC_result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.BIC <- data.frame(Mean.BIC=AUC.mean.v, Standard_Error.BIC=AUC.sd.v)

L.v.TPR.BIC.mod1 <- list()
len.p.v <- length(p.v)
for (i in 1:len.p.v){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.BIC.mod1[[i]] <- v.TPR
}

################# Fourth for Multiple Networks ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  load(paste(model.path, "/model1_p", p, ".RData", sep=""))
  True.Model <- g
  
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  
  for (i in 1:M){
    Supp.True <- True.Model$SupportDelta
    load(paste(Multiple.result.path, "/Sep_Time_Est_p", p, "_runind", i, ".RData", sep=""))
    
    ROC_table <- matrix(NA, nrow=length(Paral.result), ncol=2)
    colnames(ROC_table) <- c("TPR", "FPR")
    for (n.gam in 1:length(Paral.result)){
      Supp.Mat.sum <- matrix(0, nrow=p, ncol=p)
      for (t in 1:obo){
        Delta.List <- Paral.result[[n.gam]]
        Delta.hat <- Delta.List[[t]]
        Supp.Mat.temp <- matrix(0, nrow=p, ncol=p)
        Supp.Mat.temp[Delta.hat$blockFrob > 0] <- 1
        Supp.Mat.sum <- Supp.Mat.sum + Supp.Mat.temp
      }
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[Supp.Mat.sum >= obo/2] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.True[k, (k + 1):p])
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
      
      ROC_table[n.gam, 1] <- TPR
      ROC_table[n.gam, 2] <- FPR
    }
    
    ROC_list[[i]] <- ROC_table
    AUC.v[i] <- AUC_Func(ROC_table[, 2], ROC_table[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.Mult <- data.frame(Mean.Mult=AUC.mean.v, Standard_Error.Mult=AUC.sd.v)

L.v.TPR.Mult.mod1 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.Mult.mod1[[i]] <- v.TPR
}

################# Finally, combine result ########################

##################################### Model 2 #################################################
DFGM.result.path <- ""  # FuDGE result for model2 path

Multiple.result.path <- ""  # Multipe result for model2 path

Single.result.path <- ""  # Single result for model2 path

model.path <- ""  # model 1 saving path

source(paste(func.path, "AUC_Func.R", sep="/"))
source(paste(func.path, "Avg_ROC.R", sep="/"))
source(paste(func.path,"blockwise_Frob.R", sep="/"))

n.sample <- 100

M <- 30  # number of simulations

obo <- 15

p.v <- c(30, 60, 90, 120)

v.FPR <- seq(0, 1, 0.01)

thresh <- 0.01

##################### First for DFGM mehod #################################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(DFGM.result.path, "/Para_result_p", p, "_runind", i, ".RData", sep=""))
    ROC_list[[i]] <- Paral.result
    AUC.v[i] <- AUC_Func(Paral.result[, 2], Paral.result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- median(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.my <- data.frame(p=p.v, Mean.my=AUC.mean.v, Standard_Error.my=AUC.sd.v)

L.v.TPR.my.mod2 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.my.mod2[[i]] <- v.TPR
}

################# Second for AIC selected Single Method ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(Single.result.path, "/ThetaEst_Model1_p", p, "_runind", i, ".RData", sep=""))
    load(paste(Single.result.path, "/Otherinfo_Model1_p", p, "_runind", i, ".RData", sep=""))
    
    SX <- Otherinfo$SX
    SY <- Otherinfo$SY
    Theta.L <- Otherinfo$Theta.L
    gam.v <- Otherinfo$gam.v
    M.selected <- Otherinfo$M.selected
    
    # Plot ROC by using AIC
    N <- length(Paral.result)
    v.AIC.X <- numeric(N)
    v.AIC.Y <- numeric(N)
    
    for (j in 1:N){
      result <- Paral.result[[j]]
      gX <- result$gX
      gY <- result$gY
      
      v.AIC.X[j] <- n.sample*(-log(det(gX$ThetaMathat)) +
                                sum(diag(SX %*% gX$ThetaMathat))) + sum(gX$Support)*(M.selected^2)
      v.AIC.Y[j] <- n.sample*(-log(det(gY$ThetaMathat)) +
                                sum(diag(SY %*% gY$ThetaMathat))) + sum(gY$Support)*(M.selected^2)
    }
    
    best.index.X <- which.min(v.AIC.X)
    best.index.Y <- which.min(v.AIC.Y)
    
    Theta.hat.final.X <- Paral.result[[best.index.X]]$gX$ThetaMathat
    Theta.hat.final.Y <- Paral.result[[best.index.Y]]$gY$ThetaMathat
    
    Supp.Mat <- Theta.L$SupportDelta
    
    diff.Mat <- abs(Theta.hat.final.X - Theta.hat.final.Y)
    FrobMat <- blockwise_Frob(diff.Mat, M.selected)
    
    lambda.v <- seq(2*max(diff.Mat), 0, length.out=200)
    ROC_result <- matrix(NA, nrow=length(lambda.v), ncol=2)
    colnames(ROC_result) <- c("TPR", "FPR")
    count2 <- 0
    for (lambda in lambda.v){
      count2 <- count2 + 1
      
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[which(FrobMat >= lambda)] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.Mat[k, (k + 1):p])
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
      
      ROC_result[count2, 1] <- TPR
      ROC_result[count2, 2] <- FPR
    }
    
    ROC_list[[i]] <- ROC_result
    AUC.v[i] <- AUC_Func(ROC_result[, 2], ROC_result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.AIC <- data.frame(Mean.AIC=AUC.mean.v, Standard_Error.AIC=AUC.sd.v)

L.v.TPR.AIC.mod2 <- list()
len.p.v <- length(p.v)
for (i in 1:len.p.v){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.AIC.mod2[[i]] <- v.TPR
}

################# Third for BIC selected Single Method ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(Single.result.path, "/ThetaEst_Model1_p", p, "_runind", i, ".RData", sep=""))
    load(paste(Single.result.path, "/Otherinfo_Model1_p", p, "_runind", i, ".RData", sep=""))
    
    SX <- Otherinfo$SX
    SY <- Otherinfo$SY
    Theta.L <- Otherinfo$Theta.L
    gam.v <- Otherinfo$gam.v
    M.selected <- Otherinfo$M.selected
    
    # Plot ROC by using BIC
    N <- length(Paral.result)
    v.BIC.X <- numeric(N)
    v.BIC.Y <- numeric(N)
    
    for (j in 1:N){
      result <- Paral.result[[j]]
      gX <- result$gX
      gY <- result$gY
      
      v.BIC.X[j] <- n.sample*(-log(det(gX$ThetaMathat)) + 
                                sum(diag(SX %*% gX$ThetaMathat))) + sum(gX$Support)*log(n.sample)*(M.selected^2)
      v.BIC.Y[j] <- n.sample*(-log(det(gY$ThetaMathat)) + 
                                sum(diag(SY %*% gY$ThetaMathat))) + sum(gY$Support)*log(n.sample)*(M.selected^2)
    }
    
    best.index.X <- which.min(v.BIC.X)
    best.index.Y <- which.min(v.BIC.Y)
    
    Theta.hat.final.X <- Paral.result[[best.index.X]]$gX$ThetaMathat
    Theta.hat.final.Y <- Paral.result[[best.index.Y]]$gY$ThetaMathat
    
    Supp.Mat <- Theta.L$SupportDelta
    
    diff.Mat <- abs(Theta.hat.final.X - Theta.hat.final.Y)
    FrobMat <- blockwise_Frob(diff.Mat, M.selected)
    
    lambda.v <- seq(2*max(diff.Mat), 0, length.out=200)
    ROC_result <- matrix(NA, nrow=length(lambda.v), ncol=2)
    colnames(ROC_result) <- c("TPR", "FPR")
    count2 <- 0
    for (lambda in lambda.v){
      count2 <- count2 + 1
      
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[which(FrobMat >= lambda)] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.Mat[k, (k + 1):p])
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
      
      ROC_result[count2, 1] <- TPR
      ROC_result[count2, 2] <- FPR
    }
    
    ROC_list[[i]] <- ROC_result
    AUC.v[i] <- AUC_Func(ROC_result[, 2], ROC_result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.BIC <- data.frame(Mean.BIC=AUC.mean.v, Standard_Error.BIC=AUC.sd.v)

L.v.TPR.BIC.mod2 <- list()
len.p.v <- length(p.v)
for (i in 1:len.p.v){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.BIC.mod2[[i]] <- v.TPR
}

################# Fourth for Multiple Networks ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  load(paste(model.path, "/model2_p", p, ".RData", sep=""))
  True.Model <- g
  
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  
  for (i in 1:M){
    Supp.True <- True.Model$SupportDelta
    load(paste(Multiple.result.path, "/Sep_Time_Est_p", p, "_runind", i, ".RData", sep=""))
    
    ROC_table <- matrix(NA, nrow=length(Paral.result)+1, ncol=2)
    colnames(ROC_table) <- c("TPR", "FPR")
    for (n.gam in 1:length(Paral.result)){
      Supp.Mat.sum <- matrix(0, nrow=p, ncol=p)
      for (t in 1:obo){
        Delta.List <- Paral.result[[n.gam]]
        Delta.hat <- Delta.List[[t]]
        Supp.Mat.temp <- matrix(0, nrow=p, ncol=p)
        Supp.Mat.temp[Delta.hat$blockFrob > 0] <- 1
        Supp.Mat.sum <- Supp.Mat.sum + Supp.Mat.temp
      }
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[Supp.Mat.sum >= obo/2] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.True[k, (k + 1):p])
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
      
      ROC_table[n.gam, 1] <- TPR
      ROC_table[n.gam, 2] <- FPR
    }
    ROC_table[n.gam+1, 1] <- 1
    ROC_table[n.gam+1, 2] <- 1
    
    ROC_list[[i]] <- ROC_table
    AUC.v[i] <- AUC_Func(ROC_table[, 2], ROC_table[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.Mult <- data.frame(Mean.Mult=AUC.mean.v, Standard_Error.Mult=AUC.sd.v)

L.v.TPR.Mult.mod2 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.Mult.mod2[[i]] <- v.TPR
}

##################################### Model 3 #################################################
DFGM.result.path <- ""  # FuDGE result for model3 path

Multiple.result.path <- ""  # Multipe result for model3 path

Single.result.path <- ""  # Single result for model3 path

model.path <- ""  # model 1 saving path

source(paste(func.path, "AUC_Func.R", sep="/"))
source(paste(func.path, "Avg_ROC.R", sep="/"))
source(paste(func.path,"blockwise_Frob.R", sep="/"))

n.sample <- 100

M <- 30  # number of simulations

obo <- 15

p.v <- c(30, 60, 90, 120)

v.FPR <- seq(0, 1, 0.01)

thresh <- 0.01

##################### First for DFGM mehod #################################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(DFGM.result.path, "/Para_result_p", p, "_runind", i, ".RData", sep=""))
    ROC_list[[i]] <- Paral.result
    AUC.v[i] <- AUC_Func(Paral.result[, 2], Paral.result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- median(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.my <- data.frame(p=p.v, Mean.my=AUC.mean.v, Standard_Error.my=AUC.sd.v)

L.v.TPR.my.mod3 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.my.mod3[[i]] <- v.TPR
}

################# Second for AIC selected Single Method ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(Single.result.path, "/ThetaEst_Model1_p", p, "_runind", i, ".RData", sep=""))
    load(paste(Single.result.path, "/Otherinfo_Model1_p", p, "_runind", i, ".RData", sep=""))
    
    SX <- Otherinfo$SX
    SY <- Otherinfo$SY
    Theta.L <- Otherinfo$Theta.L
    gam.v <- Otherinfo$gam.v
    M.selected <- Otherinfo$M.selected
    
    # Plot ROC by using AIC
    N <- length(Paral.result)
    v.AIC.X <- numeric(N)
    v.AIC.Y <- numeric(N)
    
    for (j in 1:N){
      result <- Paral.result[[j]]
      gX <- result$gX
      gY <- result$gY
      
      v.AIC.X[j] <- n.sample*(-log(det(gX$ThetaMathat)) +
                                sum(diag(SX %*% gX$ThetaMathat))) + sum(gX$Support)*(M.selected^2)
      v.AIC.Y[j] <- n.sample*(-log(det(gY$ThetaMathat)) +
                                sum(diag(SY %*% gY$ThetaMathat))) + sum(gY$Support)*(M.selected^2)
    }
    
    best.index.X <- which.min(v.AIC.X)
    best.index.Y <- which.min(v.AIC.Y)
    
    Theta.hat.final.X <- Paral.result[[best.index.X]]$gX$ThetaMathat
    Theta.hat.final.Y <- Paral.result[[best.index.Y]]$gY$ThetaMathat
    
    Supp.Mat <- Theta.L$SupportDelta
    
    diff.Mat <- abs(Theta.hat.final.X - Theta.hat.final.Y)
    FrobMat <- blockwise_Frob(diff.Mat, M.selected)
    
    lambda.v <- seq(2*max(diff.Mat), 0, length.out=1000)
    ROC_result <- matrix(NA, nrow=length(lambda.v), ncol=2)
    colnames(ROC_result) <- c("TPR", "FPR")
    count2 <- 0
    for (lambda in lambda.v){
      count2 <- count2 + 1
      
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[which(FrobMat >= lambda)] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.Mat[k, (k + 1):p])
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
      
      ROC_result[count2, 1] <- TPR
      ROC_result[count2, 2] <- FPR
    }
    
    ROC_list[[i]] <- ROC_result
    AUC.v[i] <- AUC_Func(ROC_result[, 2], ROC_result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.AIC <- data.frame(Mean.AIC=AUC.mean.v, Standard_Error.AIC=AUC.sd.v)

L.v.TPR.AIC.mod3 <- list()
len.p.v <- length(p.v)
for (i in 1:len.p.v){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.AIC.mod3[[i]] <- v.TPR
}

################# Third for BIC selected Single Method ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(Single.result.path, "/ThetaEst_Model1_p", p, "_runind", i, ".RData", sep=""))
    load(paste(Single.result.path, "/Otherinfo_Model1_p", p, "_runind", i, ".RData", sep=""))
    
    SX <- Otherinfo$SX
    SY <- Otherinfo$SY
    Theta.L <- Otherinfo$Theta.L
    gam.v <- Otherinfo$gam.v
    M.selected <- Otherinfo$M.selected
    
    # Plot ROC by using BIC
    N <- length(Paral.result)
    v.BIC.X <- numeric(N)
    v.BIC.Y <- numeric(N)
    
    for (j in 1:N){
      result <- Paral.result[[j]]
      gX <- result$gX
      gY <- result$gY
      
      v.BIC.X[j] <- n.sample*(-log(det(gX$ThetaMathat)) + 
                                sum(diag(SX %*% gX$ThetaMathat))) + sum(gX$Support)*log(n.sample)*(M.selected^2)
      v.BIC.Y[j] <- n.sample*(-log(det(gY$ThetaMathat)) + 
                                sum(diag(SY %*% gY$ThetaMathat))) + sum(gY$Support)*log(n.sample)*(M.selected^2)
    }
    
    best.index.X <- which.min(v.BIC.X)
    best.index.Y <- which.min(v.BIC.Y)
    
    Theta.hat.final.X <- Paral.result[[best.index.X]]$gX$ThetaMathat
    Theta.hat.final.Y <- Paral.result[[best.index.Y]]$gY$ThetaMathat
    
    Supp.Mat <- Theta.L$SupportDelta
    
    diff.Mat <- abs(Theta.hat.final.X - Theta.hat.final.Y)
    FrobMat <- blockwise_Frob(diff.Mat, M.selected)
    
    lambda.v <- seq(2*max(diff.Mat), 0, length.out=1000)
    ROC_result <- matrix(NA, nrow=length(lambda.v), ncol=2)
    colnames(ROC_result) <- c("TPR", "FPR")
    count2 <- 0
    for (lambda in lambda.v){
      count2 <- count2 + 1
      
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[which(FrobMat >= lambda)] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.Mat[k, (k + 1):p])
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
      
      ROC_result[count2, 1] <- TPR
      ROC_result[count2, 2] <- FPR
    }
    
    ROC_list[[i]] <- ROC_result
    AUC.v[i] <- AUC_Func(ROC_result[, 2], ROC_result[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.BIC <- data.frame(Mean.BIC=AUC.mean.v, Standard_Error.BIC=AUC.sd.v)

L.v.TPR.BIC.mod3 <- list()
len.p.v <- length(p.v)
for (i in 1:len.p.v){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.BIC.mod3[[i]] <- v.TPR
}

#################### Fourth for Multiple Networks ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  load(paste(model.path, "/model3_p", p, ".RData", sep=""))
  True.Model <- g
  
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  
  for (i in 1:M){
    Supp.True <- True.Model$SupportDelta
    load(paste(Multiple.result.path, "/Sep_Time_Est_p", p, "_runind", i, ".RData", sep=""))
    
    ROC_table <- matrix(NA, nrow=length(Paral.result)+1, ncol=2)
    colnames(ROC_table) <- c("TPR", "FPR")
    for (n.gam in 1:length(Paral.result)){
      Supp.Mat.sum <- matrix(0, nrow=p, ncol=p)
      for (t in 1:obo){
        Delta.List <- Paral.result[[n.gam]]
        Delta.hat <- Delta.List[[t]]
        Supp.Mat.temp <- matrix(0, nrow=p, ncol=p)
        Supp.Mat.temp[Delta.hat$blockFrob > 0] <- 1
        Supp.Mat.sum <- Supp.Mat.sum + Supp.Mat.temp
      }
      Edgeshat <- matrix(0, nrow=p, ncol=p)
      Edgeshat[Supp.Mat.sum >= obo/2] <- 1
      
      v.estimate <- c()
      v.true <- c()
      for (k in 1:(p - 1)){
        v.estimate <- c(v.estimate, Edgeshat[k, (k + 1):p])
        v.true <- c(v.true, Supp.True[k, (k + 1):p])
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
      
      ROC_table[n.gam, 1] <- TPR
      ROC_table[n.gam, 2] <- FPR
    }
    ROC_table[n.gam+1, 1] <- 1
    ROC_table[n.gam+1, 2] <- 1
    
    ROC_list[[i]] <- ROC_table
    AUC.v[i] <- AUC_Func(ROC_table[, 2], ROC_table[, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- mean(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.Mult <- data.frame(Mean.Mult=AUC.mean.v, Standard_Error.Mult=AUC.sd.v)

L.v.TPR.Mult.mod3 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.Mult.mod3[[i]] <- v.TPR
}

################################## Finally, combine result ##############################
p.v <- c(30, 60, 90, 120)
v.FPR <- seq(0, 1, 0.01)

png(file=paste(save.path, "/Average_ROC_Curve", ".png", sep=""), width=1280, height=1280)
layout(matrix(c(rep(14, 5), 1,2,3,4,13,5,6,7,8,13,9,10,11,12,13), ncol=4, byrow=F),
       heights=c(4, 4, 4, 4, 1), widths=c(1, 10, 10, 10))

for (i in 1:length(p.v)){
  if (i == 1){
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="", type="n",
         ylab="True Postive Rate", cex.main=4, cex.lab=1.8, cex.axis=2, main="Model 1")
    par(mar=c(5.1, 4.5, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod1[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod1[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod1[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod1[[i]], lty=4, col=4, lwd=3)
    mtext(text=paste("p=", p.v[i], sep=""), side=2.2, line=5, cex=2)
  } else if (i==4){
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="False Postive Rate", type="n",
         ylab="True Postive Rate", cex.main=4, cex.lab=2, cex.axis=2)
    par(mar=c(5.1, 4.5, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod1[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod1[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod1[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod1[[i]], lty=4, col=4, lwd=3)
    mtext(text=paste("p=", p.v[i], sep=""), side=2.2, line=5, cex=2)
  } else {
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="", type="n",
         ylab="True Postive Rate", cex.main=4, cex.lab=2, cex.axis=2)
    par(mar=c(5.1, 4.5, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod1[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod1[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod1[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod1[[i]], lty=4, col=4, lwd=3)
    mtext(text=paste("p=", p.v[i], sep=""), side=2.2, line=5, cex=2)
  }
}

for (i in 1:length(p.v)){
  if (i == 1){
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="", type="n",
         ylab="", cex.main=4, cex.lab=2, cex.axis=2, main="Model 2")
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod2[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod2[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod2[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod2[[i]], lty=4, col=4, lwd=3)
  } else if (i==4){
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="False Postive Rate", type="n",
         ylab="", cex.main=4, cex.lab=2, cex.axis=2)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod2[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod2[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod2[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod2[[i]], lty=4, col=4, lwd=3)
  } else {
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="", type="n",
         ylab="", cex.main=4, cex.lab=2, cex.axis=2)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod2[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod2[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod2[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod2[[i]], lty=4, col=4, lwd=3)
  }
}

for (i in 1:length(p.v)){
  if (i == 1){
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="", type="n",
         ylab="", cex.main=4, cex.lab=2, cex.axis=2, main="Model 3")
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod3[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod3[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod3[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod3[[i]], lty=4, col=4, lwd=3)
  } else if (i==4) {
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="False Postive Rate", type="n",
         ylab="", cex.main=2.5, cex.lab=2, cex.axis=2)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod3[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod3[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod3[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod3[[i]], lty=4, col=4, lwd=3)
  } else {
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="", type="n",
         ylab="", cex.main=2.5, cex.lab=2, cex.axis=2)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    lines(v.FPR, L.v.TPR.my.mod3[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR, L.v.TPR.AIC.mod3[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.BIC.mod3[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.Mult.mod3[[i]], lty=4, col=4, lwd=3)
  }
}

par(mai=c(0,0,0,0))
plot.new()
legend(x="center",ncol=4,legend=c("FuDGE", "AIC", "BIC", "Multiple"), 
       col=1:4, lty=1:4, lwd=4, cex=3)

par(mai=c(0,0,0,0))
plot.new()

dev.off()