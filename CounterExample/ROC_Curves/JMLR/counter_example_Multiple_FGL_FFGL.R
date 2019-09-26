DFGM.result.path <- ""  # Fudge result path

Multiple.result.path <- ""  # Multiple result path

FGL.result.path <- ""  # FGL result path

FFGL.result.path <- ""  # FFGL result path

FFGL2.result.path <- ""  # FFGL2 result path

save.path <- ""  # Save path

func.path <- ""  # Function definitions path

model.path <- ""  # Model saving path

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

AUC.table.my <- data.frame(p=p.v, Mean.my=AUC.mean.v, SE.my=AUC.sd.v)

L.v.TPR.my <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.my[[i]] <- v.TPR
}

################# Second for Multiple Networks ########################
AUC.mean.v <- numeric(length(p.v))
AUC.sd.v <- numeric(length(p.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(p.v))
colnames(AUC.record) <- p.v

ROC_list_list <- list()

count <- 0
for (p in p.v){
  load(paste(model.path, "/model_c_p", p, ".RData", sep=""))
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

AUC.table.Mult <- data.frame(Mean.Mult=AUC.mean.v, SE.Mult=AUC.sd.v)

L.v.TPR.Mult <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.Mult[[i]] <- v.TPR
}

################# Third for FGL Networks with lambda=0.01 ########################
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
    load(paste(FGL.result.path, "/Para_result_p", p, "_runind", i, "_lam1_0.01", ".RData", sep=""))
    ROC_list[[i]] <- Paral.result
    AUC.v[i] <- AUC_Func(Paral.result[, 1], Paral.result[, 2])
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

AUC.table.FGLLam001 <- data.frame(p=p.v, Mean.FGL=AUC.mean.v, SE.FGL=AUC.sd.v)

L.v.TPR.FGLLam001 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.FGLLam001[[i]] <- v.TPR
}

################# Fourth for FFGL Networks with lambda=0.1 ########################
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
    load(paste(FFGL.result.path, "/Para_result_p", p, "_runind", i, "_lam1_0.1", ".RData", sep=""))
    ROC_list[[i]] <- Paral.result
    AUC.v[i] <- AUC_Func(Paral.result[, 1], Paral.result[, 2])
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

AUC.table.FFGLLam01 <- data.frame(p=p.v, Mean.FFGL=AUC.mean.v, SE.FFGL=AUC.sd.v)

L.v.TPR.FFGLLam01 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.FFGLLam01[[i]] <- v.TPR
}

################# Fifth for FFGL2 Networks with lambda=0.1 ########################
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
    load(paste(FFGL2.result.path, "/Para_result_p", p, "_runind", i, "_lam1_0.1", ".RData", sep=""))
    ROC_list[[i]] <- Paral.result
    AUC.v[i] <- AUC_Func(Paral.result[, 1], Paral.result[, 2])
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

AUC.table.FFGL2Lam01 <- data.frame(p=p.v, Mean.FFGL2=AUC.mean.v, SE.FFGL2=AUC.sd.v)

L.v.TPR.FFGL2Lam01 <- list()
for (i in 1:length(p.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.FFGL2Lam01[[i]] <- v.TPR
}

################# Finally, combine result ########################
png(file=paste(save.path, "/ROC_Counter_Example.png", sep=""), width=1280, height=1080)
layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=T),
       heights=c(5,5,1), widths=c(1,1))
for (i in 1:length(p.v)){
  if (i==1){
    p <- p.v[i]
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), type="n", xlab="", 
         ylab="True Postive Rate", main=paste("p=", p, sep=""), cex.main=4, cex.lab=3, cex.axis=2)
    lines(v.FPR, L.v.TPR.my[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR,L.v.TPR.Mult[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.FGLLam001[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.FFGLLam01[[i]], lty=4, col=4, lwd=3)
    lines(v.FPR,L.v.TPR.FFGL2Lam01[[i]], lty=5, col=5, lwd=3)
  } else if (i==2) {
    p <- p.v[i]
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="", type="n", 
         ylab="", main=paste("p=", p, sep=""), cex.main=4, cex.lab=3, cex.axis=2)
    lines(v.FPR, L.v.TPR.my[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR,L.v.TPR.Mult[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.FGLLam001[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.FFGLLam01[[i]], lty=4, col=4, lwd=3)
    lines(v.FPR,L.v.TPR.FFGL2Lam01[[i]], lty=5, col=5, lwd=3)
  } else if (i==3) {
    p <- p.v[i]
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="False Postive Rate", type="n", 
         ylab="True Postive Rate", main=paste("p=", p, sep=""), cex.main=4, cex.lab=3, cex.axis=2)
    lines(v.FPR, L.v.TPR.my[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR,L.v.TPR.Mult[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.FGLLam001[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.FFGLLam01[[i]], lty=4, col=4, lwd=3)
    lines(v.FPR,L.v.TPR.FFGL2Lam01[[i]], lty=5, col=5, lwd=3)
  } else {
    p <- p.v[i]
    par(mar=c(5.1, 3.1, 4.1, 2.1))
    plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="False Postive Rate", type="n", 
         ylab="", main=paste("p=", p, sep=""), cex.main=4, cex.lab=3, cex.axis=2)
    lines(v.FPR, L.v.TPR.my[[i]], lty=1, col=1, lwd=3)
    lines(v.FPR,L.v.TPR.Mult[[i]], lty=2, col=2, lwd=3)
    lines(v.FPR,L.v.TPR.FGLLam001[[i]], lty=3, col=3, lwd=3)
    lines(v.FPR,L.v.TPR.FFGLLam01[[i]], lty=4, col=4, lwd=3)
    lines(v.FPR,L.v.TPR.FFGL2Lam01[[i]], lty=5, col=5, lwd=3)
  }
}

par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=4, 
       legend=c("FuDGE", "Multiple", expression(paste("FGL ",lambda, "1=0.01", sep="")), 
                expression(paste("FFGL ",lambda, "1=0.1", sep="")), 
                expression(paste("FFGL2 ",lambda, "1=0.1", sep=""))), 
       col=1:5, lty=1:5, lwd=4, cex=3)

dev.off()

AUC.table <- cbind(AUC.table.my, AUC.table.Mult, AUC.table.FGLLam001, 
                   AUC.table.FFGLLam01, AUC.table.FFGL2Lam01)
round(AUC.table, 2)

write.table(AUC.table, file=paste(save.path, "AUC_tabel_Counter_Example.txt", sep="/"))
write.csv(round(AUC.table, 2), file=paste(save.path, 
                                          "AUC_tabel_Counter_Example(readable).csv", sep="/"))

