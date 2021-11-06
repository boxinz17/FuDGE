result.path <- ""
func.path <- ""
save.path <- ""

source(paste(func.path, "AUC_Func.R", sep="/"))
source(paste(func.path, "Avg_ROC.R", sep="/"))
source(paste(func.path,"blockwise_Frob.R", sep="/"))
source(paste(func.path, "F1_score.R", sep="/"))

M <- 30  # number of simulations

p <- 30

m.v <- 5:9

basis.num.v <- matrix(0, nrow=M, ncol=length(m.v))  # matrix of chosen basis numbers
basis.num.v <- data.frame(basis.num.v)
colnames(basis.num.v) <- m.v

v.FPR <- seq(0, 1, 0.01)

#############################################################################

AUC.mean.v <- numeric(length(m.v))
AUC.sd.v <- numeric(length(m.v))

AUC.record <- matrix(NA, nrow=M, ncol=length(m.v))
colnames(AUC.record) <- m.v

thresh <- 0.1

ROC_list_list <- list()

count <- 0
for (m in m.v){
  count <- count + 1
  AUC.v <- numeric(M)
  ROC_list <- list()
  for (i in 1:M){
    load(paste(result.path, "/Para_result_m", m, "_runind", i, ".Rdata", sep=""))
    ROC_list[[i]] <- Paral.result[[1]]
    AUC.v[i] <- AUC_Func(Paral.result[[1]][, 2], Paral.result[[1]][, 1])
    if (AUC.v[i] < thresh){
      AUC.record[i, count] <- NA
    } else {
      AUC.record[i, count] <- AUC.v[i]
    }
    basis.num.v[i, count] <- Paral.result[[2]]
  }
  ROC_list_list[[count]] <- ROC_list
  AUC.mean.v[count] <- median(AUC.v[AUC.v>thresh])
  AUC.sd.v[count] <- sd(AUC.v[AUC.v>thresh])
}

AUC.table.my <- data.frame(m=m.v, Mean.my=AUC.mean.v, Standard_Error.my=AUC.sd.v)

L.v.TPR.my <- list()
for (i in 1:length(m.v)){
  ROC_list <- ROC_list_list[[i]]
  v.TPR <- Avg_ROC(ROC_list=ROC_list, v.FPR=v.FPR)
  L.v.TPR.my[[i]] <- v.TPR
}

#######################################################################
setwd(save.path)

png("ROC_basis_change.png", width = 960, height = 960)
par(mar = c(7, 7, 7, 3))
plot(c(0, 0), xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate", type="n", 
     ylab="True Positive Rate", 
     main="Average ROC curve for different basis numbers", 
     cex.main=3, cex.lab=3, cex.axis=3)
for (i in 1:length(m.v)){
  lines(v.FPR, L.v.TPR.my[[i]], col=i, lwd=4)
}
legend("bottomright", legend=c("m=5", "m=6", "m=7", "m=8", "m=9"), col=1:length(m.v), cex=3, lwd=4)
dev.off()

png("Boxplot.png")
boxplot(basis.num.v, xlab="Basis number to generate data", ylab="Basis number chosen by CV", 
        main="Boxplot for basis numbers chosen by CV")
dev.off()

write.table(AUC.table.my, "AUC-table.txt", sep="\t", row.names=FALSE)
