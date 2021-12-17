func.path <- ""
save.path <- ""

source(paste(func.path,"Graph_Genertation_Example.R", sep="/"))
source(paste(func.path,"ProxAlg_DFGM.R", sep="/"))

m.v <- seq(10, 1010, 50)
signal_incomp.v <- numeric(length(m.v))  # signal strength of incomparable graphs
signal_board.v <- numeric(length(m.v))  # signal strength of boarderline case
signal_comp.v <- numeric(length(m.v))  # signal strength of comparable graphs

k_m <- 0
for (m in m.v){
  k_m <- k_m + 1
  g <- Graph_Generation_Model_Examp(m=m, alpha=1.5, beta=2.0, seed=996)
  blockwise.norm <- blockwise_Frob(g$X-g$Y, m)
  signal_incomp.v[k_m] <- blockwise.norm[1,2]
}

k_m <- 0
for (m in m.v){
  k_m <- k_m + 1
  g <- Graph_Generation_Model_Examp(m=m, alpha=1.5, beta=3.0, seed=996)
  blockwise.norm <- blockwise_Frob(g$X-g$Y, m)
  signal_board.v[k_m] <- blockwise.norm[1,2]
}

k_m <- 0
for (m in m.v){
  k_m <- k_m + 1
  g <- Graph_Generation_Model_Examp(m=m, alpha=1.5, beta=4.0, seed=996)
  blockwise.norm <- blockwise_Frob(g$X-g$Y, m)
  signal_comp.v[k_m] <- blockwise.norm[1,2]
}

save(signal_incomp.v, file=paste(save.path, "/signal_incomp.v.RData", sep=""))
save(signal_board.v, file=paste(save.path, "/signal_board.v.RData", sep=""))
save(signal_comp.v, file=paste(save.path, "/signal_comp.v.RData", sep=""))

load(file=paste(save.path, "/signal_incomp.v.RData", sep=""))
load(file=paste(save.path, "/signal_board.v.RData", sep=""))
load(file=paste(save.path, "/signal_comp.v.RData", sep=""))

setwd(save.path)
png("Signal_Strength_Example.png", width=1800, height=600)
par(mar = c(7.2, 8.2, 4.1, 2.1), mfrow=c(1,3), mgp=c(5.5, 1.7, 0))
plot(m.v, signal_incomp.v, xlab="M*", type="b", ylab="Signal Strength", 
     main="alpha=1.5, beta=2.0", cex.main=4.2, cex.lab=4.2, cex.axis=3.5, lwd=3, cex=3)
plot(m.v, signal_board.v, xlab="M*", type="b", ylab="", main="alpha=1.5, beta=3.0",
     cex.main=4.2, cex.lab=4.2, cex.axis=3.5, lwd=3, cex=3)
plot(m.v, signal_comp.v, xlab="M*", type="b", ylab="", main="alpha=1.5, beta=4.0",
     cex.main=4.2, cex.lab=4.2, cex.axis=3.5, lwd=3, cex=3)
dev.off()