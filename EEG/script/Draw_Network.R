# Draw the final picture

library("igraph")

file.path <- ""  # file path

load(paste(file.path, "/Position_list.Rdata", sep=""))  # load pos.list
load(paste(file.path, "/AdjMat.Rdata", sep=""))  # load AdjMat

p <- 64

diag(AdjMat) <- 0

node.names <- numeric(p)
for (i in 1:p){
  node.names[i] <- pos.list[[i]]
}

colnames(AdjMat) <- node.names
row.names(AdjMat) <- node.names

position.List <- list("FPZ"=c(0, 0.8), "AFZ"=c(0, 0.6), "FZ"=c(0, 0.4), FCZ=c(0, 0.2), 
                      "CZ"=c(0, 0), "CPZ"=c(0, -0.2), "PZ"=c(0, -0.4), "POZ"=c(0, -0.6), 
                      "OZ"=c(0, -0.8), "nd"=c(0,-1), "C2"=c(0.2, 0), "C4"=c(0.4, 0), 
                      "C6"=c(0.6, 0), "T8"=c(0.8, 0), "Y"=c(1, 0), "C1"=c(-0.2, 0), 
                      "C3"=c(-0.4, 0), "C5"=c(-0.6, 0), "T7"=c(-0.8, 0), "X"=c(-1, 0), 
                      "FP2"=c(0.2, 0.76), "FP1"=c(-0.2, 0.76), "AF2"=c(0.26, 0.62), 
                      "AF1"=c(-0.26, 0.62), "AF8"=c(0.45, 0.68), "AF7"=c(-0.45, 0.68), 
                      "F2"=c(0.21, 0.41), "F1"=c(-0.21, 0.41), "F4"=c(0.4, 0.45), 
                      "F3"=c(-0.4, 0.45), "F6"=c(0.55, 0.5), "F5"=c(-0.55, 0.5), 
                      "F8"=c(0.65, 0.55), "F7"=c(-0.65, 0.55), "FC2"=c(0.25, 0.21), 
                      "FC1"=c(-0.25, 0.21), "FC4"=c(0.5, 0.22), "FC3"=c(-0.5, 0.22), 
                      "FC6"=c(0.7, 0.26), "FC5"=c(-0.7, 0.26), "FT8"=c(0.9, 0.31), 
                      "FT7"=c(-0.9, 0.31), "CP2"=c(0.25, -0.21), "CP4"=c(0.5, -0.22), 
                      "CP6"=c(0.7, -0.26), "TP8"=c(0.9, -0.31), "CP1"=c(-0.25, -0.21), 
                      "CP3"=c(-0.5, -0.22), "CP5"=c(-0.7, -0.26), "TP7"=c(-0.9, -0.31), 
                      "P2"=c(0.21, -0.41), "P4"=c(0.4, -0.45), "P6"=c(0.55, -0.5), 
                      "P8"=c(0.65, -0.55), "P1"=c(-0.21, -0.41), "P3"=c(-0.4, -0.45), 
                      "P5"=c(-0.55, -0.5), "P7"=c(-0.65, -0.55), "PO2"=c(0.2, -0.62), 
                      "PO8"=c(0.45, -0.68), "PO1"=c(-0.2, -0.62), "PO7"=c(-0.45, -0.68), 
                      "O2"=c(0.2, -0.85), "O1"=c(-0.2, -0.85))

name.v <- c()

layMat <- matrix(NA, nrow=64, ncol=2)
for (i in 1:length(node.names)){
  x <- node.names[i]
  layMat[i, 1] <- unlist(position.List)[paste(x, 1, sep="")]
  layMat[i, 2] <- unlist(position.List)[paste(x, 2, sep="")]
}


jpeg(file=paste(file.path, "/Differential_Network.jpeg", sep=""), width=1280, height=1280)
net <- graph_from_adjacency_matrix(AdjMat, mode="undirected")
V(net)$label.cex <- 2.5
plot(net, edge.color="blue", vertex.size=2, edge.width=3, margin=0, layout=layMat, 
     edge.curved=1, vertex.label.dist=1, vertex.label.cex=3)
dev.off()
