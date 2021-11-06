p.v <- c(30, 60, 90, 120)

save.path <- ""
func.path <- ""

source(paste(func.path,"Graph_Genertation_Model3_zerodiag.R", sep="/"))

for (p in p.v){
  if (p==30){
    s <- 3
  } else if (p==60) {
    s <- 4
  } else if (p==90) {
    s <- 5
  } else if (p==120) {
    s <- 6
  }
  
  g <- Graph_Generation_Model3(p, m=5, s=s, seed=996)
  save(g, file=paste(save.path, "/model3_p", p, ".RData", sep=""))
}

for (p in p.v){
  load(paste(save.path, "/model3_p", p, ".RData", sep=""))
  print(is.positive.definite(g$X))
  print((sum(g$SupportX)-p)/2)
  print(is.positive.definite(g$Y))
  print((sum(g$SupportY)-p)/2)
  print(sum(g$SupportDelta)/2)
}
