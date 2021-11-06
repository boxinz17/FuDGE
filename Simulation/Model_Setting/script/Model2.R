p.v <- c(30, 60, 90, 120)

save.path <- ""
func.path <- ""

source(paste(func.path,"Graph_Genertation_Model2_zerodiag.R", sep="/"))

for (p in p.v){
  g <- Graph_Generation_Model2(p, m=5, s=4, seed=233)
  save(g, file=paste(save.path, "/model2_p", p, ".RData", sep=""))
}

for (p in p.v){
  load(paste(save.path, "/model2_p", p, ".RData", sep=""))
  print(is.positive.definite(g$X))
  print((sum(g$SupportX)-p)/2)
  print(is.positive.definite(g$Y))
  print((sum(g$SupportY)-p)/2)
  print(sum(g$SupportDelta)/2)
}
