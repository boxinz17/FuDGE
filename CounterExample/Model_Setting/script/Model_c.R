p.v <- c(30, 60, 90, 120)

save.path <- ""  # model saving path

func.path <- ""  # function definitions path

source(paste(func.path,"Graph_Genertation_Model2.R", sep="/"))

for (p in p.v){
  g <- Graph_Generation_Model2(p, m=5, s=4)
  save(g, file=paste(save.path, "/model_c_p", p, ".RData", sep=""))
}
# to load data, run for example:
# load(paste(save.path, "/model_c_p", p, ".RData", sep=""))

for (p in p.v){
  load(paste(save.path, "/model_c_p", p, ".RData", sep=""))
  print(is.positive.definite(g$X))
  print((sum(g$SupportX)-p)/2)
  print(is.positive.definite(g$Y))
  print((sum(g$SupportY)-p)/2)
  print(sum(g$SupportDelta)/2)
}