p.v <- c(30, 60, 90, 120)

save.path <- ""  # path to save the model

func.path <- ""  # path saving function definitions

source(paste(func.path,"Graph_Genertation_PowLaw_zerodiag.R", sep="/"))

for (p in p.v){
  g <- Graph_Generation_PowLaw(p, m=5, alpha=2, seed=666)
  save(g, file=paste(save.path, "/model1_p", p, ".RData", sep=""))
}
# to load data, run for example:
# load(paste(save.path, "/model_powlaw_p", p, ".RData", sep=""))

for (p in p.v){
  load(paste(save.path, "/model1_p", p, ".RData", sep=""))
  print(is.positive.definite(g$X))
  print((sum(g$SupportX)-p)/2)
  print(is.positive.definite(g$Y))
  print((sum(g$SupportY)-p)/2)
  print(sum(g$SupportDelta)/2)
}
