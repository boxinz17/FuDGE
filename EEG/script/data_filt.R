# filt data with frequency filter

file.path <- ""  # file path
library(eegkit)
load(paste(file.path, "/alco_array.Rdata", sep=""))  # load alco.array
load(paste(file.path, "/contrl_array.Rdata", sep=""))  # load contrl.array

N.alco <- 77
N.contrl <- 45
p <- 64
N.sample <- 256

alco.filtered.array <- array(0, c(N.alco, p, N.sample))
contrl.filtered.array <- array(0, c(N.contrl, p, N.sample))

for (i in 1:N.alco){
  for (j in 1:p){
    alco.filtered.array[i, j, ] <- eegfilter(alco.array[i, j, ], 
                                           Fs=N.sample, lower=8, upper=12.5)
  }
}

for (i in 1:N.contrl){
  for (j in 1:p){
    contrl.filtered.array[i, j, ] <- eegfilter(contrl.array[i, j, ], 
                                           Fs=N.sample, lower=8, upper=12.5)
  }
}

save(alco.filtered.array, file=paste(file.path, "/alco_filtered_array.Rdata", sep=""))
save(contrl.filtered.array, file=paste(file.path, "/contrl_filtered_array.Rdata", sep=""))

