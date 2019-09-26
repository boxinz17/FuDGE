# data preparation

timestart <- proc.time()

# load data and organize it for latter use
# path to store files
file.path <- ""  # file path
setwd(paste(file.path, "./eeg_full_unzip", sep=""))

load(paste(file.path, "/Position_list.Rdata", sep=""))  # load postion lists
load(paste(file.path, "/Inv_Position_list.Rdata", sep="")) # load inverse postion lists

# parameters setting
N.alco <- 77
N.contrl <- 45
p <- 64
N.sample <- 256

alco.array <- array(0, c(N.alco, p, N.sample))  # array to store alcoholic data
contrl.array <- array(0, c(N.contrl, p, N.sample))  # # array to store control data

count.alco <- 0
count.contrl <- 0

subjects <- dir()
for (sub in subjects){
  num.single.tr <- 0  # number of single obj treatments
  indvd.record <- list()  # temp list to store individual data
  setwd(paste(file.path, "./eeg_full_unzip/", sub, sep=""))
  trials <- dir()
  for (tri in trials){
    f <- file(paste(file.path, "./eeg_full_unzip/", sub, "/", tri, sep=""), open="r")
    lines <- readLines(f)
    len.lines <- length(lines)
    if (len.lines < 4){
      close(f)
      next
    }
    trt.type <- strsplit(lines[4], split=" ")[[1]][2]
    if(trt.type == "S1"){
      temp.matrix <- matrix(0, nrow=p, ncol=N.sample)  # temp matrix to store trial data
      num.single.tr <- num.single.tr + 1
      N.lines <- length(lines)
      num.vertex <- 0
      for (num.line in 1:N.lines){
        if (num.line > 4){
          if (strsplit(lines[num.line], split=" ")[[1]][1] == "#"){
            num.vertex <- num.vertex + 1
            next
          }
          num.col <- as.numeric(strsplit(lines[num.line], split=" ")[[1]][3]) + 1
          temp.matrix[num.vertex, num.col] <- as.numeric(strsplit(lines[num.line], 
                                                                  split=" ")[[1]][4])
        }
      }
      indvd.record[[num.single.tr]] <- temp.matrix
      close(f)
    } else {
      close(f)
      next
    }
  }
  
  indvd.array.temp <- array(0, c(num.single.tr, p, N.sample))
  for (i in 1:num.single.tr){
    indvd.array.temp[i, , ] <- indvd.record[[i]]
  }
  indvd.matrix <- apply(indvd.array.temp, c(2,3), mean)
  
  if (strsplit(sub, split="")[[1]][4] == 'a'){
    # data is from alcoholic group
    count.alco <- count.alco + 1
    alco.array[count.alco, , ] <- indvd.matrix
  } else if (strsplit(sub, split="")[[1]][4] == 'c'){
    count.contrl <- count.contrl + 1
    contrl.array [count.contrl, , ] <- indvd.matrix
  } else {
    stop("Not belong to a or c group")
  }
}

save(alco.array, file=paste(file.path, "/alco_array.Rdata", sep=""))
save(contrl.array, file=paste(file.path, "/contrl_array.Rdata", sep=""))

timeend <- proc.time()
runningtime <- timeend-timestart
print(runningtime[3])