# K-fold Cross Validation function
KFoldCVSplit <- function(m, u, k){
  # K-fold Cross Validation Split Function
  #
  # Args:
  #   m: the data matrix, we will split the data by row
  #   u: time points vector
  #   k: the number of folds
  #
  # Returns:
  #   A list with CV splited observation matrix and time points. 
  #   The first element is an array containing k matrix divied by m. 
  #   The second element is an array containg the splited time points.
  num.row <- dim(m)[1]
  if ((num.row %% k) != 0){
    stop("The number of rows cannot be divied by the number of folds, please check")
  } else {
    num.row.divided <- (num.row / k)
  }
  
  record.array <- array(0, dim=c(k, num.row.divided, dim(m)[2]))
  time.array <- array(0, dim=c(k, num.row.divided))
  index.remain <- c(1:num.row)
  
  for (i in c(1:k)){
    index.selected <- sample(1:((k - i + 1) * num.row.divided), num.row.divided, 
                             replace=FALSE)
    record.array[i, , ] <- array(m[index.selected, ], c(num.row.divided, dim(m)[2]))
    time.array[i, ] <- u[index.selected]
    m <- m[-index.selected, ]
  }
  
  return(list(fval=record.array, t=time.array))
}