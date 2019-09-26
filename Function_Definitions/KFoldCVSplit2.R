# K-fold Cross Validation function
KFoldCVSplit2 <- function(m, u, k){
  # K-fold Cross Validation Split Function
  #
  # Args:
  #   m: the data matrix, we will split the data by row
  #   u: time points vector
  #   k: the number of folds
  #
  # Returns:
  #   A list with CV splited observation matrix and time points. 
  #   The first element is a list containing k matrix divied by m. 
  #   The second element is a list containg the splited time points.
  num.row <- dim(m)[1]
  
  num.row.divided <- floor(num.row/k)
  
  data.list <- list()
  time.list <- list()
  
  for (i in 1:k){
    if (i == k){
      data.list[[i]] <- m[((k-1)*num.row.divided+1):num.row, ]
      time.list[[i]] <- u[((k-1)*num.row.divided+1):num.row]
    } else {
      data.list[[i]] <- m[((k-1)*num.row.divided+1):(k*num.row.divided), ]
      time.list[[i]] <- u[((k-1)*num.row.divided+1):(k*num.row.divided)]
    }
  }
  
  return(list(fval=data.list, t=time.list))
}