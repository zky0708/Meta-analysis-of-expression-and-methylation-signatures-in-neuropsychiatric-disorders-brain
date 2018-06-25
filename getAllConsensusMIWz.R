## Compute the consensus MI of every row and the metagene across multiple data sets

library("limma")
library("cafr")

getAllConsensusMIWz <- function(dataList, metageneList, bin=6, so=3, negateMI = TRUE){

  n = length(dataList)
  if(length(metageneList) != n){
    stop("lengths of lists of dataset and metagene are different!")
  } 

  N = nrow(dataList[[1]])
  attractorMatrix <- matrix(0, n, N)
  for(i in 1:n){
    attractorMatrix[i,] <- getAllMIWz(dataList[[i]], metageneList[[i]], 
                                      bin=bin, so=so, sorting=FALSE, negateMI=negateMI)
    # Set negative MI elements to be zero
    attractorMatrix[i, attractorMatrix[i,]<0] <- 0
  }
      
  dataSetSize <- sapply(dataList,ncol)
  dataSetSize <- dataSetSize/sum(dataSetSize)

  w <- apply(attractorMatrix, 2, function(x) weighted.median(x, dataSetSize))
  
  return(w)
}
