## Compute the maximum absolute difference of the top probes of two vectors

maxAbsDiffTopFeatures<-function(vectorX, vectorY, NumTopFeature=20){
  
  vectorX = sort(vectorX, decreasing = TRUE)
  vectorY = sort(vectorY, decreasing = TRUE)
  
  return( max(abs(vectorX[1:NumTopFeature] - vectorY[1:NumTopFeature])) )
}
