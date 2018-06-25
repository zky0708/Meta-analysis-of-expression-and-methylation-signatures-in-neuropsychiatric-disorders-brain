### Consensus Probe-Selection Attractor Finding Algorithm for Multiple Seeds ###

## This algorithm will output the consensus attractors for assembled data sets based on the selected seeds
## The following code is based on the work of Tai-Hsien Ou Yang https://github.com/th86/gislkit


## There are two functions to be called in the main function.
# source("maxAbsDiffTopFeatures.R")
# source("getAllConsensusMIWz.R")

findMultipleConsPSAttractors <- function(dataMatrixList, map_short, seeds = NULL, 
                                          a = 5, maxIter = 500, epsilon=1E-7, 
                                          bin = 6, so = 3, negateMI = TRUE, 
                                          verbose = TRUE, stopWhenMaxIterReached = FALSE){
  
  ## dataMatrixList: list of multiple datasets with the same number of probes (rows)
  ## map_short: mapping each probe to a gene name. Note that multiple probes may be mapped to the same gene
  ## bin, so, negateMI: parameters used in function 'getAllConsensusMIWz()'
  ## seeds: probe names selected as seeds
  
  # Disregard probes not associated with gene symbols
  map_short <- map_short[which(is.na(map_short) == FALSE)]
  map_short <- map_short[which(names(map_short) %in% rownames(dataMatrixList[[1]]))]
  probeList <- names(map_short)
  names(probeList) <- map_short

  cat(length(probeList),"probes are included in the following scanning.\n")
  for(j in 1:length(dataMatrixList) ){
    dataMatrixList[[j]] <- dataMatrixList[[j]][probeList,]
  }
  
  # Check whether the provided seeds are valid 
  # (neither empty nor not applicable in the dataMatrixList)
  if(is.null(seeds) || any(! seeds %in% probeList))
    stop("Please provide a valid seed list.")
  else
    cat("There are", length(seeds), "seeds to search...\n")

 
  newAttractors <- NULL
  num <- 1
  
  while(length(seeds) > 0){
     
    # Initialize the 'metageneList' object with the vector of probes in 'seeds'
    metageneList <- list()																												  
    for(i in 1:length(dataMatrixList))
      metageneList[[i]]=dataMatrixList[[i]][seeds[1],]
     
    # Compute the consensus mutual information (MI) vector
    mi <- getAllConsensusMIWz(dataMatrixList, metageneList, 
                              bin = bin, so = so, negateMI = negateMI)
    names(mi) <- probeList
    
    
    ## Probe Selection
    # We select one and only one probe for each gene that has the highest consensus MI score
    gene_list <- unique(names(probeList))
    probe.filtered <- rep(0, length(gene_list))
    names(probe.filtered) <- gene_list
    for(i in 1:length(gene_list))
      probe.filtered[i] <- names(which.max(mi[which(names(probeList) == gene_list[i])]))
    
    # Update MI vector with selected probes and perform normalization
    mi <- mi[probe.filtered]
    w <- mi^a / sum( mi^a )                                                         				
    w <- w[probe.filtered]
    names(mi) <- gene_list

    # Update the metagenes using the weights vector
    for(i in 1:length(dataMatrixList))
      metageneList[[i]]= w %*% dataMatrixList[[i]][probe.filtered,]
    

    iter <- 1
    if(verbose==TRUE){
      cat("Iteration", iter, "\n")
      print(mi[order(mi, decreasing=TRUE)[1:10]])
    }   
    
    premi <- mi 

    # start of the inner while
    while(iter < maxIter){
      
      ## Iterate using the updated metagenes
      mi <- getAllConsensusMIWz(dataMatrixList, metageneList, 
                                bin = bin, so = so, negateMI = negateMI) 			
      names(mi) <- probeList
      
      probe.filtered=rep(0, length(gene_list))
      names(probe.filtered)=gene_list
      for(i in 1:length(gene_list))
        probe.filtered[i] <- names(which.max(mi[which(names(probeList) == gene_list[i])]))

      mi <- mi[probe.filtered]
      w <- mi^a / sum( mi^a )                                                       				
      w <- w[probe.filtered]
      names(mi) <- gene_list

      for( i in 1:length(dataMatrixList))
        metageneList[[i]]= w %*% dataMatrixList[[i]][probe.filtered, ]                              
       
      # Check whether reach the convergence threshold, 
      # and if so, then we can claim the attractor converges
      delta <- maxAbsDiffTopFeatures(mi, premi)
      if(verbose==TRUE){
        cat("Iteration ", (iter+1), "\tDelta = ", delta, "\n", sep="")
        print(mi[order(mi, decreasing=TRUE)[1:10]])
      }
      
      if(is.na(delta) == TRUE) break
      
      if(delta < epsilon){        
        names(mi) <- probe.filtered
        mi = sort(mi, decreasing = TRUE)
       
        break
        
      } else{

        premi <- mi
        iter <- iter + 1
        
      }
      
      # Stop iterating if reaches the maximum of iterations.
      if(iter >= maxIter & stopWhenMaxIterReached == TRUE)	
        stop(paste(seeds[1], "reached the maximum of iterations."))
      
    } #End of inner while
    
    tmp <- cbind(names(mi), map_short[names(mi)], round(mi, digits = 4))
    colnames(tmp) <- paste("Attractor", num, c("probes", "genes", "scores"), sep = "_")
    newAttractors <- cbind(newAttractors, tmp)
    
    seeds <- seeds[-1]
    num <- num + 1
    cat(paste(length(seeds), "seed(s) left \n", sep = " "))

  } #End of outer while
    
  # all seeds have been searched


  rownames(newAttractors) <- NULL
  return(newAttractors)
  
}
