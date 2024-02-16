Partitiondata=function(ntrain, cvreps, cvfolds){
  #set.seed(10)
  foldsize=ntrain/cvfolds
  TrainInd=array(dim=c(cvreps, cvfolds, foldsize))  
  for (rep in 1:cvreps){
    for (fold in 1:cvfolds){
      samp=sample(1:ntrain, ntrain, replace=FALSE)
      TrainInd[rep,fold,] = samp[(foldsize*(fold-1)+1):(fold*foldsize)]
    }
  }
  return(TrainInd)  
}