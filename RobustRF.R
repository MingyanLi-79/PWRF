RobustRF <- function(OOBWeights, PredWeights, TRAINY, alpha=6, lambda = 2, method="PWLS", 
                     tol=10^-6){
  if (method%in%c("Tukey", "Huber"))
    Res <- Compute_deltas(OOBWeights=OOBWeights, TRAINY=TRAINY, alpha=alpha, method=method, tol=tol)
  if (method=="PWLS")
    Res <- Compute_deltas_PWLS(OOBWeights=OOBWeights, TRAINY=TRAINY, lambda=lambda, method=method, tol=tol)
  d <- Res[[1]]
  niter <- Res[[3]]
  D <- matrix(rep(d, nrow(PredWeights)), nrow=nrow(PredWeights), ncol=length(TRAINY), byrow=T)
  AdjWts <- D*PredWeights
  AdjWts <- AdjWts/rowSums(AdjWts)
  Pred <- as.matrix(AdjWts)%*%TRAINY
  return(list(prediction=Pred, adjustedweight=AdjWts, iteration=niter,trainingweight=d))
}
