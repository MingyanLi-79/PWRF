Compute_deltas_PWLS=function(OOBWeights, TRAINY, alpha=6, lambda = 2, method="PWLS", 
                             tol=0.001, maxiter=100){
  d <- rep(1, length(TRAINY))
  D0 <- matrix(rep(d, length(TRAINY)), nrow=length(TRAINY), ncol=length(TRAINY), byrow=T)
  AdjWts <- D0 * OOBWeights + 10^-9  # ensure weights not all 0 when alpha is very small
  AdjWts <- AdjWts / rowSums(AdjWts) # row average
  #OOBPred <- as.matrix(AdjWts) %*% TRAINY 
  OOBPred <- drop(as.matrix(AdjWts) %*% TRAINY)
  resid <- TRAINY- OOBPred
  resid[is.na(resid)] <- mean(abs(resid), na.rm=T) #if OOB is NA set residual to mean of residuals
  niter <- 1
  Change <- 1
  while(Change > tol & niter < maxiter){  #main loop for update weights for the observations
    d0 <- d
    d<- pmin(lambda/(10^-9+(resid)^2) ,1)
    D0 <- matrix(rep(d, length(TRAINY)), nrow=length(TRAINY), ncol=length(TRAINY), byrow=T)
    AdjWts <- D0*OOBWeights + 10^-9 
    AdjWts <- AdjWts / rowSums(AdjWts)
    OOBPred0 <- OOBPred
    OOBPred <- as.matrix(AdjWts) %*% TRAINY
    resid <- TRAINY - OOBPred
    resid[is.na(resid)] <- mean(abs(resid), na.rm=T)
    Change <- max(abs(OOBPred-OOBPred0))
    niter <- niter + 1
  }
  return(list(d, OOBPred, niter))
}
