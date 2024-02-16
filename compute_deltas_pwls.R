Compute_deltas_PWLS=function(OOBWeights, TRAINY, alpha=6, lambda = 2, method="PWLS", 
                             tol=0.001, maxiter=100){
  #Before going into the loop, we set up the matrices and get the original OOB
  #prediction and residuals
  d <- rep(1, length(TRAINY))
  D0 <- matrix(rep(d, length(TRAINY)), nrow=length(TRAINY), ncol=length(TRAINY), byrow=T)
  AdjWts <- D0 * OOBWeights + 10^-9  # ensure weights not all 0 when alpha is very small
  AdjWts <- AdjWts / rowSums(AdjWts) # row average
  #OOBPred <- as.matrix(AdjWts) %*% TRAINY 
  OOBPred <- drop(as.matrix(AdjWts) %*% TRAINY)
  resid <- TRAINY- OOBPred
  resid[is.na(resid)] <- mean(abs(resid), na.rm=T) #if OOB is NA set residual to mean of residuals
  # initialize iteration
  niter <- 1
  Change <- 1
  while(Change > tol & niter < maxiter){  #main loop for update weights for the observations
    d0 <- d # this is the same logic from the function from RFLOWESS package, this is also
    #         why the function is initially returning d0 instead of d, again, the
    #         difference is not that big anyways.
    ## Method 1
    # d <- pmin(lambda/(10^-9+c((resid)^2)*OOBWeights) ,1) # will be a matrix
    # AdjWts <- d*OOBWeights + 10^-9 # add in very small buffer so that we don't divide by 0
    ## Method 2
    d<- pmin(lambda/(10^-9+(resid)^2) ,1) # PWLS loss
    D0 <- matrix(rep(d, length(TRAINY)), nrow=length(TRAINY), ncol=length(TRAINY), byrow=T)
    AdjWts <- D0*OOBWeights + 10^-9 
    AdjWts <- AdjWts / rowSums(AdjWts)
    OOBPred0 <- OOBPred
    OOBPred <- as.matrix(AdjWts) %*% TRAINY
    resid <- TRAINY - OOBPred
    resid[is.na(resid)] <- mean(abs(resid), na.rm=T)   # in case OOB NA
    Change <- max(abs(OOBPred-OOBPred0))
    niter <- niter + 1
  }
  return(list(d, OOBPred, niter))
}