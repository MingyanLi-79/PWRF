cvfun_ts <- function(train, validate, lams_PWRF, alphas_Huber, alphas_Tukey,
                     reps = 4, trim=0.2, plot = TRUE){ # i = 1; j = 1
  MSRE_PW=MSRE_Tukey=MSRE_Huber=array(0,dim=c(length(lams_PWRF),reps))
  for (i in 1:length(lams_PWRF)){
    for (j in 1:reps){
      rrf <- regression_forest(X=train[,-dim(train)[2]], Y=train[,dim(train)[2]],num.trees = 500)
      OOBWeights<- as.matrix(get_forest_weights(rrf))
      PredWeights <- as.matrix(get_forest_weights(rrf, validate[,-dim(validate)[2]]))
        
      LPred_i <- RobustRF(OOBWeights, PredWeights, TRAINY=train[,dim(train)[2]], lambda=lams_PWRF[i], method="PWLS", tol=10^-6)
      MSRE_PW[i,j]= mean( (validate[,dim(validate)[2]]-LPred_i[[1]])^2,trim=trim)
      LPred_i <- RobustRF(OOBWeights, PredWeights, TRAINY=train[,dim(train)[2]], alpha=alphas_Huber[i],  method="Huber", tol=10^-6)
      MSRE_Huber[i,j]= mean( (validate[,dim(validate)[2]]-LPred_i[[1]])^2,trim=trim)
      LPred_i <- RobustRF(OOBWeights, PredWeights, TRAINY=train[,dim(train)[2]], alpha=alphas_Tukey[i], method="Tukey", tol=10^-6)
      MSRE_Tukey[i,j]= mean( (validate[,dim(validate)[2]]-LPred_i[[1]])^2,trim=trim)
    }# j = 2
  }
  
  cvr_PW=apply(MSRE_PW, 1, mean)
  lamopt_PW=lams_PWRF[which.min(cvr_PW)]
  cvr_Huber=apply(MSRE_Huber, 1, mean)
  alpha_opt_Huber=alphas_Huber[which.min(cvr_Huber)]
  cvr_Tukey=apply(MSRE_Tukey, 1, mean)
  alpha_opt_Tukey=alphas_Tukey[which.min(cvr_Tukey)]
  
  if (plot){
    par(mfrow=c(1,3))
    plot(lams_PWRF,cvr_PW,xlab='lambda',ylab='cve', main="PW-RF")
    abline(v=lams_PWRF[which.min(cvr_PW)],lty=3,col='gray')
    plot(alphas_Huber,cvr_Huber,xlab='alpha',ylab='cve', main="Huber-RF")
    abline(v=alphas_Huber[which.min(cvr_Huber)],lty=3,col='gray')
    plot(alphas_Tukey,cvr_Tukey,xlab='alpha',ylab='cve', main="Tukey-RF")
    abline(v=alphas_Tukey[which.min(cvr_Tukey)],lty=3,col='gray')
  }
  return(list(lams_PWRF=lams_PWRF,alphas_Huber=alphas_Huber, alphas_Tukey = alphas_Tukey,
              cvr_PW=cvr_PW,cvr_Huber=cvr_Huber, cvr_Tukey=cvr_Tukey, 
              lamopt=lamopt_PW, alpha_opt_Huber=alpha_opt_Huber,alpha_opt_Tukey=alpha_opt_Tukey))
}