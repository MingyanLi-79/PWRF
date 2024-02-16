# *****************************
#     Date Modified: 12/15    *
# *****************************
cvfun=function(dt, lam_list, alpha_list, alpha_list_Tukey, 
               cvreps=4, cvfolds=5, trim=0.2, plot=FALSE){
  source("Partitiondata.R")
  source("compute_deltas_pwls.R")
  source("RobustRF.R")
  # set j =1, k = 1, dt = train
  n=dim(dt)[1]
  TrainInd=Partitiondata(ntrain=n, cvreps=cvreps, cvfolds=cvfolds)
  MSRE_PW=MSRE_Tukey=MSRE_Huber=array(0,dim=c(length(lam_list),cvreps,cvfolds))
  #three dim
  # original training matrix data with y and x as columns
  for (i in 1:length(lam_list)){
    for (j in 1:cvreps){ # added "in"
      Trainrep_j=TrainInd[j,,] # get a cvfold id 
      for (k in 1:cvfolds){ # Changed from "k=1: cvfolds" to "k in 1:cvfolds"
        valid_jk_id=as.vector(Trainrep_j[k,])
        train_jk_id=as.vector(Trainrep_j[-k,])
        train=data.frame(dt[train_jk_id,]) 
        test=data.frame(dt[valid_jk_id,]) 
        
        rrf <- regression_forest(X=train[,-dim(train)[2]], Y=train[,dim(train)[2]],num.trees = 500)
        OOBWeights<- as.matrix(get_forest_weights(rrf))
        PredWeights <- as.matrix(get_forest_weights(rrf, test[,-dim(test)[2]]))
        
        LPred_i <- RobustRF(OOBWeights, PredWeights, TRAINY=train[,dim(train)[2]], lambda=lam_list[i], method="PWLS", tol=10^-6)  #RFLOWESS calculations, capitalized RobustRF, capitalized Y, changed "PWRF" to "PWLS"
        MSRE_PW[i,j,k]= mean( (test[,dim(test)[2]]-LPred_i[[1]])^2,trim=trim)
        LPred_i <- RobustRF(OOBWeights, PredWeights, TRAINY=train[,dim(train)[2]], alpha=alpha_list[i],  method="Huber", tol=10^-6)  #RFLOWESS calculations, capitalized RobustRF
        #LPred_i <- LOWESSPred(OOBWeights, PredWeights, TRAINY=train$Y, alpha=alpha_list[i], method="Huber", tol=10^-6)  #RFLOWESS calculations, changed alpha_lis to alpha_list
        MSRE_Huber[i,j,k]= mean( (test[,dim(test)[2]]-LPred_i[[1]])^2,trim=trim)
        LPred_i <- RobustRF(OOBWeights, PredWeights, TRAINY=train[,dim(train)[2]], alpha=alpha_list_Tukey[i], method="Tukey", tol=10^-6)  #RFLOWESS calculations, capitalized RobustRF, changed alpha_lis to alpha_list
        MSRE_Tukey[i,j,k]= mean( (test[,dim(test)[2]]-LPred_i[[1]])^2,trim=trim) # Made all Y's capital, changed to Tukey
      }
    }
  }
  cvr_PW=apply(MSRE_PW, 1, mean)
  lamopt_PW=lam_list[which.min(cvr_PW)]
  cvr_Huber=apply(MSRE_Huber, 1, mean)
  alpha_opt_Huber=alpha_list[which.min(cvr_Huber)] # changed to alpha_list
  cvr_Tukey=apply(MSRE_Tukey, 1, mean)
  alpha_opt_Tukey=alpha_list_Tukey[which.min(cvr_Tukey)] # changed to alpha_list, might need to use a broader list for alpha_list and lamlist since it's optimized at the extreme
  if (plot){
    par(mfrow=c(1,3))
    plot(lam_list,cvr_PW,xlab='lambda',ylab='cve', main="PW-RF")
    abline(v=lam_list[which.min(cvr_PW)],lty=3,col='gray')
    plot(alpha_list,cvr_Huber,xlab='alpha',ylab='cve', main="Huber-RF")
    abline(v=alpha_list[which.min(cvr_Huber)],lty=3,col='gray')
    plot(alpha_list_Tukey,cvr_Tukey,xlab='alpha',ylab='cve', main="Tukey-RF")
    abline(v=alpha_list_Tukey[which.min(cvr_Tukey)],lty=3,col='gray')
  }
  return(list(lam_list=lam_list,alpha_list=alpha_list,cvr_PW=cvr_PW,cvr_Huber=cvr_Huber, 
              cvr_Tukey=cvr_Tukey, lamopt=lamopt_PW, alpha_opt_Huber=alpha_opt_Huber,
              alpha_opt_Tukey=alpha_opt_Tukey))
}