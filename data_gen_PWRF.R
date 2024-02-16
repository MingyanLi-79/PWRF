# *****************************
#     Date Modified: 12/16    *
# *****************************
# --------------- Generate X & Y, contaminate train Y --------------------
data_gen <- function(nobs, distX = "normal", p, m = NULL,
                     train_size, sd_gam = NULL, example = c(1, 2, 3, 4)) {
  # -- Dimensionality of covariates
  if (example == 1|example == 2) {dimX = 6} else {dimX = 10}
  # -- X
  if (distX == "normal") {
    if (example == 4) {
      X = MASS::mvrnorm(n = nobs, mu=c(rep(0,dimX)), Sigma=toeplitz(0.7^seq.int(0, dimX-1)))
      } else {
        X = matrix(rnorm(nobs*dimX), ncol = dimX)
      }
    TRAINX = X[1:train_size,]
    TESTX = X[(train_size+1):nobs,]
  }
  # -- epsilon -> uncontaminated noise
  eps = rnorm(nobs)
  # ---- Example 1 Y ----
  if (example == 1){
    # -- Y in different signal-to-noise ratio
    I = matrix(rep(0, nobs*7), nrow = nobs)
    I[(X[,1] <= 0 & X[,2] <= 0), 1] <- 1
    I[(X[,1] <= 0 & X[,2] > 0 & X[,4] <= 0), 2] <- 1
    I[(X[,1] <= 0 & X[,2] > 0 & X[,4] > 0 & X[,6] <= 0), 3] <- 1
    I[(X[,1] <= 0 & X[,2] > 0 & X[,4] > 0 & X[,6] > 0), 4] <- 1
    I[(X[,1] <= 0 & X[,3] <= 0), 5] <- 1
    I[(X[,1] <= 0 & X[,3] <= 0 & X[,5] <= 0), 6] <- 1
    I[(X[,1] <= 0 & X[,3] <= 0 & X[,5] > 0), 7] <- 1 # ;I
    coeff = seq(1:7)
    # -- Y
    Y = m*(I%*%coeff) + eps
    # end of example 1
  # ---- Example 2 Y ----
  } else if (example == 2){
    # -- Y in different signal-to-noise ratio
    I = matrix(rep(0, nobs*7), nrow = nobs)
    I[, 1] <- X[,1]
    I[, 2] <- (X[,2])^2
    I[(X[,3] > 0), 3] <- 1
    I[, 4] <- log(abs(X[,1]))*X[,3]
    I[, 5] <- X[,2]*X[,4]
    I[(X[,5] > 0), 6] <- 1
    I[, 7] <- exp(X[,6]) # ;I
    coeff = c(1, 0.707, 1, 0.873, 0.894, 2, 0.464)
    # -- Y
    Y = m*(I%*%coeff) + eps
    # end of example 2
  # ---- Example 3&4 Y ----
  } else {
    # -- Y
    Y = apply(X^2, 1, sum) + eps
    # end of example 3&4
  # ---- Example 4 Y ----
  } 
  # split train and test
  TRAINY = Y[1:train_size]
  TESTY = Y[(train_size+1):nobs]
  # Contaminate Train only with mean shift
  con_id = sort(sample(1:train_size, ceiling(p*train_size)))
  TRAINY[con_id] = TRAINY[con_id] + 3*max(TRAINY)
  return(list(TRAINX = TRAINX, TRAINY = TRAINY, TESTX = TESTX, TESTY = TESTY))
}


# -----------------------------------------------------------------------------------------------
dt <- data_gen(nobs = 15, p = 0.2, train_size = 10, example = 4)
  
  # , sd_gam = 5
  
  