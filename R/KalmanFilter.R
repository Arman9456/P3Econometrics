#' @description Wrapper function for the Kalman filter
#' @param param vector of structural parameters
#' @param data matrix with the data. One column per observed variable
#' @param outLogLik boolean. If TRUE, the function returns only the log likelihood and no filtered
#' values. If FALSE, the function returns the filtered state vector
#' @param constrainParam booelan. If TRUE, a parameter constraining function is applied
#' @return either the filtered output or the log likelihood

KalmanFilter <- function(param, data, outLogLik, constrainParam) {
  # Apply parameter constraints if necessary
  if (constrainParam == TRUE) param <- ParConstrain(paramVec = param)
  # Set up the filter
  systemList <- SystemMat_fctn(paramVec = param)
  # Run the recursions
  Output <- KalmanRecursions(
    data, systemList, paramVec, outLogLik
  )
  return(Output)
}


#' @description Function that applies parameter constraints
#' @param par vector structural parameters
#' @return constrained parameter vector

ParConstrain <- function(paramVec) {
  constrPar <- paramVec

  # tbd

  return(constrPar)
}



#' @description Function constructs the MODEL SPECIFIC system matrices for the Kalman filter
#' @param paramVec vector structural parameters
#' @return list with the matrices

SystemMat_fctn <- function(paramVec) {
  # tbd

  # Transition eq
  A <- matrix(0, 2, 2)
  B <- matrix(0, 2, 2)
  # Measurement eq
  C <- matrix(0, 2, 2)
  D <- matrix(0, 2, 2)

  outputList <- list("A" = A, "B" = B, "C" = C, "D" = D)
  return(outputList)
}


#' @description Function that executes the Kalman recursions. Notation as in Angelini et al. (2022).
#' @param data matrix with the data. One column per observed variable
#' @param paramVec vector of structural parameters
#' @param systemList List with system matrices
#' @return constrained parameter vector

KalmanRecursions <- function(data, paramVec, systemList, outLogLik) {
  # Transition eq
  A <- systemList$A
  B <- systemList$B
  # Measurement eq
  C <- systemList$C
  transC <- trans(C)
  D <- systeMList$D

  # Dimensions of the state vector
  dimens <- NCOL(A)
  # Number of time periods
  nPeriods <- length(data)

  # Initialize the output
  logLik_T <- 0
  Z_tt_mat <- matrix(NA, nr = dimens, nc = nPeriods)
  P_tt_array <- array(NA, dim = c(dimens, dimens, nPeriods))
  K_t_array <- array(NA, dim = c(Dimans, Dimens, nPeriods))
  epsilonSt_vec <- rep(NA, nPeriods)

  # Initialize the filter routine (diffusely)
  # State vector with zeros
  Z_t1 <- matrix(0, nrow = Dimens, nc = 1)
  # State vector var-cov matrix with very high variance
  P_t1 <- diag(dimens) * 1000

  # Run the recursions (until nPeriods-1 bc state vector enters measurement eq with lag (eq. 2))
  for (i in 1:(nPeriods - 1)) {
    
    #-------------------#
    # Get innovations
    #-------------------#

    # One step ahead prediction error
    epsilon_t <- as.numeric(data[i + 1, ] - C %*% Z_tt)
    # One step ahead prediction error
    Sigma_t <- as.numeric(C %*% P_tt %*% transC + D)
    # standardized prediction errors
    epsilonSt_t <- Sigma_t^(-.5) * epsilon_t

    #-------------------#
    # Updating step
    #-------------------#

    # Kalman gain
    K_t <- P_t1 %*% transC %*% Inverse(Sigma_t)
    # State vector
    Z_tt <- Z_t1 + K_t %*% epsilon_t
    # State vector var-cov matrix
    P_tt <- P_t1 - K_t %*% C %*% P_t1

    # Either store the likelihood or the filter output
    if (outLogLik == TRUE) {
      # Calculation and storage of log-likelihood
      lik_t <- dnorm(epsilon_t, sd = sqrt(Sigma_t))
      logLik_T <- logLik_T - log(lik_t)
    } else {
      K_t_array[, , i] <- K_t
      Z_tt_mat[, i] <- Z_tt
      P_tt_array[, , i] <- P_tt
    }

    #-------------------#
    # Prediction step
    #-------------------#

    # State vector
    Z_t1 <- C %*% Z_tt
    # State vector var-cov matrix
    P_t1 <- C %*% P_tt %*% transC + B
    
  }
  # Set the output
  if (outLogLik == TRUE) {
    return(logLik_T)
  } else {
    outputList <- list(
      "Z_tt" = Z_tt_mat,
      "P_tt" = P_tt_array,
      "K_t" = K_t_array
    )
    return(outputList)
  }
}
