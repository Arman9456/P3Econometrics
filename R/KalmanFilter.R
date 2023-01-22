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
    data = data, paramVec = param, systemList = systemList, outLogLik = outLogLik
  )
  return(Output)
}

dataMat <- matrix(rnorm(10), nc = 2)
paramVec <- rep(0, 4)

#' @description Function that applies parameter constraints
#' @param spec character indicating the model specification
#' @param par vector structural parameters
#' @return constrained parameter vector

ParConstrain <- function(paramVec){
  constrPar <- paramVec

  # tbd

  return(constrPar)
}



#' @description Function constructs the MODEL SPECIFIC system matrices for the Kalman filter
#' @param spec character indicating the model specification
#' @param paramVec vector structural parameters
#' @return list with the matrices

SystemMat_fctn <- function(paramVec) {
  # tbd
set.seed(2)
  # Transition eq
  A <- matrix(rnorm(25), nr = 5)
  B <- diag(5)
  # Measurement eq
  C <- matrix(rnorm(10), nr = 2, nc = 5)
  D <- diag(2)

  outputList <- list("A" = A, "B" = B, "C" = C, "D" = D)
  return(outputList)
}


#' @description Function that executes the Kalman recursions. Notation as in Angelini et al. (2022).
#' @param data matrix with the data. One column per observed variable
#' @param paramVec vector of structural parameters
#' @param systemList List with system matrices
#' @return constrained parameter vector

RKalmanRecursions <- function(data, paramVec, systemList, outLogLik) {

  # Transition eq
  A <- systemList$A
  transA <- Transp(A)
  B <- systemList$B
  # Measurement eq
  C <- systemList$C
  transC <- Transp(C)
  D <- systemList$D

  # Dimensions
  transDim <- NCOL(A)
  obsDim <- NCOL(data)
  # Number of time periods
  nPeriods <- NROW(data)

  # Initialize the output
  logLik_vec <- rep(NA, nPeriods)
  Z_tt_mat <- matrix(NA, nc = transDim, nr = nPeriods)
  P_tt_array <- array(NA, dim = c(transDim, transDim, nPeriods))
  K_t_array <- array(NA, dim = c(transDim, obsDim, nPeriods))
  e_hat_mat <- matrix(NA, nc = obsDim, nr = nPeriods)

  # Initialize the filter routine (diffusely)
  # State vector with zeros
  Z_t1 <- matrix(0, nrow = transDim, nc = 1)
  # State vector var-cov matrix with very high variance
  P_t1 <- diag(transDim) * 1000

  # Run the recursions (until nPeriods-1 bc state vector enters measurement eq with lag (eq. 2))
  for (i in 1:(nPeriods)) {
    #-------------------#
    # Get innovations
    #-------------------#

    # One step ahead prediction error
    epsilon_t <- data[i, ] - C %*% Z_t1
    # One step ahead prediction error
    Sigma_t <- C %*% P_t1 %*% transC + D
    Sigma_inv_t <- Inverse(Sigma_t)
    # standardized prediction errors eq. (13)
    # Bootstrap algorithm step 1
    e_hat_t <- diag((diag(Sigma_t))^(-.5)) %*% epsilon_t

    #-------------------#
    # Updating step
    #-------------------#
    
    # Kalman gain
    K_t <- P_t1 %*% transC %*% Sigma_inv_t
    # State vector
    Z_tt <- Z_t1 + K_t %*% epsilon_t
    # State vector var-cov matrix
    P_tt <- P_t1 - K_t %*% C %*% P_t1

    # Either store the likelihood or the filter output
    if (outLogLik == TRUE) {
      # Calculation and storage of the log likelihood
      logLik_vec[i] <- .5 * log(2 * pi) - .5 * log(abs(det(Sigma_t))) + (-.5 * Transp(epsilon_t) %*% Sigma_inv_t %*% epsilon_t)
    } else {
      K_t_array[, , i] <- K_t
      Z_tt_mat[i, ] <- Z_tt
      P_tt_array[, , i] <- P_tt
      e_hat_mat[i, ] <- e_hat_t
    }

    #-------------------#
    # Prediction step
    #-------------------#

    # State vector
    Z_t1 <- A %*% Z_tt
    # State vector var-cov matrix
    P_t1 <- A %*% P_tt %*% transA + B
  }
  # Set the output
  if (outLogLik == TRUE) {
    return(-sum(logLik_vec))
  } else {
    # Note that first innovation is dropped (eq 13)
    e_hat_mat[1, ] <- NA
    # Construct the output list
    outputList <- list(
      "Z_tt" = Z_tt_mat,
      "P_tt" = P_tt_array,
      "K_t" = K_t_array,
      "e_hat" = e_hat_mat
    )
    return(outputList)
  }
}
