#' @description Wrapper function for the Kalman filter
#' @param spec character indicating the model specification
#' @param param vector of structural parameters
#' @param data matrix with the data. One column per observed variable
#' @param outLogLik boolean. If TRUE, the function returns only the log likelihood and no filtered
#' values. If FALSE, the function returns the filtered state vector
#' @param constrainParam booelan. If TRUE, a parameter constraining function is applied
#' @return either the filtered output or the log likelihood

KalmanFilter <- function(spec, param, data, outLogLik, constrainParam) {
  # Apply parameter constraints if necessary
  if (constrainParam == TRUE) param <- ParConstrain(spec = spec, paramVec = param)
  # Set up the filter
  systemList <- SystemMat_fctn(spec = spec, paramVec = param)
  # Run the recursions
  Output <- KalmanRecursions(
    data, systemList, paramVec, outLogLik
  )
  return(Output)
}


#' @description Function that applies parameter constraints
#' @param spec character indicating the model specification
#' @param par vector structural parameters
#' @return constrained parameter vector

ParConstrain <- function(spec, paramVec) {
  constrPar <- paramVec

  # tbd

  return(constrPar)
}



#' @description Function constructs the MODEL SPECIFIC system matrices for the Kalman filter
#' @param spec character indicating the model specification
#' @param paramVec vector structural parameters
#' @return list with the matrices

SystemMat_fctn <- function(spec, paramVec) {
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

  # Dimensions
  transDim <- NCOL(A)
  obsDim <- NCOL(data)
  # Number of time periods
  nPeriods <- length(data)

  # Initialize the output
  logLik_mat <- matrix(0, nc = obsDim, nr = nPeriods)
  Z_tt_mat <- matrix(NA, nc = transDim, nr = nPeriods)
  P_tt_array <- array(NA, dim = c(transDim, transDim, nPeriods))
  K_t_array <- array(NA, dim = c(Dimans, transDim, nPeriods))
  epsilon_mat <- rep(NA, nc = obsDim, nr = nPeriods)
  Sigma_array <- rep(NA, dim = c(obsDim, obsDim, nPeriods))

  # Initialize the filter routine (diffusely)
  # State vector with zeros
  Z_t1 <- matrix(0, nrow = transDim, nc = 1)
  # State vector var-cov matrix with very high variance
  P_t1 <- diag(transDim) * 1000

  # Run the recursions (until nPeriods-1 bc state vector enters measurement eq with lag (eq. 2))
  for (i in 1:(nPeriods - 1)) {
    #-------------------#
    # Get innovations
    #-------------------#

    # One step ahead prediction error
    epsilon_t <- data[i + 1, ] - C %*% Z_tt
    # One step ahead prediction error
    Sigma_t <- C %*% P_tt %*% transC + D
    Sigma_inv_t <- Inverse(Sigma_t)
    
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
      logLik_mat[i, ] <- .5 * log(2 * pi) - .5 * log(abs(det(Sigma_t))) + (-.5 * epsilon %*% Sigma_inv_t %*% epsilon)
    } else {
      K_t_array[, , i] <- K_t
      Z_tt_mat[i,] <- Z_tt
      P_tt_array[, , i] <- P_tt
      epsilon_mat[i,] <- epsilon_t
      Sigma_array[,,i] <- Sigma_t
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
    return(-sum(logLik_mat))
  } else {
    # standardized prediction errors eq. (13) 
    # Bootstrap algorithm step 1 (note that first innovation is dropped)
    standFactor <- 1 / (nPeriods - 1) * apply(epsilon_mat[,-1], 2, sum)
    epsilonCenter_t <- epsilon_mat - standFactor
    epsilonCenter_list <- apply(epsilonCenter_t, 1, as.matrix)
    
    # Check if this really works
    
    e_hat_mat <- apply(Sigma_array^(-.5), c(1, 2), function(x){
      x %*% epsilonCenter_list[[as.integer(substring(deparse(substitute(x)), 2))]]
      })
    e_hat_mat[,1] <- NA
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

