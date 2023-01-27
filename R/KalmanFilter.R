#' @description Wrapper function for the Kalman filter. Executes the Kalman recursions
#' @param param vector of structural parameters
#' @param data matrix with the data. One column per observed variable
#' @param outLogLik boolean. If TRUE, the function returns only the log likelihood and no filtered
#' values. If FALSE, the function returns the filtered state vector
#' @param constrainParam booelan. If TRUE, a parameter constraining function is applied
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return either a list of the filtered output or a double of the log likelihood

KalmanFilter <- function(param, data, outLogLik, constrainParam, dgp1) {
  # Run the recursions
  if (!is.matrix(data)) data <- matrix(data, nc = 1)
  Output <- KalmanRecursions(
    paramVec = param, data = data, outLogLik = outLogLik, constrainParam = constrainParam,
    dgp1 = dgp1
  )
  # If the function only returns a likelihood value, extract this value from the output list
  if (outLogLik == TRUE) Output <- Output[[1]]
  return(Output)
}


#' @description Function that finds the ML estimates via numerical optimization using a range of different initial parameters
#' @param thetaMat Matrix with different initial values. One set of values per row
#' @param data matrix with the data. One column per observed variable
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return matrix with the log likelihood and ML estimates for every initial parameter vector (still constrained parameters)

GridParamOptim <- function(thetaMat, data, dgp1) {
  if (!is.matrix(data)) data <- matrix(data, nc = 1)
  
  # Wrap the optimization in a tryCatch in case the routine fails
  resultsOptim <- apply(thetaMat, 1, function(theta, data) {
    return(c(
      tryCatch(
        optim(theta,
          fn = KalmanRecursions, hessian = FALSE, method = "BFGS",
          control = list(maxit = 1e5, reltol = 1e-06),
          data = data, outLogLik = TRUE, constrainParam = TRUE, dgp1 = dgp1
        ),
        error = function(e) NA
      ),
      theta
    ))
  }, data = data) %>%
    t()
  # Collect parameters from the list of optim results
  if (!is.list(resultsOptim)) stop("Optimization failed for all initial values\n")
  optimUnfltrd <- sapply(resultsOptim, function(x, nparams) {
    if (is.list(x)) {
      outputVec <- c(-x$value, x$par)
    } else {
      outputVec <- rep(NA, nparams)
    }
  }, nparams = NCOL(thetaMat) + 1) %>%
    t()
  
  if (dgp1 == TRUE) {
    colnames(optimUnfltrd) <- c("ll", "phi_1", "phi_2", "SdEta", "SdU", "SdE", "SdEpsilon")
  } else {
    colnames(optimUnfltrd) <- c("ll", "SdEta")
  }

  # Sort according to the log likelihood values
  optimUnfltrd <- optimUnfltrd[order(-optimUnfltrd[, 1]), ]
  # Discard optimization runs where the routine failed
  optimMat <- optimUnfltrd[!is.na(optimUnfltrd[, 1]), ]
  return(optimMat)
}


#' @description Function that finds the ML estimates for vector of initial values
#' @param data matrix with the data. One column per observed variable
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return vector with the log likelihood and ML estimates

ParamOptim <- function(theta, data, dgp1) {
  if (!is.matrix(data)) data <- matrix(data, nc = 1)
  
  # Wrap the optimization in a tryCatch in case the routine fails
  optimResult <- c(
    tryCatch(optim(theta,
      fn = KalmanRecursions, hessian = FALSE, method = "BFGS",
      control = list(maxit = 1e5, reltol = 1e-06),
      data = data, outLogLik = TRUE, constrainParam = TRUE, dgp1 = dgp1
    ), error = function(e) NA)
  )
  if (!is.list(optimResult)) {
    outputVec <- rep(NA, 1 + length(theta))
  } else {
    if (optimResult$convergence != 0) warning("Convergence of the likelihood maximation likely not archieved \n")
    outputVec <- c(-optimResult$value, ParConstrain(optimResult$par, dgp1 = dgp1))
  }
  
  if (dgp1 == TRUE) {
    names(outputVec) <- c("ll", "phi_1", "phi_2", "SdEta", "SdU", "SdE", "SdEpsilon")
  } else {
    names(outputVec) <- c("ll", "SdEta")
  }
  return(outputVec)
}


#' @description Function that computes SE based on the Hessian matrix, corrected for the optimization parameter
#' constraining function
#' @param theta vector of ML estimates (still constrained)
#' @param hessianMat Hessian of theta
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return Vector of the estimated standard errors

SEfctn <- function(theta, hessianMat, dgp1) {
  correctionMat <- InvParConstrain_fctn(paramVec = theta, dgp1 = dgp1)
  CovMat <- correctionMat %*% Inverse(hessianMat) %*% Transp(correctionMat)
  SEVec <- sqrt(diag(CovMat))
  names(SEVec) <- names(theta)
  return(SEVec)
}


#' @description Function that reverses the ParConstrain_fctn to correct for the SE calculation
#' @param paramVec vector of constrained parameters (still constrained)
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return a diagonal matrix with the correction factors

InvParConstrain_fctn <- function(paramVec, dgp1) {
  correctionVec <- rep(1, length(paramVec))
  if (dgp1 == TRUE) {
    correctionVec[1] <- numDeriv::jacobian(function(x) {
      2 * x / (1 + abs(x))
    }, x = paramVec[1])
    # this function triggers a warning since it returns the correction factor for only variable (phi_2)
    # even though it takes two parameters as its input
    suppressWarnings({
    correctionVec[2] <- numDeriv::jacobian(function(x) {
      phi_1 <- x[1]
      phi_2 <- x[2]
      -sum(
        ((1 - abs(2 * phi_1 / (1 + abs(phi_1)))) * phi_2) / (1 + abs(phi_2)),
        2 * phi_1 / (1 + phi_1)
      )
    }, x = c(paramVec[1] / 2, paramVec[2]))
    })
    correctionVec[3:6] <- diag(numDeriv::jacobian(exp, paramVec[3:6]))
    # Build the final correction matrix
    correctionMat <- diag(correctionVec)
  } else {
    correctionVec[1] <- diag(numDeriv::jacobian(exp, paramVec[1]))
    # Build the final correction matrix
    correctionMat <- matrix(correctionVec)
  }
  return(correctionMat)
}
