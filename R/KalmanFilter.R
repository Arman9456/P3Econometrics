#' @description Wrapper function for the Kalman filter
#' @param param vector of structural parameters
#' @param data matrix with the data. One column per observed variable
#' @param outLogLik boolean. If TRUE, the function returns only the log likelihood and no filtered
#' values. If FALSE, the function returns the filtered state vector
#' @param constrainParam booelan. If TRUE, a parameter constraining function is applied
#' @return either a list of the filtered output or a double of the log likelihood

KalmanFilter <- function(param, data, outLogLik, constrainParam) {
  # Apply parameter constraints if necessary
  if (constrainParam == TRUE) param <- ParConstrain(paramVec = param)
  # Run the recursions
  if (!is.matrix(data)) data <- matrix(data, nc = 1)
  Output <- KalmanRecursions(
    paramVec = param, data = data, outLogLik = outLogLik
  )
  # If the function only returns a likelihood value, extract this value from the output list
  if (outLogLik == TRUE) Output <- Output[[1]]
  return(Output)
}


#' @description Function that finds the ML estimates via numerical optimization with a range of initial parameters
#' @param thetaMat Matrix with different initial values. One set of values per row
#' @param data matrix with the data. One column per observed variable
#' @return matrix with the log likelihood and ML estimates for every initial parameter vector

GridParamOptim <- function(thetaMat, data) {
  if (!is.matrix(data)) data <- matrix(data, nc = 1)
  thetaMat <- t(apply(thetaMat, 1, ParConstrain))
  resultsOptim <- apply(thetaMat, 1, function(theta, data) {
    #return(c(
     # tryCatch(
        optim(theta,
        fn = KalmanRecursions, hessian = FALSE, method = "Nelder-Mead",
        control = list(maxit = 1e5, reltol = 1e-06),
        data = data, outLogLik = TRUE,
      )
      #, error = function(e) NA),
     # theta
  #  ))
  }, data = data) %>%
    t()
  # Collect parameters from the list of optim results
  if (!is.list(resultsOptim)) stop("Optimization failed for all initial values\n")
  optimUnfltrd <- sapply(resultsOptim, function(x, nparams) {
    if (is.list(x)) {
      outputVec <- c(-x$value, ParConstrain(x$par))
    } else {
      outputVec <- rep(NA, nparams)
    }
  }, nparams = NCOL(thetaMat) + 1) %>%
    t()
   colnames(optimUnfltrd) <- c("ll", "phi_1", "phi_2", "SdEta", "SdU", "SdE")

  # Sort according to the log likelihood values
  optimUnfltrd <- optimUnfltrd[order(-optimUnfltrd[, 1]), ]
  # Discard optimization runs where the routine failed
  optimMat <- optimUnfltrd[!is.na(optimUnfltrd[, 1]), ]
  return(optimMat)
}


#' @description Function that finds the ML estimates for one set of initial values
#' @param data matrix with the data. One column per observed variable
#' @return vector with the log likelihood and ML estimates

ParamOptim <- function(theta, data) {
  theta <- ParConstrain(paramVec = theta)
  if (!is.matrix(data)) data <- matrix(data, nc = 1)
  optimResult <- c(
    tryCatch(optim(theta,
      fn = KalmanRecursions, hessian = FALSE, method = "Nelder-Mead",
      control = list(maxit = 1e5, reltol = 1e-06),
      data = data, outLogLik = TRUE,
    ), error = function(e) NA),
    theta
  )
  if (!is.list(optimResult)) {
    outputVec <- rep(NA, 1 + length(theta))
  } else {
    if (optimResult$convergence != 0) warning("Convergence of the likelihood maximation likely not archieved \n")
    outputVec <- c(optimResult$value, ParConstrain(optimResult$par))
  }
  names(outputVec) <- c("ll", "phi_1", "phi_2", "SdEta", "SdU", "SdE")

  return(outputVec)
}


#' @description Function that applies parameter constraints
#' @param spec character indicating the model specification
#' @param par vector structural parameters
#' @return constrained parameter vector

ParConstrain <- function(paramVec) {
  constrPar <- paramVec
  # Constrain AR parameters of the cycle to result in a stable process
  phi_1 <- 2 * paramVec[1] / (1 + abs(paramVec[1]))
  phi_2 <- -(1 - abs(phi_1)) * paramVec[2] / (1 + abs(paramVec[2])) - abs(phi_1)
  constrPar[1] <- phi_1
  constrPar[2] <- phi_2
  # Constrain system St.Devs to be positive
  constrPar[3:5] <- exp(paramVec[3:5])
  return(constrPar)
}
