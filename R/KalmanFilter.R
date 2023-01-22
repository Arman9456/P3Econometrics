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
    paramVec = param, data = data, systemList = systemList, outLogLik = outLogLik
  )
  return(Output)
}

dataMat <- matrix(rnorm(10), nc = 2)
paramVec <- rep(1, 4)


#' @description Function that finds the ML estimates via numerical optimization with a range of initial parameters
#' @param thetaMat Matrix with different initial values. One set of values per row
#' @param data matrix with the data. One column per observed variable
#' @return matrix with the log likelihood and ML estimates for every initial parameter vector

thetaMat <- matrix(rnorm(20), nc = 2)

GridParamOptim <- function(thetaMat, data) {
  resultsOptim <- apply(thetaMat, 1, function(theta, data, systemList) {
    systemList <- SystemMat_fctn(paramVec = theta)
    return(c(
      tryCatch(optim(theta,
        fn = KalmanRecursions, hessian = FALSE, method = "BFGS",
        control = list(maxit = 100, reltol = 1e-06),
        data = data, systemList = systemList, outLogLik = TRUE,
      ), error = function(e) NA),
      theta
    ))
  }, data = data, systemList = systemList) %>%
    t()
  # Collect parameters from the list of optim results
  if (!is.list(resultsOptim)) stop("Optimization failed for all initial values\n")
  optimUnfltrd <- sapply(resultsOptim, function(x, nparams) {
    if (is.list(x)) {
      outputVec <- c(-x$value, x$par)
    } else {
      oututVec <- rep(NA, nparams)
    }
  }, nparams = NCOL(thetaMat) + 1) %>%
    t()
  colnames(optimUnfltrd) <- c("ll", "alpha", "beta")
  # Sort according to the log likelihood values
  optimUnfltrd <- optimUnfltrd[order(-optimUnfltrd[, 1]), ]
  # Discard optimization runs where the routine failed
  optimMat <- optimUnfltrd[!is.na(optimUnfltrd[, 1]), ]
  return(optimMat)
}


#' @description Function that finds the ML estimates for one set of initial values
#' @param data matrix with the data. One column per observed variable
#' @return vector with the log likelihood and ML estimates

theta <- rnorm(2)

ParamOptim <- function(theta, data) {
  systemList <- SystemMat_fctn(paramVec = theta)
  optimResult <- c(
    tryCatch(optim(theta,
      fn = KalmanRecursions, hessian = FALSE, method = "BFGS",
      control = list(maxit = 100, reltol = 1e-06),
      data = data, systemList = systemList, outLogLik = TRUE,
    ), error = function(e) NA),
    theta
  )
  if (!is.list(optimResult)) {
    outputVec <- rep(NA, 1 + length(theta))
  } else {
    if (optimResult$convergence != 0) warning("Convergence of the likelihood maximation likely not archieved \n")
    outputVec <- c(optimResult$value, optimResult$par)
  }
  names(outputVec) <- c("ll", "alpha", "beta")
  return(outputVec)
}


#' @description Function that applies parameter constraints
#' @param spec character indicating the model specification
#' @param par vector structural parameters
#' @return constrained parameter vector

ParConstrain <- function(paramVec) {
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

  # Just random values to check the filter
  A <- matrix(rnorm(25), nr = 5)
  B <- diag(5)
  # Measurement eq
  C <- matrix(rnorm(10), nr = 2, nc = 5)
  D <- diag(2)

  outputList <- list("A" = A, "B" = B, "C" = C, "D" = D)
  return(outputList)
}
