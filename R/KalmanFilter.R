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
  # Set up the filter
  systemList <- SystemMat_fctn(paramVec = param)
  # Run the recursions
  Output <- KalmanRecursions(
    paramVec = param, data = data, systemList = systemList, outLogLik = outLogLik
  )
  # If the function only returns a likelihood value, extract this value from the output list
  if (outLogLik == TRUE) Output <- Output[[1]]
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
  colnames(optimUnfltrd) <- c("ll", "phi_1", "phi_2", "SdEpsilon", "SdEta", "SdU", "SdE")
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
  colnames(optimUnfltrd) <- c("ll", "phi_1", "phi_2", "SdEpsilon", "SdEta", "SdU", "SdE")
  return(outputVec)
}


#' @description Function that applies parameter constraints
#' @param spec character indicating the model specification
#' @param par vector structural parameters
#' @return constrained parameter vector

ParConstrain <- function(paramVec) {
  constrPar <- paramVec
  # Constrain AR parameters of the cycle to result in a stable process
  phi_1 <- 2 * par[1] / (1 + abs(par[1]))
  phi_2 <- -(1 - abs(phi_1)) * par[2] / (1 + abs(par[2])) - abs(phi_1)
  constrPar[6] <- phi_1
  constrPar[7] <- phi_2
  # Constrain system St.Devs to be positive
  constrPar[3:5] <- exp(-paramVec[3:5]) 
  return(constrPar)
}


#' @description Function constructs the MODEL SPECIFIC system matrices for the Kalman filter
#' @param spec character indicating the model specification
#' @param paramVec vector structural parameters
#' @return list with the matrices

SystemMat_fctn <- function(paramVec) {
  phi_1 <- paramVec[1]
  phi_2 <- paramvec[2]
  # Transition eq
  A <- rbind(c(1, 1, 0, 0),
             c(0, 1, 0, 0),
             c(0, 0, phi_1, phi_2),
             c(0, 0, 1, 0))
  B <- rbind(diag(3), rep(0, 3))
  # Measurement eq
  C <- matrix(c(1, 0, 1, rep(0, NCOL(A) - 3)), nr = 1)

  outputList <- list("A" = A, "B" = B, "C" = C)
  return(outputList)
}
