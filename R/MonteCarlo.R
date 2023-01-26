#' @description Executes the whole Monte Carlo simulation study
#' @param nSim number of Monte Carlo simulations
#' @param nPeriods length of the simulated time series
#' @param nBoot number of Bootstrap draws
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return Array. Each matrix has the support for the CDF as its first column and the values of the density d (eq 23) on all subsequent columns.
#' (one per estimated structural parameter)

RunMonteCarlo <- function(nSim, nPeriods, nBoot, dgp1) {
  yRandom <- GenSamples(N = nSim, BID = nPeriods, dgp1 = dgp1)
  xSeq <- seq(-10, 10, .01)
  monteCarloOutput <- apply(yRandom, 2, MonteCarloRoutine, nBoot = nBoot, dgp1 = dgp1, CDFsupport = xSeq)
  browser()
  # Bring the output into a proper form
  # theta mean
  thetaMeanVec <- apply(sapply(monteCarloOutput, function(x) x$theta), 1, mean)
  # theta se
  thetaSeVec <- apply(sapply(monteCarloOutput, function(x) x$thetaSE), 1, mean)
  # median theta CI length
  thetaCIlength <- apply(sapply(monteCarloOutput, function(x) x$thetaCIlength), 1, median)
  # median theta CI length
  thetaCI <- matrix(apply(sapply(monteCarloOutput, function(x) x$thetaCI), 1, median), nr = 2)
  colnames(thetaCI) <- names(thetaMeanVec)
  # theta star mean
  thetaStarMeanVec <- apply(sapply(monteCarloOutput, function(x) x$theta_star), 1, mean)
  monteCarloOutput[[2]][[2]]
  
  
  return(Output)
}


#' @description Runs one individual Monte Carlo simulation
#' @param dataVec simulated sample path
#' @param nBoot number of Bootstrap draws
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @param CDFsupport numerical vector. Values for which the CDF is evaluated
#' @return the test statistic d for each element of the parameter vector

MonteCarloRoutine <- function(dataVec, nBoot, dgp1, CDFsupport) {
  dataMat <- as.matrix(dataVec)
  # Run a parameter gridsearch for the original theta
  nGrid <- 50
  if (dgp1 == TRUE) {
    thetaMat <- matrix(
      c(
        runif(nGrid, 0, 1), # phi_1
        runif(nGrid, 0, 1), # phi_2
        log(runif(nGrid, .5, 2)), # eta
        log(runif(nGrid, 0, .1)), # u
        log(runif(nGrid, 0, 1)), # e
        log(runif(nGrid, 0, 1.5)) # epsilon
      ),
      nr = nGrid
    )
  } else {
    thetaMat <- matrix(c(runif(nGrid, 0, 1)),
      nr = nGrid
    )
  }
  Output <- GridParamOptim(thetaMat = thetaMat, data = dataMat, dgp1 = dgp1)
  theta <- Output[1, -1]
  # Calculate the Hessian matrix for the best fitting parameter vector
  Hessian <- optimHess(theta, KalmanRecursions, data = dataMat, outLogLik = TRUE, constrainParam = TRUE, dgp1 = dgp1)
  thetaSE <- SEfctn(theta = theta, hessianMat = Hessian, dgp1 = dgp1)
  # Calculate the 90% CI for theta
  thetaCIlength <- (1.64 * thetaSE) * 2
  thetaCImat <- rbind(thetaCIlength / 2 + theta,
                      thetaCIlength / 2 - theta)
  # Get the filter output
  filterOutput <- KalmanFilter(param = theta, data = dataMat, outLogLik = F, constrainParam = T, dgp1 = dgp1)
  bootstrapOutput <- BootstrapRoutine(B = nBoot, data = dataMat, filterOutput = filterOutput, theta = theta, CDFsupport = CDFsupport)
  # Construct the output list
  outputList <- list(
    "d_mat" = bootstrapOutput$d_statistic,
    "theta" = theta,
    "thetaSE" = thetaSE,
    "thetaCIlength" = thetaCIlength,
    "thetaCI" = thetaCImat,
    "theta_star" = bootstrapOutput$theta_star
  )
  return(outputList)
}


#' @description Generates the sample paths for the Monte Carlo simulation study
#' @param N number of drawn sample paths
#' @param BID Number of observations to be drawn
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return matrix with the drawn samples. One column per sample path.

GenSamples <- function(N, BID, dgp1) {
  set.seed(123)
  # Specify the parameters
  phi_1 <- 0.1
  phi_2 <- 0.12
  phi_3 <- 0
  SdEta <- 1
  SdU <- .01
  SdE <- .5
  SdEpsilon <- .05
  paramVec <- c(phi_1, phi_2, phi_3, SdEpsilon, SdEta, SdU, SdE)
  # Generate the samples
  ymat <- sapply(1:N, function(x, BI, BID, paramVec) {
    phi_1 <- paramVec[1]
    phi_2 <- paramVec[2]
    SdEpsilon <- paramVec[4]
    SdEta <- paramVec[5]
    SdU <- paramVec[6]
    SdE <- paramVec[7]
    # Create vectors
    nPeriods <- BI + BID + 3
    y <- rep(0, nPeriods)
    mu <- rep(0, nPeriods)
    alpha <- rep(0, nPeriods)
    cyc <- rep(0, nPeriods)
    epsilon <- rnorm(nPeriods, sd = SdEpsilon)
    eta <- rnorm(nPeriods, sd = SdEta)
    u <- rnorm(nPeriods, sd = SdU)
    e <- rnorm(nPeriods, sd = SdE)
    # Construct the components
    if (dgp1 == TRUE) {
      for (i in 4:nPeriods) {
        alpha[i] <- alpha[i - 1] + u[i]
        mu[i] <- alpha[i] + mu[i - 1] + eta[i]
        cyc[i] <- phi_1 * cyc[i - 1] + phi_2 * cyc[i - 2] + e[i]
        y[i] <- mu[i] + cyc[i] + epsilon[i]
      }
    } else {
      for (i in 4:nPeriods) {
        y[i] <- .01 + y[i - 1] + eta[i]
      }
    }
    # Discard the burn in period
    y_out <- y[-c(1:(3 + BI))]
    return(y_out)
  }, BI = BI, BID = BID, paramVec = paramVec)
  return(ymat)
}
