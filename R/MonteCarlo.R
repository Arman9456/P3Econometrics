#' @description Executes the whole Monte Carlo simulation study
#' @param nSim number of Monte Carlo simulations
#' @param nPeriods length of the simulated time series
#' @param nBoot number of Bootstrap draws
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return Array. Each matrix has the support for the CDF as its first column and the values of the density d (eq 23) on all subsequent columns.
#' (one per estimated structural parameter)

RunMonteCarlo <- function(nSim, nPeriods, nBoot, dgp1) {
  # Generate all the sample paths
  yRandom <- GenSamples(N = nSim, BID = nPeriods, dgp1 = dgp1)
  # Support of the empirical distribution functions
  xSeq <- seq(-20, 20, .01)
  # Run the individual simulations
  monteCarloOutput <- apply(yRandom, 2, MonteCarloRoutine, nBoot = nBoot, dgp1 = dgp1, CDFsupport = xSeq)
  # Save the output
  fileName <- ifelse(dgp1 == TRUE, "DGP_1", "DGP_2")
  write_rds(monteCarloOutput, file = paste0(getwd(), "/Output_", fileName, ".rds"))

  # # Bring the output into a proper form
  # # theta mean
  thetaMeanMat <- sapply(monteCarloOutput, function(x) x$theta)
  if (!is.matrix(thetaMeanMat)) thetaMeanMat <- matrix(thetaMeanMat, nr = 1)
  thetaMeanVec <- apply(thetaMeanMat, 1, mean)
  # theta se
  thetaSeMat <- sapply(monteCarloOutput, function(x) x$thetaSE)
  if (!is.matrix(thetaSeMat)) thetaSeMat <- matrix(thetaSeMat, nr = 1)
  thetaSeVec <- apply(thetaSeMat, 1, mean, na.rm = TRUE)
  # median theta CI length
  thetaCIlengtMat <- sapply(monteCarloOutput, function(x) x$thetaCIlength)
  if (!is.matrix(thetaCIlengtMat)) thetaCIlengtMat <- matrix(thetaCIlengtMat, nr = 1)
  thetaCIlength <- apply(thetaCIlengtMat, 1, median, na.rm = T)
  # median theta CI length
  thetaCI <- matrix(apply(sapply(monteCarloOutput, function(x) x$thetaCI), 1, median, na.rm = TRUE), nr = 2)
  colnames(thetaCI) <- names(thetaMeanVec)
  # theta star mean
  thetaStarMeanMat <- sapply(monteCarloOutput, function(x) x$theta_star)
  if (!is.matrix(thetaStarMeanMat)) thetaStarMeanMat <- matrix(thetaStarMeanMat, nr = 1)
  thetaStarMeanVec <- apply(thetaStarMeanMat, 1, mean)

  # Test the final test density for normality
  # Doornik Hansen test
  DHtestVec <- sapply(monteCarloOutput, function(x) {
    d_df <- as.data.frame(x$d_mat[, -1])
    testResult <- DH.test(d_df, Y.names = NULL)
    pValue <- testResult$multi[, 3]
    rejection <- ifelse(pValue <= .05, 0, 1)
    return(rejection)
  })
  DHfrequency <- mean(DHtestVec)
  # Jarque-Bera test
  JBtestMat <- sapply(monteCarloOutput, function(x) {
    d_mat <- as.matrix(x$d_mat[, -1])
    rejectionVec <- apply(d_mat, 2, function(dVec) {
      pValue <- JBtest(dVec)
      rejection <- ifelse(pValue <= .05, 0, 1)
      return(rejection)
    })
    return(rejectionVec)
  })
  JBfrequency <- apply(as.matrix(JBtestMat), 1, mean)
  # Shapiro-Wilk
  SWtestMat <- sapply(monteCarloOutput, function(x) {
    d_mat <- as.matrix(x$d_mat[, -1])
    rejectionVec <- apply(d_mat, 2, function(dVec) {
      pValue <- shapiro.test(dVec[-1])$p.value
      rejection <- ifelse(pValue <= .05, 0, 1)
      return(rejection)
    })
    return(rejectionVec)
  })
  SWfrequency <- apply(as.matrix(SWtestMat), 1, mean, na.rm = T)

  # Construct the Bootstrap Percentile CI
  Q_vec <- sapply(monteCarloOutput, function(x) x$Qstat)
  theta_star_CI <- BootCI(Q = Q_vec, mean = thetaMeanVec, se = thetaSeVec, alpha = .9, nPeriods = nPeriods)
  # Calculate the frequency with which the bootstrap theta falls in between the interval boundaries
  coverageOutputMat <- apply(thetaStarMeanMat, 2, function(x, CI_mat) {
    tempOutput <- sapply(1:length(x), function(i, param, CI_mat) {
      if (CI_mat[2, i] <= param[i] & param[i] <= CI_mat[1, i]) {
        indicator <- 1
      } else {
        indicator <- 0
      }
    }, param = x, CI_mat = CI_mat)
    names(tempOutput) <- names(x)
    return(tempOutput)
  }, CI_mat = theta_star_CI)
  if (!is.matrix(coverageOutputMat)) coverageOutputMat <- matrix(coverageOutputMat, nr = 1)
  coverageOutputVec <- apply(coverageOutputMat, 1, mean)

  # Construct the output
  Output <- list(
    "thetaMean" = thetaMeanVec,
    "thetaSE" = thetaSeVec,
    "thetaCIlength" = thetaCIlength,
    "thetaCI" = thetaCI,
    "thetaStarMean" = thetaStarMeanVec,
    "DHfrequency" = DHfrequency,
    "JBfrequency" = JBfrequency,
    "SWfrequency" = SWfrequency,
    "PercCoverage" = coverageOutputVec
  )
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

  # Run a quick grid search for the original parameter vector theta
  nGrid <- 50
  if (dgp1 == FALSE) {
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
    thetaMat <- matrix(c(runif(nGrid, 0, 1)), # eta
      nr = nGrid
    )
  }
  # Run the grid search
  Output <- GridParamOptim(thetaMat = thetaMat, data = dataMat, dgp1 = dgp1)
  thetaConst <- Output[1, -1]
  # Note that the routine outputs parameters that are still unconstrained
  theta <- ParConstrain(thetaConst, dgp1 = dgp1)

  # Calculate the Hessian matrix for the best fitting parameter vector
  Hessian <- optimHess(thetaConst, KalmanRecursions, data = dataMat, outLogLik = TRUE, constrainParam = TRUE, dgp1 = dgp1)

  if (is.singular.matrix(Hessian)) {
    thetaSE <- rep(0, length(theta))
  } else {
    # Based on the Hessian, compute the standard errors of the sample theta
    thetaSE <- SEfctn(theta = thetaConst, hessianMat = Hessian, dgp1 = dgp1)
  }

  # Calculate the 90% CI for theta
  thetaCIlength <- (1.64 * thetaSE) * 2
  thetaCImat <- rbind(
    theta + thetaCIlength / 2,
    theta - thetaCIlength / 2
  )
  # Get the filter output
  filterOutput <- KalmanFilter(param = thetaConst, data = dataMat, outLogLik = F, constrainParam = T, dgp1 = dgp1)

  # Based on the sample theta and the filter output, run the bootstrap
  bootstrapOutput <- BootstrapRoutine(B = nBoot, data = dataMat, filterOutput = filterOutput, theta = thetaConst, CDFsupport = CDFsupport)

  # Construct the output list
  outputList <- list(
    "d_mat" = bootstrapOutput$d_statistic,
    "theta" = theta,
    "thetaSE" = thetaSE,
    "thetaCIlength" = thetaCIlength,
    "thetaCI" = thetaCImat,
    "theta_star" = bootstrapOutput$theta_star,
    "Qstat" = bootstrapOutput$Qstat
  )
  return(outputList)
}


#' @description Generates the sample paths for the Monte Carlo simulation study
#' @param N number of drawn sample paths
#' @param BID Number of observations to be drawn
#' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
#' @return matrix with the drawn samples. One column per sample path.

GenSamples <- function(N, BID, dgp1) {
  set.seed(2)
  # Specify the parameters
  phi_1 <- 0.1
  phi_2 <- 0.12
  phi_3 <- 0.7
  SdEta <- 1
  SdU <- .01
  SdE <- .5
  SdEpsilon <- .05
  paramVec <- c(phi_1, phi_2, phi_3, SdEpsilon, SdEta, SdU, SdE)
  BI <- 500

  # Generate the samples
  ymat <- sapply(1:N, function(x, BI, BID, paramVec) {
    phi_1 <- paramVec[1]
    phi_2 <- paramVec[2]
    phi_3 <- paramVec[3]
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
    if (dgp1 == FALSE) {
      for (i in 4:nPeriods) {
        alpha[i] <- alpha[i - 1] + u[i]
        mu[i] <- alpha[i] + mu[i - 1] + eta[i]
        cyc[i] <- phi_1 * cyc[i - 1] + phi_2 * cyc[i - 2] + phi_3 * cyc[i - 3] + e[i]
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


#' @description Constructs the bootstrap CI with which the boot coverage may be computed
#' @param Q vector with the bootstrap test statistic
#' @param mean vector with the mean of the sample theta
#' @param se vector with the SE of the sample theta
#' @param alpha confidence level
#' @param nPeriods number of time periods observed
#' @return matrix with two rows for the lower and upper limits

BootCI <- function(Q, mean, se, alpha, nPeriods) {
  c_star <- quantile(Q, prob = c((1 - alpha) / 2, 1 - (1 - alpha) / 2), na.rm = TRUE)
  lowerCI <- mean - c_star[2] * se / sqrt(nPeriods)
  upperCI <- mean - c_star[1] * se / sqrt(nPeriods)
  return(rbind(upperCI, lowerCI))
}
