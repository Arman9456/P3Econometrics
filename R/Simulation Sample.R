#' @description draws observations. Y is a noisy measurement of a persistent I(2) component and a
#' stationary AR(3) cycle
#' @param N number of drawn sample paths
#' @param BI Burn in period for the draws
#' @param BID Number of observations to be drawn
#' @return matrix with the drawn samples. One column per sample path.

GenSamples <- function() {
  # Coefficients cycle AR(3) equation & initial values which we choose s.t. the process is stationary
  phi_1 <- 0.1
  phi_2 <- 0.12
  phi_3 <- 0
  SdEpsilon <- .05
  SdEta <- .6
  SdU <- .01
  SdE <- .5
  paramVec <- c(phi_1, phi_2, phi_3, SdEpsilon, SdEta, SdU, SdE)

  N <- 1000 # number of simulations
  BI <- 500 # number of iterations for burn-in
  BID <- 1000 # number of observations

  # Generate a sample of y
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
    alpha <- rep(0.01, nPeriods)
    cyc <- rep(0, nPeriods)
    epsilon <- rnorm(nPeriods, sd = SdEpsilon)
    eta <- rnorm(nPeriods, sd = SdEta)
    u <- rnorm(nPeriods, sd = SdU)
    e <- rnorm(nPeriods, sd = SdE)
    # Construct the components
    for (i in 4:nPeriods) {
      alpha[i] <- alpha[i - 1] + u[i]
      mu[i] <- alpha[i] + mu[i - 1] + eta[i]
      cyc[i] <- phi_1 * cyc[i - 1] + phi_2 * cyc[i - 2] + phi_3 * cyc[i - 3] + e[i]
      y[i] <- mu[i] + cyc[i] + epsilon[i] 
    }
    # Discard the burn in period
    y_out <- y[-c(1:(3 + BI))]
    return(y_out)
  }, BI = BI, BID = BID, paramVec = paramVec)

  return(ymat)
}

# test = GenSamples()
