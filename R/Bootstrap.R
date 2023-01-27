#' @description Function that runs the bootstrap routine
#' @param B number of bootstrap draws
#' @param data matrix with the data. One column per observed variable
#' @param filterOutput output from the KalmanFilter function based on the sample theta
#' @param theta sample maximum likelihood parameter vector (still unconstrained)
#' @param CDFsupport numerical vector. Values for which the CDF is evaluated
#' @return the test statistic d for each element of the parameter vector

BootstrapRoutine <- function(B, data, filterOutput, theta, CDFsupport) {
  # To prevent errors
  B <- max(B, 3)
  # Get some parameters
  nPeriods <- NROW(data)
  dgp1 <- filterOutput$DGP1
  # Generate the bootstrap samples and the associated parameter estimates. The sample theta serves as
  # the initial values for the parameter optimization
  theta_star_matRaw <- sapply(1:B, function(x) {
    y_star <- GenBootObs(data, filterOutput)
    theta_star <- ParamOptim(theta = theta, data = y_star, dgp1 = dgp1)[-1]
  }) %>%
    t()
  # In case only one parameter is estimated, this matrix has to be transposed
  if (NROW(theta_star_matRaw) == 1) theta_star_matRaw <- t(theta_star_matRaw)
  theta_star <- theta_star_matRaw[!is.na(theta_star_matRaw[, 1]), ]
  if (NCOL(theta_star_matRaw) == 1) theta_star <- as.matrix(theta_star)
  nColResults <- NROW(theta_star)

  # In case the optimization failed somewhere, the routine is repeated until we have B theta_star estimates
  while (nColResults < B) {
    theta_star_matRaw <- sapply(1:(B / 3), function(x){
      y_star <- GenBootObs(data, filterOutput)
      theta_star <- ParamOptim(theta = theta, data = y_star, dgp1 = dgp1)[-1]
    }) %>%
      t()
    # In case only one parameter is estimated, this matrix has to be transposed
    if (NROW(theta_star_matRaw) == 1) theta_star_matRaw <- t(theta_star_matRaw)
    theta_star_sub <- theta_star_matRaw[!is.na(theta_star_matRaw[, 1]), ]
    if (NCOL(theta_star_matRaw) == 1) theta_star_sub <- as.matrix(theta_star_sub)
    theta_star <- rbind(theta_star, theta_star_sub)
    nColResults <- NROW(theta_star)
  }
  # Collect the results
  theta_star <- as.matrix(theta_star[1:B, ])
  theta_star_mean <- apply(theta_star, 2, mean)
  thetaML <- as.numeric(ParConstrain(theta, dgp1 = dgp1))

  # Compute the test statistic for the bootstrap CI
  Q_stat <- BootstrapQstat(theta = thetaML, theta_star_mat = theta_star, nPeriods = nPeriods)

  # Compute the standardized distance from sample to bootstrap theta
  W_T_star <- sqrt(nPeriods) * (theta_star - thetaML)

  # Compute the empirical distribution (G_star) of the difference between the bootstrap and the sample theta
  G_starRaw <- apply(W_T_star, 2, function(W, CDFsupport) {
    G_star <- sapply(CDFsupport, function(x, W_star) {
      # eq (17)
      value <- length(which(W_star <= x)) / length(W_star)
      return(c(x, value))
    }, W_star = W)
    return(t(G_star))
  }, CDFsupport = CDFsupport)
  # Bring the matrix into a proper format (first column is x, the remaining columns are the function values)
  G_star <- cbind(CDFsupport, G_starRaw[-c(1:length(CDFsupport)), ])

  # Pull the values of the standard normal CDF
  cdfNorm <- pnorm(CDFsupport)

  # Construct the final test distribution d eq (23)
  V_hat <- cdfNorm * (1 - cdfNorm)
  d_mat <- apply(as.matrix(G_star[, -1]), 2, function(x) {
    d <- sqrt(B) * V_hat^(-.5) * (x - cdfNorm)
  })
  d_output <- cbind(CDFsupport, d_mat)

  return(list(
    "d_statistic" = d_output,
    "theta_star" = theta_star_mean,
    "Qstat" = Q_stat
  ))
}


#' @description generates bootstrapped observations eq. (14) (Bootstrap algorithm step 3)
#' @param data matrix with the data. One column per observed variable
#' @param filterOutput output from the KalmanFilter function
#' @return matrix with the bootstrapped observation vector

GenBootObs <- function(data, filterOutput) {
  dimObs <- NCOL(data)
  # Draw the errors
  e_hat_star <- GenBootErr(forecastErr = filterOutput$e_hat)
  C <- filterOutput$C_theta
  Z_hat_star <- filterOutput$Z_tt
  Sigma_sqrt <- filterOutput$Sigma_sqrt
  # Standardize the errors
  e_hat_stand <- sapply(1:NROW(e_hat_star), function(x, Sigma, e_hat_star) {
    e_hat_star[x, ] %*% Sigma[, , x]
  }, Sigma = Sigma_sqrt, e_hat_star = e_hat_star)

  # Construct the observation vector
  y_star_mat <- apply(Z_hat_star, 1, function(x) C %*% as.matrix(x))[-1] + e_hat_stand
  y_star_final <- rbind(data[1, ], as.matrix(y_star_mat))

  return(y_star_final)
}


#' @description Function that draws a bootstrap sample from the standardized Kalman innovations of size T - 1 where T is the observational period
#' (Bootstrap algorithm step 2)
#' @param forecastErr matrix with the standardized one-step ahead forecasting error from the Kalman filter. One row per period
#' @return matrix with the bootstrapped sample

GenBootErr <- function(forecastErr) {
  nBoot <- NROW(forecastErr) - 1
  e_star_mat <- apply(forecastErr, MARGIN = 2, sample, size = nBoot, replace = TRUE)
  return(e_star_mat)
}


#' @description Function constructing the test statistic for the Equal-tailed Percentile CI - Interval
#' @param theta sample maximum likelihood parameter vector
#' @param theta_star_mat matrix with the bootstrap theta vectors
#' @param nPeriods number of time periods considered
#' @return the test statistic

BootstrapQstat <- function(theta, theta_star_mat, nPeriods) {
  theta_star_sd <- apply(theta_star_mat, 2, sd)
  theta_star_mean <- apply(theta_star_mat, 2, mean)
  Q_stat <- sqrt(nPeriods) * (theta_star_mean - theta) / theta_star_sd
  return(Q_stat)
}
