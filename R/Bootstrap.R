#' @description generates bootstrapped observations eq. (14) (Bootstrap algorithm step 3)
#' @param data matrix with the data. One column per observed variable
#' @param filterOutput output from the KalmanFilter function
#' @param dimObs number of variables in the observed vector
#' @param bootErr matrix of bootstrapped innovations. Output from the GenBootErr function
#' @return list with the bootstrapped state and observation vector

GenBootObs <- function(data, filterOutput, dimObs, bootErr) {
  K_gain <- filterOutput$K_t
  # Construct state vectors
  epsilonCenter_list <- apply(bootErr, 1, as.matrix)
  Z_hat_star_mat <- apply(K_gain, c(1, 2), function(x) {
    x %*% bootErr[, as.integer(substring(deparse(substitute(x)), 2))]
  })
  Z_hat_star_mat[1, ] <- filterOutput$Z_tt[1, 1]
  # Construct the observation vector
  identMat <- diag(dimObs)
  y_star_mat <- identMat %*% bootErr
  y_star_mat[1, ] <- data[1, ]
  # Construct the output
  outputList <- list("Z_hat_star" = Z_hat_star_mat, "y_star" = y_star_mat)
  return(outputList)
}


#' @description Function that draws a bootstrap sample from the standardized Kalman innovations of size T - 1 where T is the observational period
#' (Bootstrap algorithm step 2)
#' @param forecastErr matrix with the standardized one-step ahead forecasting error from the Kalman filter. One row per period
#' @return matrix with the bootstrapped sample

GenBootErr <- function(forecastErr) {
  nBoot <- NCOL(forecastErr) - 1
  e_star_mat <- apply(forecastErr, MARGIN = 2, sample, size = nBoot, replace = TRUE)
  return(Transp(e_star_mat))
}


#' @description Function that computes the test statistic
#' @param data matrix with the data. One column per observed variable
#' @param theta parameter vector of ML estimates
#' @param thetaBootMat matrix with the bootstrapped thetas (one row per bootstrap)
#' @return matrix with the bootstrapped sample

DistThetastarDiff <- function(data, theta, thetaBootMat) {
  nPeriods <- NROW(data)
  # Quantify the bootstrap deviations in the paramter vector
  W_star_mat <- sqrt(nPeriods) * (thetaBootMat - theta)
  # Approximate the empirical distribution

  # subtract normal

  # tbd
}
