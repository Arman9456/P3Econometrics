#' @description Function that runs the bootstrap routine
#' @param B number of bootstrap draws
#' @param data matrix with the data. One column per observed variable
#' @param filterOutput output from the KalmanFilter function
#' @param theta non bootstrap maximum likelihood parameter vector
#' @return the test statistic d for each element of the parameter vector

BootstrapRoutine <- function(B, data, filterOutput, theta){
  nPeriods <- NROW(data)
  # Compute the bootstrap ML estimates
  theta_star_matRaw <- sapply(1:B, function(x){
    y_star <- GenBootObs(data, filterOutput)
    theta_star <- ParamOptim(theta = rep(0, 5), data = y_star)[-1]
  }) %>%
    t()
  theta_star <- theta_star_matRaw[!is.na(theta_star_matRaw[,1]),]
  # In case the optimization failed somewhere
  while (NROW(theta_star) < B){
    theta_star_matRaw <- sapply(1:(B/3), function(x){
      y_star <- GenBootObs(data, filterOutput)
      theta_star <- ParamOptim(theta = rep(0, 5), data = y_star)[-1]
    }) %>%
      t()
    theta_star <- rbind(theta_star, theta_star_matRaw[!is.na(theta_star_matRaw[,1]),])
  }
  theta_star <- theta_star[1:B,]
  # Compute the distance to the non bootstrap theta
  W_T_star <- sqrt(nPeriods) * (theta_star - theta)
  # Approximate the EDF G_star with support from -5 to 5 (range is arbitrary. Still tbd)
  xSeq <- seq(-10, 10, .01)
  G_starRaw <- apply(W_T_star, 2, function(W, xSeq){
   G_star <- sapply(xSeq, function(x, W_star){
     # eq (17)
      value <- length(which(W_star <= x)) / length(W_star)
      return(c(x, value))
    }, W_star = W)
   return(t(G_star))
  }, xSeq = xSeq)
  # Bring the matrix into a proper format (first column is x, the remaining columns are the function values)
  G_star <- cbind(xSeq, G_starRaw[-c(1:length(xSeq)),])
  # plot(G_star[,6] ~ G_star[,1], type = "l")
  
  # Everything after here is inefficient to compute for every bootstrap iteration individually.
  # Better to end function here, compute the rest once and apply to all monte carlo sample draws
  
  # Pull the values of the standard normal CDF
  cdfNorm <- pnorm(xSeq)
  # eq (23)
  V_hat <- cdfNorm * (1 - cdfNorm)
  d_mat <- apply(G_star[,-1], 2, function(x){
    d <- sqrt(B) * V_hat^(-.5) * (x - cdfNorm)
  })
  d_output <- cbind(xSeq, d_mat)
  #plot(d_mat[, 3], type = "l")
  return(d_output)
}


#' @description generates bootstrapped observations eq. (14) (Bootstrap algorithm step 3)
#' @param data matrix with the data. One column per observed variable
#' @param filterOutput output from the KalmanFilter function
#' @return list with the bootstrapped state and observation vector

GenBootObs <- function(data, filterOutput) {
  dimObs <- NCOL(data)
  # Draw the errors
  e_hat_star <- GenBootErr(forecastErr = filterOutput$e_hat)
  # Construct the observation vector
  y_star_mat <- e_hat_star
  y_star_mat[1, ] <- data[1, ]
  return(y_star_mat)
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

