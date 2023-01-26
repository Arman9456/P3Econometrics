#' @description Function that runs the bootstrap routine
#' @param B number of bootstrap draws
#' @param data matrix with the data. One column per observed variable
#' @param filterOutput output from the KalmanFilter function
#' @param theta non bootstrap maximum likelihood parameter vector
#' @param CDFsupport numerical vector. Values for which the CDF is evaluated
#' @return the test statistic d for each element of the parameter vector

BootstrapRoutine <- function(B, data, filterOutput, theta, CDFsupport){
  # To prevent errors
  B <- max(B, 3)
  # Get some parameters
  nPeriods <- NROW(data)
  dgp1 <- filterOutput$DGP1
  # Compute the bootstrap ML estimates
  theta_star_matRaw <- sapply(1:B, function(x){
    y_star <- GenBootObs(data, filterOutput)
    theta_star <- ParamOptim(theta = theta, data = y_star, dgp1 = dgp1)[-1]
  }) %>%
    t()
  # In case the matrix is transposes when only one parameter is estimates
  if (NROW(theta_star_matRaw) == 1) theta_star_matRaw <- t(theta_star_matRaw) 
  theta_star <- theta_star_matRaw[!is.na(theta_star_matRaw[,1]),]
  if (NCOL(theta_star_matRaw) == 1) theta_star <- as.matrix(theta_star) 
  nColResults <- NROW(theta_star)
  
  # In case the optimization failed somewhere
  while (nColResults < B){
    theta_star_matRaw <- sapply(1:(B/3), function(x){
      y_star <- GenBootObs(data, filterOutput)
      theta_star <- ParamOptim(theta = theta, data = y_star, dgp1 = dgp1)[-1]
    }) %>%
      t()
    # In case the matrix is transposes when only one parameter is estimates
    if (NROW(theta_star_matRaw) == 1) theta_star_matRaw <- t(theta_star_matRaw)
    theta_star_sub <- theta_star_matRaw[!is.na(theta_star_matRaw[,1]),]
    if (NCOL(theta_star_matRaw) == 1) theta_star_sub <- as.matrix(theta_star_sub) 
    theta_star <- rbind(theta_star, theta_star_sub)
    nColResults <- NROW(theta_star)
  }
  theta_star <- as.matrix(theta_star[1:B,])
  theta_star_mean <- apply(theta_star, 2, mean)
  
  # Compute the distance to the non bootstrap theta
  W_T_star <- sqrt(nPeriods) * (theta_star - as.numeric(ParConstrain(theta, dgp1 = dgp1)))
  
  
  G_starRaw <- apply(W_T_star, 2, function(W, CDFsupport){
   G_star <- sapply(CDFsupport, function(x, W_star){
     # eq (17)
      value <- length(which(W_star <= x)) / length(W_star)
      return(c(x, value))
    }, W_star = W)
   return(t(G_star))
  }, CDFsupport = CDFsupport)
  # Bring the matrix into a proper format (first column is x, the remaining columns are the function values)
  G_star <- cbind(CDFsupport, G_starRaw[-c(1:length(CDFsupport)),])
  # Pull the values of the standard normal CDF
  cdfNorm <- pnorm(CDFsupport)
  # eq (23)
  V_hat <- cdfNorm * (1 - cdfNorm)
  d_mat <- apply(as.matrix(G_star[,-1]), 2, function(x){
    d <- sqrt(B) * V_hat^(-.5) * (x - cdfNorm)
  })
  d_output <- cbind(CDFsupport, d_mat)
  return(list("d_statistic" = d_output, "theta_star" = theta_star_mean))
}


#' @description generates bootstrapped observations eq. (14) (Bootstrap algorithm step 3)
#' @param data matrix with the data. One column per observed variable
#' @param filterOutput output from the KalmanFilter function
#' @return list with the bootstrapped state and observation vector

GenBootObs <- function(data, filterOutput) {
  dimObs <- NCOL(data)
  # Draw the errors
  e_hat_star <- GenBootErr(forecastErr = filterOutput$e_hat)
  C <- filterOutput$C_theta
  Z_hat_star <- filterOutput$Z_tt
  Sigma_sqrt <- filterOutput$Sigma_sqrt
  
  e_hat_stand <- rep(NA, NROW(e_hat_star))
  for (i in 1:length(e_hat_stand)){
    Sigma_sqrt_mat <- Sigma_sqrt[,,i]
    e_hat_val <- e_hat_star[i,]
    e_hat_stand[i] <- c(Sigma_sqrt_mat * e_hat_star[i,])
  }
  
  browser()
  
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
  set.seed(123)
  nBoot <- NROW(forecastErr) - 1
  e_star_mat <- apply(forecastErr, MARGIN = 2, sample, size = nBoot, replace = TRUE)
  return(e_star_mat)
}

