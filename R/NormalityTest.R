#' @description Executes the univariate Jarque and Bera test for normality
#' @param dist vector with the distribution to test
#' @return p-Value of the test

JBtest <- function(dist) {
  if (!is.numeric(dist)) {
    if (is.matrix(dist) & NCOL(dist) == 1) {
      dist <- dist[, 1]
    } else {
      stop("Input distribution must be a matrix with one column or a numeric vector \nRows")
    }
  }
  nElem <- length(dist)
  # Compute the moments
  mean <- mean(dist)
  mean_sq <- sum((dist - mean)^2) / nElem
  mean_3 <- sum((dist - mean)^3) / nElem
  mean_4 <- sum((dist - mean)^4) / nElem
  b_1 <- (mean_3 / mean_sq^(3 / 2))^2
  b_2 <- (mean_4 / mean_sq^2)
  # Calculate the nCols-Value
  testStat <- nElem * b_1 / 6 + nElem * (b_2 - 3)^2 / 24
  pValue <- 1 - pchisq(testStat, df = 2)
  if (is.na(pValue) | is.infinite(pValue)) {
    pValue <- 0
  }
  return(pValue)
}