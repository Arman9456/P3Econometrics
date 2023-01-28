
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
  mean_sq <- sum((dist - m1)^2) / nElem
  mean_3 <- sum((dist - m1)^3) / nElem
  mean_4 <- sum((dist - m1)^4) / nElem
  b_1 <- (mean_3 / mean_sq^(3 / 2))^2
  b_2 <- (mean_4 / mean_sq^2)
  # Calculate the nCols-Value
  testStat <- nElem * b_1 / 6 + nElem * (b_2 - 3)^2 / 24
  pValue <- 1 - pchisq(testStat, df = 2)

  return(pValue)
}



DHtest <- function(dist) {
  if (!is.matrix(dist)) dist <- as.matrix(dist)
  nRows <- nrow(dist)
  nCols <- NCOL(dist)
  Var <- var(dist)
  S.d <- diag(Var)
  VarCorrec <- diag(S.d^(-0.5))
  C <- VarCorrec %*% Var %*% VarCorrec
  L <- diag((eigen(C)$values)^-0.5)
  H <- eigen(C)$vectors
  y.i <- H %*% L %*% Transp(H) %*% VarCorrec %*% Transp(dist - sapply(dist, mean))
  B.1 <- apply(y.i, 1, function(x) {
    skew(x, "moments")
  })
  B.2 <- apply(y.i, 1, function(x) {
    kurt(x, "moments")
  })
  del <- (nRows - 3) * (nRows + 1) * (nRows^2 + (15 * nRows) - 4)
  a <- ((nRows - 2) * (nRows + 5) * (nRows + 7) * (nRows^2 + (27 * nRows) - 70)) / (6 *
    del)
  c <- ((nRows - 7) * (nRows + 5) * (nRows + 7) * (nRows^2 + (2 * nRows) - 5)) / (6 *
    del)
  k <- ((nRows + 5) * (nRows + 7) * (nRows^3 + 37 * nRows^2 + (11 * nRows) - 313)) / (12 *
    del)
  alpha <- a + B.1^2 * c
  chi <- (B.2 - 1 - B.1^2) * 2 * k
  Z.2 <- (((chi / (2 * alpha))^(1 / 3)) - 1 + (1 / (9 * alpha))) *
    ((9 * alpha)^(0.5))
  del <- (nRows - 3) * (nRows + 1) * (nRows^2 + (15 * nRows) - 4)
  a <- ((nRows - 2) * (nRows + 5) * (nRows + 7) * (nRows^2 + (27 * nRows) - 70)) / (6 *
    del)
  c <- ((nRows - 7) * (nRows + 5) * (nRows + 7) * (nRows^2 + (2 * nRows) - 5)) / (6 *
    del)
  k <- ((nRows + 5) * (nRows + 7) * (nRows^3 + 37 * nRows^2 + (11 * nRows) - 313)) / (12 *
    del)
  alpha <- a + B.1^2 * c
  chi <- (B.2 - 1 - B.1^2) * 2 * k
  Z.2 <- (((chi / (2 * alpha))^(1 / 3)) - 1 + (1 / (9 * alpha))) *
    ((9 * alpha)^(0.5))
  beta <- (3 * (nRows^2 + (27 * nRows) - 70) * (nRows + 1) * (nRows + 3)) / ((nRows -
    2) * (nRows + 5) * (nRows + 7) * (nRows + 9))
  w2 <- -1 + ((2 * (beta - 1))^0.5)
  del <- 1 / ((log(sqrt(w2)))^0.5)
  y <- B.1 * ((((w2 - 1) / 2) * (((nRows + 1) * (nRows + 3)) / (6 * (nRows -
    2))))^0.5)
  Z.1 <- del * (log(y + (y^2 + 1)^0.5))
  E <- t(Z.1) %*% Z.1 + t(Z.2) %*% Z.2
  pValue <- pchisq(E, 2 * nCols, lower.tail = FALSE)

  return(pValue)
}
