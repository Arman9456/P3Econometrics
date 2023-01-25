# Main code#
source("R/Init.R")
Init()

# Test the filter
yRandom <- GenSamples()
dataMat <- as.matrix(yRandom[, 1])


logLik <- KalmanFilter(param = rep(0, 6), data = dataMat, outLogLik = T, constrainParam = T)
Output <- KalmanFilter(param = rep(0, 6), data = dataMat, outLogLik = F, constrainParam = T)
# Test the parameter optimization
ParamOptim(theta = rep(0, 5), data = dataMat)
# Test the paramter grid search
N <- 5
thetaMat <- matrix(c(runif(N, 0, 1),
  runif(N, 0, 1),
  -log(runif(N, 0, .1)),
  -log(runif(N, 0, .1)),
  -log(runif(N, 0, 1.5))),
  nr = N)

Output <- GridParamOptim(thetaMat = thetaMat, data = dataMat)

# Issue with parameter identification

filterOutput <- KalmanFilter(param = Output[1,-1], data = dataMat, outLogLik = F, constrainParam = T)


plot(filterOutput$Z_tt[, 1] + filterOutput$Z_tt[, 3], type = "l")
lines(dataMat, col = "red")

# Test the bootstrap
y_star <- GenBootObs(data = dataMat, filterOutput = filterOutput)

d_mat <- BootstrapRoutine(B = 10, data = dataMat, filterOutput = filterOutput, theta = Output[1,-1])
plot(d_mat[,2] ~ d_mat[,1], type = "l")
