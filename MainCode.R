# Main code#
rm(list = ls())

source("R/Init.R")
Init()

# Test the filter
yRandom <- GenSamples(dgp1 = F, N = 1, BID = 500)
dataMat <- as.matrix(yRandom[, 1])

KalmanFilter(param = c(0), data = dataMat, outLogLik = T, constrainParam = T, dgp1 = F)

# Test the parameter optimization
ParamOptim(theta = rep(0, 1), data = dataMat, dgp1 = FALSE)
# Test the paramter grid search
N <- 10
thetaMat <- matrix(c(runif(N, 0, 1),
  runif(N, 0, 1),
  log(runif(N, .5, 2)), # eta
  log(runif(N, 0, .1)), # u
  log(runif(N, 0, 1)), # e
  log(runif(N, 0, 1.5))), # epsilon
  nr = N)

thetaMat <- matrix(c(log(runif(N, .5, 2))), # eta
                   nr = N)

Output <- GridParamOptim(thetaMat = thetaMat, data = dataMat, dgp1 = F)
apply(as.matrix(Output[, -1]), 1, ParConstrain, dgp1 = F) %>% t()

filterOutput <- KalmanFilter(param = Output[1, -1], data = dataMat, outLogLik = F, constrainParam = T, dgp1 = F)

plot(filterOutput$Z_tt[, 1] + filterOutput$Z_tt[, 3], type = "l")
plot(filterOutput$Z_tt[, 1], type = "l")
lines(dataMat, col = "red")

filterOutput$Z_tt[, 2]
# Test the bootstrap
y_star <- GenBootObs(data = dataMat, filterOutput = filterOutput)

d_mat <- BootstrapRoutine(B = 5, data = dataMat, filterOutput = filterOutput, theta = Output[1,-1])
plot(d_mat[,2] ~ d_mat[,1], type = "l")

#iterationOutput <- MonteCarloRoutine(dataVec = yRandom[, 1], nBoot = 10, dgp1 = T, CDFsupport = seq(-10, 10, .01))
MonteCarloOutput <- RunMonteCarlo(nSim = 5, nPeriods = 100, nBoot = 10, dgp1 = F)
