# Set up the project
rm(list = ls())
source("R/Init.R")
Init()

# Run the Monte Carlo simulation for DGP 1
MonteCarloOutput_DGP1 <- RunMonteCarlo(nSim = 2000, nPeriods = 100, nBoot = 20, dgp1 = TRUE)


# Run the Monte Carlo simulation for DGP 2
MonteCarloOutput_DGP2 <- RunMonteCarlo(nSim = 2000, nPeriods = 100, nBoot = 20, dgp1 = FALSE)
