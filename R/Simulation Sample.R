#### Constructing Monte Carlo Samples ######

GenSamples <- function() {
# Coefficients cycle AR(3) equation & initial values which we choose s.t. the process is stationary
c1= 0.1
c2 = 0.12
c3 = 0.7

N = 1000 # number of simulations
BI = 1000 # number of iterations for burn-in
BID = 500 # number of observations to discard for burn-in
ymat = matrix(0, nrow = BID, ncol = N) # create matrix of dimension 500 x 1000, one column of 500 observations for each simulation

for(k in 1:N){
# Create vectors for burn in
y = rep(0,BI+1)
mu = rep(0,BI+1)
alpha = rep(0,BI+1)
cyc = rep(0,BI+3) # runs from cyc[4] to cyc[BI+3]
# initial values cycle equation
cyc[1] = rnorm(1,mean=0,sd=1)
cyc[2] = rnorm(1,mean=0,sd=1)
cyc[3] = rnorm(1,mean=0,sd=1)

# Initial values drift, trend, observation
y[1] = rnorm(1,mean=0,sd=1)
mu[1] = rnorm(1,mean=0,sd=1)
alpha[1] = rnorm(1,mean=0,sd=1)

  # Generate a sample of 500 observations of y
  for(i in 1:BI){
    alpha[i+1] = alpha[i] + rnorm(1,mean=0,sd=1)
    mu[i+1] = alpha[i+1] + mu[i] + rnorm(1,mean=0,sd=1)
    cyc[i+3] = c1*cyc[i+2] + c2*cyc[i+1] + c3*cyc[i] + rnorm(1,mean=0,sd=1)
    y[i+1] = mu[i+1] + cyc[i+4] + rnorm(1,mean=0,sd=1)
  }

# y values after discarding first 500 observations due to burn-in period
y = y[(BID+1):BI]
ymat[,k] = y # for each simulation we get 500 observations
}

return(ymat)
}

#test = GenSamples()
