# Main code 
#install.packages()

library(matlib)
library(expm)

B=5 # nr of bootstrap iterations
y = # original sample to be used in bootstrap of dimension T * n_y, so n_y variables per time point
T = 100 # observations of y

n_z = # number of state space variables
n_y = # number of observed variables
#########apply kalman filter###############


eps_hat = # initialize empty innovation residual matrix of  dimension (T-1) * n_y  
cov_eps_hat =   # estimated covariance matrices
eps_hat_c = eps_hat - (1/(T-1))*sum(eps_hat) # vector of centered residuals

sigma_epsilon_hat_sqrt_inv = inv(sqrtm(cov_eps_hat))

res_hat = sigma_epsilon_hat_sqrt_inv * eps_hat_c # vector of standardized residuals from index 1 to T-1 which are supposed to be t=2 to T

for(i in 1:T-1){
  eps_hat[,i+1] = y[,i+1] - C %*% Z_hat[...] # we start from index i+1 for eps_hat since we need to start at eps_hat_2
}

for(i in 1:T-1){
  eps_hat[i+1] = y[i+1] - C %*% Z_hat[i] # we start from index i+1 for eps_hat since we need to start at eps_hat_2
}


for (i in 1:B){
  star_indices = ceiling(runif(n=T-1,min=1,max=T-1)) #draw indices with replacement from original sample (t=2,...t=T is equivalent to t=1,...,t=T-1 in our case)
  e_star = res_hat[star_indices] # draw with replacement from sample to optain bootstrap standardized innovations
  
  # initialize starting values###
  y_star = rep(NA,T) # create empty bootstrap y vector to fill up
  Z_star = rep(NA,)
  
  Z_star[1] = Z_hat[1]
  y_star[1] = y[1]
  ########################
  
  for(i in 1:(T-1)){
    Z_star[i+1] = A %*% Z_star[i] + K_gain * sqrtm(cov_eps_hat)
  }
  
  
}