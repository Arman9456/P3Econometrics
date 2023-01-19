# Main code#
#install.packages()

library(matlib)
library(expm)

# DGP parameters #
sigma = 1
alpha = 0.75
beta = 0.99
omega = 1
rho_z = 0.90
theta = 6
sigma_z = 0.30
sigma_r = 0.20
kappa =( ((1-alpha)*(1-alpha*beta))/alpha )*( (omega+sigma)/(sigma*(omega+theta)))
phi_pi = 1.5
phi_x = 0.125
rho_r = 0.75




B=5 # nr of bootstrap iterations
y = # original sample to be used in bootstrap of dimension T * n_y, so n_y variables per time point
T = 100 # observations of y

n_z = # number of state space variables
n_y = 5# number of observed variables
C = #n_y x n_z matrix of coefficients for measurement equation
A = # n_z x n_z matrix of coefficients of state equation
Z = # (T-1) x n_z matrix of state variables over all time points (rows)
B = # n_z x n_w matrix of coefficients for shocks of state equation
D = # n_z x n_w matrix of coefficients for shocks of measurement equation


##### Monte carlo simulation #########  
for (i in 1:N){
#########Apply kalman filter and get a sample###############

  ####### apply bootstrap in this loop
}

eps_hat = matrix(0,nrow=(T-1),ncol=n_y)# initialize empty innovation residual matrix of  dimension (T-1) * n_y  

cov_eps_hat =  C_hat %*% P_t %*% t(C_hat) + D_hat %*% cov_w_hat %*% t(D_hat) # estimated covariance matrix at time t of dimension n_y x n_y 
eps_hat_c =  matrix(0,nrow=(T-1),ncol=n_y) # (T-1)*n_y matrix of centered residuals

sigma_epsilon_hat_sqrt_inv =  diag(inv(sqrtm(cov_eps_hat))) # list of dimension (T-1) x  (n_y x n_y) so T-1 matrices of dimension (n_y x n_y)

res_hat = sigma_epsilon_hat_sqrt_inv * eps_hat_c 

eps_hat_sum = t(rep(0,T)) # 1x n_y vector of sums over time

# fill up eps_hat and sum over time to get centered residuals
for(i in 1:T-1){
  eps_hat[i+1,] = y[i+1,] - C %*% Z_hat[...] # we start from index i+1 for eps_hat since we need to start at eps_hat_2
  eps_hat_sum = eps_hat_sum + eps_hat[i+1,] # for each time T, we add the 1xn_y vectors until we get a 1xn_y vector of the sum
  }

# get centered residuals and standardized innovations
for(i in 1:T-1){
  eps_hat_c[i+1,] = eps_hat[i+1,] - (1/(T-1))*eps_hat_sum # we make the centered residuals
  res_hat = sigma_epsilon_hat_sqrt_inv[i+1[... how to make list of matrices]] %*% t(eps_hat_c[i+1,]) # should be multiplication of n_y x n_y matrix with n_y x 1 eps_hat_c[i+1,] for each t
}


for (i in 1:B){
  star_indices = ceiling(runif(n=T-1,min=2,max=T)) #draw indices with replacement from original sample 
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