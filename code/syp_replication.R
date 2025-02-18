##### Writing the code from scratch.
# Only necessary library is quantreg
# Clear everything
rm(list = ls())

# Load library
library(quantreg)

# Fix seed
set.seed(8675309)

################################################################################
## Data Generating Process Parameters
N <- 38
T <- 535
scale_factor_variance <- 0
error_variance <- 1
################################################################################
## DGP
# alpha_i drawn independently from N(0, 1)
alphas = matrix(rnorm(N), nrow = 1)
# beta_t drawn independently from N(0, 1)
betas = matrix(rnorm(T), nrow = 1)
# x_t drawn independently from N(0, 1)
xs = matrix(scale_factor_variance * rnorm(T), nrow = 1)
# v_it drawn independently from N(0, 1)
vs = matrix(rnorm(T * N), nrow = T, ncol = N)

# Gamma is exponential of x, want a TxN matrix for elementwise multiplication 
gammas = exp(xs)
gamma_build = matrix(gammas, nrow = T, ncol = N, byrow = TRUE)

# Combine to create observed variable y_it
y_it = t(betas) %*% alphas + gamma_build * vs
################################################################################
# From here, we have to get a starting value for the algorithm. We'll use PCA  
# for this. There are several potential rotations, I used errors-in-variables.
# It is also worth considering the rotation used in Ando and Bai (2020). 

## Here choose quantile of interest and number of factors
qtl <- 0.25
k_tau <- 1
################################################################################
## PCA
pca_mat = y_it %*% t(y_it)
f_tilde = as.matrix(eigen(pca_mat)$vectors[,1:k_tau] * sqrt(T))
l_tilde = t(y_it) %*% f_tilde / T 

## Errors in variables rotation
l_1 = l_tilde[1:k_tau,]
f_eiv = f_tilde %*% t(l_1)
l_eiv = l_tilde %*% solve(l_1)

## Ando and Bai (2020) rotation
# Their initialization procedure includes known regressors, so will look 
# slightly different. They take PCA F, and betas on a regression of F on Z for L
# Z <- AY-XB
# VEC <- eigen(Z%*%t(Z))$vectors; 
# F <- sqrt(N)*(VEC)[,1:R]; 
# L <- t(solve(t(F)%*%F)%*%t(F)%*%Z);
# this is equivalent to just PCA

################################################################################
## Estimation Parameters
diff = Inf
tol = 1e-3
iter = 0
maxiter = 100

# Use Ando Bai Rotation?
ab_rot = TRUE

# 
################################################################################
## Specify starting values
# Add intercept initial guesses
f_int = matrix(0, nrow = T)
l_int = matrix(0, nrow = N)
if (ab_rot) {
  f_prev = cbind(f_int, f_tilde)
  l_prev = cbind(l_int, l_tilde)
} else {
  f_prev = cbind(f_int, f_eiv)
  l_prev = cbind(l_int, l_eiv)
}

f_step = f_prev
l_step = l_prev
################################################################################
## Utility function for calculating the value of the objective function
get_check <- function(data, estimate, quantile) {
  residual = data - estimate
  check = residual * (quantile - 1 * (residual <= 0))
  return(sum(check))
}
## Utility function for rotating estimates according to AB
rotate_ab <- function(f_est, l_est, T_) {
  # Rotation step
  qr_f = qr(f_est)
  rf_tau = qr.R(qr_f)
  lam_mid = rf_tau %*% t(l_est)
  qr_l = qr(lam_mid)
  # Get rotated estimates
  f_rot = sqrt(T_) * qr.Q(qr_f) %*% qr.Q(qr_l)
  l_rot = t(qr.R(qr_l))
  
  return(list(f_est, l_est))
}
################################################################################
## Iteration

if (ab_rot) {
  while (diff > tol) {
    iter = iter + 1
    
    # AB update lambda first given F 
    # Update loadings
    for (i in 1:N) {
      l_step[i, ] = rq(y_it[, i] ~ f_step[, 2], tau = qtl)[[1]]
    }
    # Update factors
    for (j in 1:T) {
      f_step[j, ] = rq(y_it[j, ] ~ l_step[, 2], tau = qtl)[[1]]
    }
    
    # Rotation procedure
    rot_est = rotate_ab(f_step, l_step, T)
    
    f_step = rot_est[[1]]
    l_step = rot_est[[2]]
    
    # Check convergence
    # We first want the common component
    fl_est_prev = f_prev[, 2] %*% t(l_prev[, 2])
    fl_est = f_step[, 2] %*% t(l_step[, 2])
    
    # Get value of objective function
    obj_val = get_check(y_it, fl_est, qtl)
    
    # Convergence criterion will be euclidean norm between common components in
    # successive iterations
    diff = sqrt(sum((fl_est - fl_est_prev)^2))
  
    cat("Iteration ", iter, 
        ", Objective Function Value: ", obj_val, 
        ", Convergence Criterion: ", diff, "\n", sep = "")
    
    # Update values
    f_prev = f_step
    l_prev = l_step
  }
} else {
  # In this instance, we only need to update the loadings (k_tau + 1) to N as 
  # the rotation specifies the others remain fixed
  while (diff > tol) {
    iter = iter + 1
    
  }
}
















