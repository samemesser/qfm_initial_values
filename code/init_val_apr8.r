########################################
## Initial Values Exploration
# Find seed using seed finder, explore what 
# changing initial estimation values does
########################################
# Preliminaries
rm(list = ls())

# Packages
library(quantreg)

set.seed(787558)

### Data generation
## This seed gives good fit of first factor and bad fit of second at 0.25
# Data size
t_ <- 200
n_ <- 50

# Generate data
beta <- matrix(rnorm(t_), nrow = 1)
alpha <- matrix(rnorm(n_), nrow = 1)

# Second factor if must be positive
gamma_c <- matrix(rchisq(t_, df = 1), nrow = 1)
eta_c <- matrix(rchisq(n_, df = 1), nrow = 1)

e_it <- matrix(rnorm(t_ * n_), nrow = t_, ncol = n_)

# Location scale data
data <- t(beta) %*% alpha + t(gamma_c) %*% eta_c * e_it

### Testing
## Set function parameters, load function
source("./code/fit_qfm_ivdiff.R")
# Set number of factors and quantile to explore
qtl <- 0.25
k_tau <- 2
# For convergence
tol <- 1e-6
maxiter <- 100

## Get starting values for estimation procedure
# Begin with PCA
f_tilde <- as.matrix(eigen(data %*% t(data))$vectors[, 1:k_tau] * sqrt(t_))
lam_tilde <- t(data) %*% f_tilde / t_

fit_model <- fit_qfm_ivdiff(data, qtl, k_tau, f_tilde, lam_tilde, tol, maxiter)

# Retrieve estimated factors
est_f <- fit_model[[1]]

# Rotation-invariant measure of fit - output depends on dgp
beta_fit <- summary(lm(t(beta) ~ est_f))$r.squared

# Location scale model
gamma_fit <- summary(lm(t(gamma_c) ~ est_f))$r.squared

# Print results
cat("Converged in ", fit_model[[3]], " iterations.\nFirst Factor fit: ", beta_fit, "\nSecond Factor fit: ", gamma_fit, "\n", sep = "" ) #nolint
