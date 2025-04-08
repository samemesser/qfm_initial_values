########################################
## Initial Values Exploration
# In this file, I do some experimentation with one dataset and see
# what effect the initial value can have
########################################
# Preliminaries
rm(list = ls())

# Load fit_qfm() function
source("./code/fit_qfm.R")

########################################
# Prepare data
########################################
set.seed(40563)
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

## Standard parameters
tol <- 1e-6
maxiter <- 100

## Model fit
fit_model <- fit_qfm(data, qtl, k_tau, tol, maxiter)

# Retrieve estimated factors
est_f <- fit_model[[1]]

# Rotation-invariant measure of fit - output depends on dgp
beta_fit <- summary(lm(t(beta) ~ est_f))$r.squared


fit_qfm <- function(data, qtl, k_tau, tol, maxiter, verbose = FALSE) {
  # Helper function to get objective function values
  get_check <- function(data, estimate, quantile) {
    residual <- data - estimate
    check <- residual * (quantile - 1 * (residual <= 0))
    sum(check)
  }
  ## Get data dimensions
  t_ <- dim(data)[1]
  n_ <- dim(data)[2]
  ## Get starting values for estimation procedure
  # Begin with PCA
  f_tilde <- as.matrix(eigen(data %*% t(data))$vectors[, 1:k_tau] * sqrt(t_))
  lam_tilde <- t(data) %*% f_tilde / t_

  # Impose restriction that first two loadings are identity
  # From Bai and Ng (2013), use lambda, inverse first K(tau) elements
  lam_1 <- lam_tilde[1:k_tau, 1:k_tau]

  # F_hat is then F_tilde * t(lam_1)
  f_hat <- f_tilde %*% t(as.matrix(lam_1))
  # lam_hat is lam_tilde * inverse of lam_1
  lam_hat <- lam_tilde %*% solve(lam_1)

  # Impose loadings so that we don't end up with small computational errors
  for (j in 1:k_tau) {
    lam_hat[j, ] <- c(rep(0, j - 1), 1, rep(0, k_tau - j))
  }
  # Starting values are f_hat, lam_hat
  ##################################
  ## Iterative portion of the process
  # Initialize for looping - need 2 copies of each
  f_prev <- f_hat
  l_prev <- lam_hat
  f_step <- f_prev
  l_step <- l_prev

  # While loop initialization - infitie difference, 0 iterations completed to start
  diff <- Inf
  iter <- 0

  while (diff > tol && iter < maxiter) {
    # Increment
    iter <- iter + 1
    # Step factors
    # QR y ~ lambda for each t
    f_step <- t(apply(data, 1, function(x) rq(x ~ l_prev + 0, tau = qtl)$coefficients))
    # Dimensions are wonky when k_tau = 1, this should fix
    if (k_tau == 1) {
      f_step <- t(f_step)
    }
    # Step loadings
    # QR y ~ loading for each n from k_tau + 1 onwards
    l_mid <- t(apply(data[, (k_tau + 1):n_], 2, function(x) rq(x ~ f_step + 0, tau = qtl)$coefficients))
    # Fix dimensions
    if (k_tau == 1) {
      l_mid <- t(l_mid)
    }
    l_step <- rbind(l_prev[1:k_tau, ], l_mid)
    # Check distance
    diff <- sqrt(sum((f_step %*% t(l_step) - f_prev %*% t(l_prev))^2) / (n_ * t_))

    # Update starting
    f_prev <- f_step
    l_prev <- l_step

    if (verbose) {
      # Check obj function value - for use when running inline
      estimate <- f_step %*% t(l_step)
      obj_val <- get_check(data, estimate, qtl)
      cat("Iteration ", iter, ", Convergence Criterion: ", diff, ", Objective Function Value: ", obj_val, "\n", sep = "") #nolint
    }
  }
  # Returns list, first value is matrix of estimated factors,
  # second is matrix of estimated loadings,
  # third is number of iterations required
  list(f_hat = f_step, l_hat = l_step, iter = iter)
}