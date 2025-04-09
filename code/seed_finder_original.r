########################################
## Seed Finder
# Finds seeds for dgp that produce desired
# fits based on usual estimation procedure
########################################
# Preliminaries
rm(list = ls())

# Packages
library(quantreg)
# Load fit_qfm() function
source("./code/fit_qfm.R")

########################################
# Prepare data
########################################

seeds_to_test <- floor(runif(1000) * 1000000)

for (i in seq_along(seeds_to_test)) {
  set.seed(seeds_to_test[i])

  cat("Trying seed: ", seeds_to_test[i], "\n", sep = "")

  # Initial exploration at 25th percentile
  qtl <- 0.25
  k_tau <- 2

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

  # Location scale model
  gamma_fit <- summary(lm(t(gamma_c) ~ est_f))$r.squared

  # Print results
  cat("Converged in ", fit_model[[3]], " iterations.\nFirst Factor fit: ", beta_fit, "\nSecond Factor fit: ", gamma_fit, "\n", sep = "" ) #nolint

  if (beta_fit > 0.95 && gamma_fit > 0.75 && fit_model[[3]] < 50) {
    cat("#######################################\n")
    cat("Seed found! Use ", seeds_to_test[i], "\n", sep = "")
    break
  }
}