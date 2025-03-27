###### simulate_qfm.R
# This file generates data for and simulates a location-scale QFM from Sagner
######
## Preliminaries
# Libraries
library(tidyverse)
library(lubridate)
library(quantreg)
library(doParallel)
library(foreach)
library(doRNG)

# Clear environment
rm(list = ls())

# Set seed for replicability
set.seed(67189426)

dgp_mar27 <- function(n_, t_) {
  ## DGP example 3 from Chen et al.
  # alpha_i drawn independently from N(0, 1)
  alphas <- matrix(rnorm(n_), nrow = 1)
  # beta_t drawn independently from N(0, 1)
  betas <- matrix(rnorm(t_), nrow = 1)
  # x_t drawn independently from N(0, 1)
  gammas <- matrix(rchisq(t_, df = 1), nrow = 1)
  # y_t drawn independently from N(0, 1) as loading for volatility factor
  etas <- matrix(rchisq(n_, df = 1), nrow = 1)
  # v_it drawn independently from N(0, 1)
  vs <- matrix(rnorm(t_ * n_), nrow = t_, ncol = n_)

  # Gamma is exponential of x, want a TxN matrix for elementwise multiplication
  #gammas <- exp(xs)
  # gamma_build <- matrix(gammas, nrow = t_, ncol = n_, byrow = TRUE)

  # Eta is an exponential of y
  #etas <- exp(ys)
  # Combine to create observed variable y_it
  y_it <- t(betas) %*% alphas + t(gammas) %*% etas * vs

  list(y_it, alphas, betas, etas, gammas)
}

dgp_ex_3 <- function(n_, t_, scale_factor_variance, error_variance) {
  ## DGP example 3 from Chen et al.
  # alpha_i drawn independently from N(0, 1)
  alphas <- matrix(rnorm(n_), nrow = 1)
  # beta_t drawn independently from N(0, 1)
  betas <- matrix(rnorm(t_), nrow = 1)
  # x_t drawn independently from N(0, 1)
  xs <- matrix(scale_factor_variance * rnorm(t_), nrow = 1)
  # y_t drawn independently from N(0, 1) as loading for volatility factor
  ys <- matrix(scale_factor_variance * rnorm(n_), nrow = 1)
  # v_it drawn independently from N(0, 1)
  vs <- matrix(rnorm(t_ * n_), nrow = t_, ncol = n_)

  # Gamma is exponential of x, want a TxN matrix for elementwise multiplication
  gammas <- exp(xs)
  # gamma_build <- matrix(gammas, nrow = t_, ncol = n_, byrow = TRUE)

  # Eta is an exponential of y
  etas <- exp(ys)
  # Combine to create observed variable y_it
  y_it <- t(betas) %*% alphas + t(gammas) %*% etas * vs

  list(y_it, alphas, betas, etas, gammas)
}

get_check <- function(data, estimate, quantile) {
  residual <- data - estimate
  check <- residual * (quantile - 1 * (residual <= 0))
  sum(check)
}


########################
## Set up
n_ <- 50
t_ <- 200

n_sim <- 10

out <- rep(NA, 2)

for (i in 1:n_sim) {
  #data <- dgp_ex_3(n_, t_, 0.5, 1)
  data <- dgp_mar27(n_, t_)
  y_it <- data[[1]]
  true_alphas <- t(data[[2]])
  true_betas <- t(data[[3]])
  true_etas <- t(data[[4]])
  true_gammas <- t(data[[5]])

  ########################
  ## Sagner (2019) Estimation
  qtl <- 0.25
  k_tau <- 2

  #######
  ## Get initial values
  # Matrix to take eigenvalues from
  pca_mat <- y_it %*% t(y_it)
  f_tilde <- as.matrix(eigen(pca_mat)$vectors[, 1:k_tau] * sqrt(t_))

  # Get lambda from PCA
  lam_tilde <- t(y_it) %*% f_tilde / t_

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

  # Initialize for looping
  f_prev <- f_hat
  l_prev <- lam_hat

  f_step <- f_prev
  l_step <- l_prev
  
  # While loop initialization
  diff <- Inf
  tol <- 1e-6
  maxiter <- 100
  iter <- 0

  while (diff > tol & iter < maxiter) {
    iter <- iter + 1
    # Step factors
    # QR y ~ lambda for each t
    f_step <- t(apply(y_it, 1, function(x) rq(x ~ l_prev + 0, tau = qtl)$coefficients))

    # Step loadings
    # QR y ~ loading for each n from k_tau + 1 onwards
    l_step <- rbind(l_prev[1:k_tau, ],
                    t(apply(y_it[, (k_tau + 1):n_], 2, function(x) rq(x ~ f_step + 0, tau = qtl)$coefficients)))
    estimate <- f_step %*% t(l_step)
    # Check obj function
    obj_val <- get_check(y_it, estimate, qtl)
    # Check distance
    diff <- sqrt(sum((f_step %*% t(l_step) - f_prev %*% t(l_prev))^2) / (n_ * t_))

    # Update starting
    f_prev <- f_step
    l_prev <- l_step

    #cat("Iteration ", iter, ", Convergence Criterion: ", diff, ", Objective Function Value: ", obj_val, "\n", sep = "")
  }

  qpc_cor_beta <- abs(cor(f_step[, 1], true_betas))
  qpc_cor_gamma <- abs(cor(f_step[, 2], true_gammas))

  qpc_r2_beta <- summary(lm(true_betas ~ f_step))$r.squared
  qpc_r2_gamma <- summary(lm(true_gammas ~ f_step))$r.squared

  tmp <- c(qpc_r2_beta, qpc_r2_gamma)

  out <- rbind(out, tmp)

  cat("Iteration ", i, " complete.\n", sep = "")
}

#out <- c(qpc_cor_beta, qpc_cor_gamma, qpc_r2_beta, qpc_r2_gamma)

#cat("QPC Correlation (Factor 1): ", round(out[1], 4), "\n",
#    "QPC Correlation (Factor 2): ", round(out[2], 4), "\n",
#    "QPC R Squared (Factor 1): ", round(out[3], 4), "\n",
#    "QPC R Squared (Factor 2): ", round(out[4], 4), "\n", sep = "")

# mean(l_step[, 2])

View(out)