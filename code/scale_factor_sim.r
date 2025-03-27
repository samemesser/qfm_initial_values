# Only necessary library is quantreg
# Clear everything
rm(list = ls())

# Load library
library(quantreg)
library(foreach)
library(doParallel)
library(doRNG)

# Fix seed
set.seed(8675309)

dgp_2 <- function(n_, t_, scale_factor_variance, error_variance) {
  ## DGP 2
  # alpha_i drawn independently from N(0, 1)
  alphas <- matrix(rnorm(n_), nrow = 1)
  # beta_t drawn independently from N(0, 1)
  betas <- matrix(rnorm(t_), nrow = 1)
  # x_t drawn independently from N(0, 1)
  xs <- matrix(scale_factor_variance * rnorm(t_), nrow = 1)
  # v_it drawn independently from N(0, 1)
  vs <- matrix(rnorm(t_ * n_), nrow = t_, ncol = n_)

  # Gamma is exponential of x, want a TxN matrix for elementwise multiplication
  gammas <- exp(xs)
  gamma_build <- matrix(gammas, nrow = t_, ncol = n_, byrow = TRUE)

  # Combine to create observed variable y_it
  y_it <- t(betas) %*% alphas + gamma_build * vs

  list(y_it, betas, gammas)
}

fit_qfm <- function(data, qtl, k_tau) {
  # Data dimensions, T rows, N columns
  t_ <- dim(data)[1]
  n_ <- dim(data)[2]
  # Initialization of intercepts
  b_est <- apply(data, 2, function(y) rq(y ~ 1, tau = qtl)$coefficients)
  ## PCA for initialization of factors and loadings
  z_res <- sweep(data, 2, b_est)
  f_est <- as.matrix(eigen(z_res %*% t(z_res))$vectors[, 1:k_tau] * sqrt(t_))
  l_est <- t(z_res) %*% f_est / t_
  # eiv rotation so we don't have to rerotate every loop
  # doesn't seem to save computation time and may be giving wrong answer
  #   l1 <- l_est[1:k_tau, 1:k_tau] #nolint
  #   f_est <- f_est %*% t(as.matrix(l1)) #nolint
  #   l_est <- l_est %*% solve(l1) #nolint
  # Save these for looping updates
  b_est_old <- b_est
  fl_est_old <- f_est %*% t(l_est)
  ## Parameters for iteration
  conv_crit <- Inf
  tol <- 1e-6
  iter <- 0
  maxiter <- 1000

  ## Iteration
  while (conv_crit > tol && iter < maxiter) {
    iter <- iter + 1
    # Update factors
    f_est <- t(apply(data, 1, function(y) rq(y - b_est ~ l_est + 0, tau = qtl)$coefficients)) #nolint - VSCode
    # Rotate
    qr_f <- qr(f_est)
    f_est <- sqrt(t_) * qr.Q(qr_f) %*% qr.Q(qr(qr.R(qr_f) %*% t(l_est)))
    # f_est <- sqrt(t_) * qr.Q(qr(f_est)) %*% qr.Q(qr(qr.R(qr(f_est)) %*% t(l_est))) #nolint - VSCode
    # Step 6
    # Update loadings
    coefs <- t(apply(data, 2, function(y) rq(y ~ f_est, tau = qtl)$coefficients)) #nolint - VSCode
    l_est <- coefs[, 2:(k_tau + 1)]
    b_est <- coefs[, 1]
    # Combine to common component
    fl_est <- f_est %*% t(l_est)
    # Check for convergence
    conv_crit <- sum((b_est_old - b_est)^2) / (n_) +
      sum((fl_est_old - fl_est)^2) / (t_)
    # Update stored values
    b_est_old <- b_est
    fl_est_old <- fl_est
  }

  list(fl_est, f_est, l_est, b_est, iter, conv_crit)
}

data <- dgp_2(n_, t_, 0.5, 1)[[1]]
output <- fit_qfm(data, qtl, k_tau)

###########################
## Simulation Code
###########################
## Simulation controls
# Which quantile
qtl <- c(0.25, 0.50, 0.75)
k_tau <- 2
n_sims <- 100


# Set of T and N over which to run simulation
n_vals <- c(50)
t_vals <- c(200)

# Make df that has the particular qtl x n x t combo for each iteration
sim_sets <- expand.grid(n_vals, t_vals, qtl)
sim_sets <- sim_sets[rep(seq_len(nrow(sim_sets)), each = n_sims), ]

# Start cluster
ncore <- detectCores()
cl <- makeCluster(ncore - 1, type = "PSOCK") #ncore - 1 usually
registerDoParallel(cl)

sim_results <- foreach(i = 1:(dim(sim_sets)[1]),
                       .combine = rbind,
                       .export = c("sim_sets", "k_tau", "dgp_2", "fit_qfm"),
                       .packages = c("quantreg")
) %dorng% {
  # prepare output row
  out <- rep(NA, 6)

  # Get quantile and t and n for this iteration
  n_ <- sim_sets[i, 1]
  t_ <- sim_sets[i, 2]
  qtl <- sim_sets[i, 3]

  data <- dgp_2(n_, t_, 0.5, 1)
  fit_model <- fit_qfm(data[[1]], qtl, k_tau)
  
  ff_fit <- abs(cor(fit_model[[2]][, 1], t(data[[2]])))
  sf_r2 <- summary(lm(t(data[[3]]) ~ fit_model[[2]]))$r.squared
  sf_fit <- abs(cor(fit_model[[2]][, 2], t(data[[3]])))

  out[1] <- n_
  out[2] <- t_
  out[3] <- qtl
  out[4] <- ff_fit
  out[5] <- sf_r2
  out[6] <- sf_fit
  out
}

# End cluster
stopCluster(cl)