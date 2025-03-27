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
n_ <- 50
t_ <- 200
scale_factor_variance <- 1
error_variance <- 1
################################################################################
## DGP
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
################################################################################
# From here, we have to get a starting value for the algorithm. We'll use PCA
# for this. There are several potential rotations, I used errors-in-variables.
# It is also worth considering the rotation used in Ando and Bai (2020).

## Here choose quantile of interest and number of factors
qtl <- 0.25
k_tau <- 2
################################################################################
# General strategy is as follows:
# 1) Estimate qr of y on 1 to get intercept starting value
# 2) Compute 1 x coefficient and subtract from y, call it z
# 3) With this, do PCA + QR rotation to get starting values for F and lambda
# 4) Take intercept + lambda as given, update F using qr without intercept
# 5) Do QR rotation on F
# 6) Update L taking F as given, using qr with intercept
# 7) Convergence then based on how common component + intercepts differ across
# iterations
################################################################################
# Specify data
y_obs <- y_it
# Initialize empty matrices for storage
xb_est <- matrix(0, nrow = t_, ncol = n_)
fl_est <- matrix(0, nrow = t_, ncol = n_)
b_est <- matrix(0, nrow = 1, ncol = n_)
# Define x as a ones vector
x_obs <- matrix(1, nrow = t_)

# Step 1
for (i in 1:n_) {
  y <- y_obs[, i]
  fit <- rq(y ~ x_obs + 0, tau = qtl)
  b_est[, i] <- fit$coefficients
  xb_est[, i] <- x_obs %*% b_est[, i]
}

# Step 2
z_res <- y_obs - xb_est

# Step 3 - PCA
pca_mat <- z_res %*% t(z_res)
f_est <- as.matrix(eigen(pca_mat)$vectors[, 1:k_tau] * sqrt(t_))
l_est <- t(z_res) %*% f_est / t_

b_est_old <- b_est
fl_est_old <- fl_est

################################################################################
## Estimation Parameters
conv_crit <- Inf
tol <- 1e-3
iter <- 0
maxiter <- 100

################################################################################
## Utility function for calculating the value of the objective function
get_check <- function(data, estimate, quantile) {
  residual <- data - estimate
  check <- residual * (quantile - 1 * (residual <= 0))
  sum(check)
}
## Utility function for rotating estimates according to AB
rotate_ab <- function(f_est, l_est, t_) {
  # Rotation step
  qr_f <- qr(f_est)
  rf_tau <- qr.R(qr_f)
  lam_mid <- rf_tau %*% t(l_est)
  qr_l <- qr(lam_mid)
  # Get rotated estimates
  f_rot <- sqrt(t_) * qr.Q(qr_f) %*% qr.Q(qr_l)
  l_rot <- t(qr.R(qr_l))

  list(f_rot, l_rot)
}
################################################################################
## Iteration
while (conv_crit > tol) {
  iter <- iter + 1
  # Step 4
  # Update factors
  for (i in 1:t_) {
    y <- y_obs[i, ]
    # Taking loadings as given
    x <- l_est
    fit <- rq(y - xb_est[i, ] ~ x + 0, tau = qtl)
    f_est[i, ] <- fit$coefficients
  }
  # Rotate factor
  qr_f <- qr(f_est)
  f_est <- sqrt(t_) * qr.Q(qr_f)
  qr_l <- qr(qr.R(qr_f) %*% t(l_est))
  f_est <- sqrt(t_) * qr.Q(qr_f) %*% qr.Q(qr_l)

  # Using new factors, update loadings and betas
  # For each asset, quantile reg of y on x and f
  for (j in 1:n_) {
    y <- y_obs[, j]
    x <- cbind(x_obs, f_est)
    # With intercept here
    fit <- rq(y ~ x + 0, tau = qtl)
    b_est[, j] <- fit$coefficients[1]
    l_est[j, ] <- fit$coefficients[2:(k_tau + 1)]
    xb_est[, j] <-  x_obs %*% b_est[, j]
    fl_est[, j] <- f_est %*% l_est[j, ]
  }
  y_fit <- xb_est + fl_est
  # Get value of objective function
  obj_val <- get_check(y_obs, y_fit, qtl)

  # Check for convergence
  conv_crit <- sum((b_est_old - b_est)^2) / (length(b_est)) +
    sum((fl_est_old - fl_est)^2) / (length(fl_est))
  # Update stored values
  b_est_old <- b_est
  fl_est_old <- fl_est

  cat("Iteration ", iter,
      ", Objective Function Value: ", obj_val,
      ", Convergence Criterion: ", conv_crit, "\n", sep = "")
}

################################################################################
## Compare to true factors
abs(cor(f_est[, 1], t(betas)))
abs(cor(f_est[, 2], t(gammas)))

## Rotation-agnostic measure of fit
summary(lm(t(betas) ~ f_est))$r.squared
summary(lm(t(gammas) ~ f_est))$r.squared
