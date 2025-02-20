# This is a cleaned up_reg version of the rep_reglication code
# from Ando and Bai (2021)
library(quantreg)

# Select seed
set.seed(8675309)

############################################
# Data Construction
############################################

# n_fac - number of common factors
n_fac <- 5
# p_reg - number of common regressors
p_reg <- 8
# e_var - Variance of the error term
e_var <- 1
# t_count - number of time periods
t_count <- 100
# n_asset - number of assets
n_asset <- 100

# tau - chosen quantile
tau <- 0.05
# Random quantiles TxN
u_ <- matrix(runif(t_count * n_asset, 0, 1), nrow = t_count, ncol = n_asset)
# l_true - true factor loadings
# f_true - true factors
l_true <- matrix(runif(n_asset * n_fac, -1, 1), nrow = n_asset, ncol = n_fac)
f_true <- matrix(runif(t_count * n_fac, 0, 2), nrow = t_count, ncol = n_fac)

# Empty matricies to fill with true loadings and factors
true_fl <- matrix(0, nrow = t_count, ncol = n_asset)
obs_fl <- matrix(0, nrow = t_count, ncol = n_asset)

# observed and true common components
for (i in 1:t_count){
  for (j in 1:n_asset){
    obs_b <- l_true[j, ] + 0.1 * u_[i, j]
    true_b <- l_true[j, ] + 0.1 * tau
    if (u_[i, j] <= 0.2) {
      obs_fl[i, j] <- f_true[i, 1:3] %*% obs_b[1:3]
      true_fl[i, j] <- f_true[i, 1:3] %*% true_b[1:3]
    } else if (0.2 <= u_[i, j] && u_[i, j] <= 0.8) {
      obs_fl[i, j] <- f_true[i, 1:4] %*% obs_b[1:4]
      true_fl[i, j] <- f_true[i, 1:4] %*% true_b[1:4]
    } else {
      obs_fl[i, j] <- f_true[i, 1:5] %*% obs_b[1:5]
      true_fl[i, j] <- f_true[i, 1:5] %*% true_b[1:5]
    }
  }
}

obs_x <- matrix(runif(p_reg * n_asset * t_count, 0, 1),
                nrow = n_asset * t_count)

for (i in 1:t_count) {
  for (j in 1:n_asset) {
    a <- f_true[i, 1]^2
    b <- l_true[j, 1] + 0.01 * u_[i, j]
    obs_x[(t_count * (j - 1) + 1):(t_count * j), 1] <-
      obs_x[(t_count * (j - 1) + 1):(t_count * j), 1] + 0.02 * a + 0.02 * b
    a <- f_true[i, 2]^2
    b <- l_true[j, 2] + 0.01 * u_[i, j]
    obs_x[(t_count * (j - 1) + 1):(t_count * j), 3] <-
      obs_x[(t_count * (j - 1) + 1):(t_count * j), 3] - 0.01 * a + 0.02 * b
    a <- f_true[i, 3]^2
    b <- l_true[j, 3] + 0.01 * u_[i, j]
    obs_x[(t_count * (j - 1) + 1):(t_count * j), 5] <-
      obs_x[(t_count * (j - 1) + 1):(t_count * j), 5] + 0.02 * a + 0.02 * b
  }
}

xb <- matrix(0, nrow = t_count, ncol = n_asset)
txb <- xb

for (i in 1:t_count) {
  for (j in 1:n_asset) {
    x <- obs_x[(t_count * (j - 1) + 1):(t_count * j), ]
    obs_b <- c(-1, 1, -1, 1, -1, 1, -1, -1) + j / n_asset + 0.1 * u_[i, j]
    true_b <- c(-1, 1, -1, 1, -1, 1, -1, -1) + j / n_asset + 0.1 * tau
    xb[i, j] <- x[i, 1:p_reg] %*% obs_b[1:p_reg]
    txb[i, j] <- x[i, 1:p_reg] %*% true_b[1:p_reg]
  }
}

tb <- matrix(0, p_reg, n_asset)
for (j in 1:n_asset) {
  tb[, j] <- c(-1, 1, -1, 1, -1, 1, -1, -1) + j / n_asset + 0.1 * tau
}


true_y <- txb + true_fl + qnorm(tau, 0, e_var)

obs_y <- xb + obs_fl + qnorm(u_, 0, e_var)

# sum(true_y - TXBFL)
# sum(obs_y - AY)

###########
# Above is correct
###########

#===================================#
#Estimation
#===================================#
# AB skip the data-driven method for selecting n_fac because the number
# is known from the DGP

if (tau == 0.05) {
  n_fac <- 3
}
if (tau == 0.50) {
  n_fac <- 4
}
if (tau == 0.95) {
  n_fac <- 5
}

# Initialization - empty matrices for estimation
# xb_est - T x N
# fl_est - T x N
# b_est - (P + 1) x N - includes intercept
xb_est <- matrix(0, nrow = t_count, ncol = n_asset)
fl_est <- matrix(0, nrow = t_count, ncol = n_asset)
b_est <- matrix(0, nrow = p_reg + 1, ncol = n_asset)

# Starting guess for betas
for (j in 1:n_asset) {
  y <- obs_y[, j]
  x <- obs_x[(t_count * (j - 1) + 1):(t_count * j), ]
  # Start with quantile regression of y on x with an intercept
  fit <- rq(y ~ x, tau = tau)
  # Store coefficents and xb to subtract from observed y to initialize
  b_est[, j] <- fit$coefficients
  xb_est[, j] <- cbind(1, x) %*% b_est[, j]
}

# Do PCA + QR rotation to get starting value for F and L
z_res <- obs_y - xb_est
vec <- eigen(z_res %*% t(z_res))$vectors
f_est <- sqrt(t_count) * vec[, 1:n_fac]
l_est <- t(solve(t(f_est) %*% f_est) %*% t(f_est) %*% z_res)
fl_est <- f_est %*% t(l_est)

b_est_old <- b_est
fl_est_old <- fl_est

################################################################
# Looping
for (ITE in 1:100) {
  for (i in 1:t_count) {
    y <- obs_y[i, ]
    x <- l_est
    # First step: update factor scores
    # Regress y - estimated intercept on loadings for each t
    fit <- rq(y - xb_est[i, ] ~ x + 0, tau = tau)
    f_est[i, ] <- fit$coefficients[1:n_fac]
  }
  # Rotate factor
  qr_f <- qr(f_est)
  f_est <- sqrt(t_count) * qr.Q(qr_f)
  qr_l <- qr(qr.R(qr_f) %*% t(l_est))
  f_est <- sqrt(t_count) * qr.Q(qr_f) %*% qr.Q(qr_l)

  # Using new factors, update loadings and betas
  # For each asset, quantile reg of y on x and f
  for (j in 1:n_asset) {
    y <- obs_y[, j]
    x <- cbind(obs_x[(t_count * (j - 1) + 1):(t_count * j), ], f_est)
    # With intercept here
    fit <- rq(y ~ x, tau = tau)
    b_est[, j] <- fit$coefficients[1:(p_reg + 1)]
    l_est[j, ] <- fit$coefficients[(p_reg + 2):(p_reg + n_fac + 1)]
    xb_est[, j] <- cbind(1, obs_x[(t_count * (j - 1) + 1):(t_count * j), ]) %*%
      b_est[, j]
    fl_est[, j] <- f_est %*% l_est[j, ]
  }
  
  browser()

  # Check for convergence
  conv_crit <- sum((b_est_old - b_est)^2) / (length(b_est)) +
    sum((fl_est_old - fl_est)^2) / (length(fl_est))
  # Update stored values
  b_est_old <- b_est
  fl_est_old <- fl_est
  print(conv_crit)
  # End looping if under tolerance
  if (conv_crit <= 1e-3) {
    break
  }
}


print(b_est)
