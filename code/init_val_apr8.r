########################################
## Initial Values Exploration
# Find seed using seed finder, explore what 
# changing initial estimation values does
########################################
# Preliminaries
rm(list = ls())

# Packages
library(quantreg)
library(tidyverse)
library(xtable)

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
## PCA
pca_f <- as.matrix(eigen(data %*% t(data))$vectors[, 1:k_tau] * sqrt(t_))
pca_l <- t(data) %*% pca_f / t_
fit_model <- fit_qfm_ivdiff(data, qtl, k_tau, pca_f, pca_l, tol, maxiter)

# Retrieve estimated factors
est_f <- fit_model$f_hat

# Rotation-invariant measure of fit - output depends on dgp
beta_fit <- summary(lm(t(beta) ~ est_f))$r.squared

# Location scale model
gamma_fit <- summary(lm(t(gamma_c) ~ est_f))$r.squared

# Print results
cat("Converged in ", fit_model$iter, " iterations.\nFirst Factor fit: ", beta_fit, "\nSecond Factor fit: ", gamma_fit, "\n", sep = "" ) #nolint

checked <- c("No Noise", beta_fit, gamma_fit, fit_model$iter, fit_model$obj_val)

## PCA + Noise

# Find seed for noise that is highly effective
num_seeds <- 10
to_try <- floor(runif(num_seeds) * 1000000)

for (i in seq_along(to_try)) {
  out <- rep(NA, 5)
  out[1] <- to_try[i]
  set.seed(to_try[i])

  cat("Trying seed ", i, "/", num_seeds, ": ", to_try[i], "\n", sep = "")

  pca_start_f <- as.matrix(eigen(data %*% t(data))$vectors[, 1:k_tau] * sqrt(t_))
  pca_start_l <- t(data) %*% pca_start_f / t_
  noise_f <- matrix(rnorm(prod(dim(pca_start_f))), nrow = dim(pca_start_f)[1], ncol = dim(pca_start_f)[2])
  noise_l <- matrix(rnorm(prod(dim(pca_start_l))), nrow = dim(pca_start_l)[1], ncol = dim(pca_start_l)[2])

  pca_noise_f <- pca_start_f + noise_f
  pca_noise_l <- pca_start_l + noise_l

  fit_model <- fit_qfm_ivdiff(data, qtl, k_tau, pca_noise_f, pca_noise_l, tol, maxiter)

  # Retrieve estimated factors
  est_f <- fit_model$f_hat

  # Rotation-invariant measure of fit - output depends on dgp
  beta_fit <- summary(lm(t(beta) ~ est_f))$r.squared

  # Location scale model
  gamma_fit <- summary(lm(t(gamma_c) ~ est_f))$r.squared

  # Print results
  cat("Converged in ", fit_model$iter, " iterations.\nFirst Factor fit: ", beta_fit, "\nSecond Factor fit: ", gamma_fit, "\n", sep = "" ) #nolint
  cat("Objective Function Value: ", fit_model$obj_val, "\n", sep = "")

  out[2] <- beta_fit
  out[3] <- gamma_fit
  out[4] <- fit_model$iter
  out[5] <- fit_model$obj_val

  checked <- rbind(checked, out)
}

rownames(checked) <- NULL

ivdiff_results <- as.data.frame(checked) %>%
  rename(noise_seed = V1, beta_fit = V2, gamma_fit = V3, iter = V4, obj_val = V5) %>%
  mutate(beta_fit = as.numeric(beta_fit),
         gamma_fit = as.numeric(gamma_fit),
         iter = as.integer(iter),
         obj_val = as.numeric(obj_val))

View(ivdiff_results)

out <- xtable(ivdiff_results, align = c("l", rep("c", 5)), digits = 4)

print(out,
      file = "./out/ivdiff_results.tex",
      include.rownames = FALSE,
      include.colnames = FALSE,
      caption.placement = "top",
      add.to.row = list(pos = list(0),
                        command = c(
                          "\nNoise Seed & First Factor Fit & Second Factor Fit & Iterations & Objective Function Value \\\\\n")), #nolint
      floating = FALSE,
      hline.after = 0)