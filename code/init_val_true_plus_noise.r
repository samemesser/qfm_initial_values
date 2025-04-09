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

# Seeds used to make examples in apr 9 version of loaction_scale_notes.pdf

# set.seed(787558) # Good FF fit, bad SF fit
# set.seed(376386) # Bad FF and SF fit
set.seed(325055) # Good FF and SF fit

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

## True values
true_start_f <- cbind(t(beta), t(gamma_c))
true_start_l <- cbind(t(alpha), t(eta_c))
# noise_f <- matrix(rnorm(prod(dim(true_start_f)), sd = 0.05), nrow = dim(true_start_f)[1], ncol = dim(true_start_f)[2])
# noise_l <- matrix(rnorm(prod(dim(true_start_l)), sd = 0.05), nrow = dim(true_start_l)[1], ncol = dim(true_start_l)[2])

true_noise_f <- true_start_f #+ noise_f
true_noise_l <- true_start_l #+ noise_l

fit_model <- fit_qfm_ivdiff(data, qtl, k_tau, true_noise_f, true_noise_l, tol, maxiter)

# Retrieve estimated factors
est_f <- fit_model$f_hat

# Rotation-invariant measure of fit - output depends on dgp
beta_fit <- summary(lm(t(beta) ~ est_f))$r.squared

# Location scale model
gamma_fit <- summary(lm(t(gamma_c) ~ est_f))$r.squared

# Print results
cat("Converged in ", fit_model$iter, " iterations.\nFirst Factor fit: ", beta_fit, "\nSecond Factor fit: ", gamma_fit, "\n", sep = "" ) #nolint
cat("Objective Function Value: ", fit_model$obj_val, "\n", sep = "")

out <- rep(NA, 5)
out[1] <- "True"
out[2] <- beta_fit
out[3] <- gamma_fit
out[4] <- fit_model$iter
out[5] <- fit_model$obj_val

checked <- rbind(checked, out)

rownames(checked) <- NULL

ivdiff_results <- as.data.frame(checked) %>%
  rename(init_val = V1, beta_fit = V2, gamma_fit = V3, iter = V4, obj_val = V5) %>%
  mutate(beta_fit = as.numeric(beta_fit),
         gamma_fit = as.numeric(gamma_fit),
         iter = as.integer(iter),
         obj_val = as.numeric(obj_val))

View(ivdiff_results)

# out <- xtable(ivdiff_results, align = c("l", rep("c", 5)), digits = 4)

# print(out,
#       file = "./out/ivdiff_results_goodbf.tex",
#       include.rownames = FALSE,
#       include.colnames = FALSE,
#       caption.placement = "top",
#       add.to.row = list(pos = list(0),
#                         command = c(
#                           "\nNoise Seed & First Factor Fit & Second Factor Fit & Iterations & Objective Function Value \\\\\n")), #nolint
#       floating = FALSE,
#       hline.after = 0)