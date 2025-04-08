## This code runs a brief simulation study with 3 different DGPs

## quantreg necessary for the function
library(tidyverse)
library(quantreg)
library(doParallel)
library(foreach)
library(doRNG)

## Preliminaries
rm(list = ls())

## Working directory needs to be set to the same location as fit_qfm.R
setwd(".")
source("./code/fit_qfm.R")

# Replicability
set.seed(67189426)

## Data generation
n_ <- 50
t_ <- 200

## Option for DGP - 1f, 2f, loc-scale
# dgp_spec <- "1f"
# dgp_spec <- "2f"
dgp_spec <- "loc_scale"

# Always work with same random vectors
beta <- matrix(rnorm(t_), nrow = 1)
alpha <- matrix(rnorm(n_), nrow = 1)

# Second factor if can be negative
gamma_n <- matrix(rnorm(t_), nrow = 1)
eta_n <- matrix(rnorm(n_), nrow = 1)

# Second factor if must be positive
gamma_c <- matrix(rchisq(t_, df = 1), nrow = 1)
eta_c <- matrix(rchisq(n_, df = 1), nrow = 1)

e_it <- matrix(rnorm(t_ * n_), nrow = t_, ncol = n_)

if (dgp_spec == "1f") {
  # One factor model
  data <- t(beta) %*% alpha + e_it
} else if (dgp_spec == "2f") {
  # Two factor model
  data <- t(beta) %*% alpha + t(gamma_n) %*% eta_n + e_it
} else if (dgp_spec == "loc_scale") {
  # Location scale model
  data <- t(beta) %*% alpha + t(gamma_c) %*% eta_c * e_it
}

## Standard parameters
tol <- 1e-6
maxiter <- 1000

## Estimation location
qtl <- 0.25
k_tau <- 2


## Model fit
fit_model <- fit_qfm(data, qtl, k_tau, tol, maxiter)

# Retrieve estimated factors
est_f <- fit_model[[1]]

# Rotation-invariant measure of fit - output depends on dgp
beta_fit <- summary(lm(t(beta) ~ est_f))$r.squared
# DGP variable outcomes
if (dgp_spec == "1f") {
  # One factor model
  gamma_fit <- NA
} else if (dgp_spec == "2f") {
  # Two factor model
  gamma_fit <- summary(lm(t(gamma_n) ~ est_f))$r.squared
} else if (dgp_spec == "loc_scale") {
  # Location scale model
  gamma_fit <- summary(lm(t(gamma_n) ~ est_f))$r.squared
}
# Print results
cat("Beta fit: ", beta_fit, "\nGamma fit: ", gamma_fit, "\nConverged in: ", fit_model$iter, " iterations.\n", sep = "") #nolint 