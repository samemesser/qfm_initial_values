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

# Set of quantiles
qtl <- c(0.25, 0.50, 0.75)
# Think we can always estimate 2 factors in a 1f model
k_tau <- 2
n_sims <- 5

# Set of T and N over which to run simulation
n_vals <- c(50)
t_vals <- c(200)

# DGP specifications
dgp_specs <- c("1f", "2f", "loc_scale")
# Make df that has the particular qtl x n x t x DGP combo for each iteration
sim_sets <- expand.grid(n_vals, t_vals, qtl, dgp_specs)
sim_sets <- sim_sets[rep(seq_len(nrow(sim_sets)), each = n_sims), ]

# Time
start <- proc.time()[3]
# Start cluster
ncore <- detectCores()
cl <- makeCluster(ncore - 1, type = "PSOCK") #ncore - 1 usually
registerDoParallel(cl)

sim_results <- foreach(i = 1:(dim(sim_sets)[1]),
                       .combine = rbind,
                       .export = c("sim_sets", "k_tau", "fit_qfm"),
                       .packages = c("quantreg")
) %dorng% {
  # Prepare output row
  out <- rep(NA, 7)
  # Get quantile and t and n for this iteration
  n_ <- sim_sets[i, 1]
  t_ <- sim_sets[i, 2]
  qtl <- sim_sets[i, 3]
  dgp_spec <- sim_sets[i, 4]
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
  maxiter <- 100

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
    gamma_fit <- summary(lm(t(gamma_c) ~ est_f))$r.squared
  }
  out[1] <- n_
  out[2] <- t_
  out[3] <- qtl
  out[4] <- beta_fit
  out[5] <- gamma_fit
  out[6] <- fit_model$iter
  out[7] <- dgp_spec

  out
}

sim_results <- as.data.frame(sim_results)
names(sim_results) <- c("n_", "t_", "qtl", "beta_fit", "gamma_fit", "iter", "dgp")
sim_results$dgp <- recode(sim_results$dgp,
                          `1` = "1f",
                          `2` = "2f",
                          `3` = "loc_scale")

# End cluster
stopCluster(cl)

agg_result <- sim_results %>%
  group_by(n_, t_, qtl, dgp) %>%
  summarise(mean_ff_fit = mean(beta_fit),
            mean_sf_fit = mean(gamma_fit),
            long_converge = mean((iter > 50) & (iter < 100)),
            no_converge = mean(iter == 100),
            prop_gt_50 = mean(gamma_fit > 0.5),
            prop_ff_wrong = mean(beta_fit < 0.9))

cat("Simulations complete in ", proc.time()[3] - start, " seconds.\n", sep = "")

write.csv(sim_results, "./out/individual_results.csv")
write.csv(agg_result, "./out/agg_results.csv")

View(agg_result)    