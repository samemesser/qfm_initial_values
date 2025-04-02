## This code runs a short example of the current QPC procedure.

## quantreg necessary for the function
library(quantreg)

## Working directory needs to be set to the same location as fit_qfm.R
setwd(".")
source("./code/fit_qfm.R")

set.seed(8675309)

qtl <- 0.5
k_tau <- 1
tol <- 1e-6
maxiter <- 1000
n_ <- 50
t_ <- 200
f_true <- matrix(rnorm(t_), nrow = 1)
l_true <- matrix(rnorm(n_), nrow = 1)
e_it <- matrix(rnorm(t_ * n_), nrow = t_, ncol = n_)
data <- t(f_true) %*% l_true + e_it
fit_model <- fit_qfm(data, qtl, k_tau, tol, maxiter)
