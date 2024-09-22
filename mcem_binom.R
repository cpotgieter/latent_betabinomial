# =============================================================================
# File: mcem_binom.R
# Author: Cornelis J Potgieter
#
# -----------------------------------------------------------------------------
# Description:
#   This script implements the `run.mcem.binom` function, which applies the 
#   MCEM algorithm for parameter estimation in a Binomial variation of the 
#   count-time model. The function iteratively imputes latent variables and 
#   updates model parameters to approximate the maximum likelihood estimates.
#
#   - `run.mcem.binom`: Main function executing the MCEM algorithm.
#
# Inputs:
#     - Y       : Matrix of count data (n x J), where n is the number of observations
#                 and J is the number of items.
#     - logT10  : Matrix of log-time data (n x J), corresponding to the counts in Y.
#     - N       : Vector of item lengths (length J).
#     - J       : Integer representing the number of items.
#     - k.in    : Integer or vector controlling the number of imputations in the
#                 E-step for each MCEM run (default is 5).
#     - reps.in : Integer or vector specifying the number of iterations for each
#                 MCEM run (default is 2).
#     - ests.in : List containing initial estimates for model parameters. If not
#                 provided, MOM estimators are used.
#
# Outputs:
#   - MCEM.ests  : A list containing matrices of estimated parameters across
#                 MCEM iterations:
#                   - a       : Estimated parameter 'a' for each passage.
#                   - b       : Estimated parameter 'b' for each passage.
#                   - alpha   : Estimated precision parameter 'alpha' for each passage.
#                   - beta    : Estimated parameter 'beta' for each passage.
#                   - numwords.p : Vector of number of words per passage (input N).
#                   - vartau  : Estimated variance parameter 'vartau' across iterations.
#                   - rho     : Estimated intra-class correlation 'rho' across iterations.
#
# Usage:
#   1. Source this script to make the `run.mcem.binom` function available:
#        source("path/to/mcem_binom.R")
#
#   2. (Optional) Compute initial parameter estimates using MOM functions:
#        ests.in <- mom.binom(Y, logT10, N, J)
#
#   3. Run the MCEM algorithm:
#        MCEM_results <- run.mcem.binom(Y, logT10, N, J, k.in = 5, reps.in = 2, ests.in)
#
#   4. Access the estimated parameters:
#        estimated_a <- MCEM_results$a
#        estimated_b <- MCEM_results$b
#        # And so on for other parameters
#
# Notes:
#   - Ensure that the input matrices `Y` and `logT10` have the same dimensions.
#   - Set missing data to NA before passing matrices to the functions.
#   - The function includes internal helper functions for negative log-likelihood
#     computation, rejection sampling, imputation, and optimization steps.
#
# =============================================================================

# -------------------------------
# Function Definitions
# -------------------------------

# Main function implementing the MCEM algorithm for the binomial case
# Parameters: 
#   Y: Matrix of count data
#   logT10: Matrix of log-scale speed data
#   N: Vector with number of words in each passage
#   J: Number of passages
#   k.in: Controls the number of imputations per iteration (default is 5)
#   reps.in: Controls the number of MCEM algorithm iterations (default is 2)
#   ests.in: Initial values for model parameters (estimates)

run.mcem.binom <- function(Y, logT10, N, J, k.in=5, reps.in=2, ests.in) {
  
  # ---------------------------------------------------------------------------
  # Internal Function: neg_logratio_binom
  # Purpose:
  #   Computes the negative log-likelihood for the binomial model. This function
  #   is utilized in the rejection sampling step of the MCEM sampler to evaluate
  #   the suitability of sampled latent variables.
  #
  # Inputs:
  #   - data : A matrix where each row corresponds to a different parameter:
  #           Row 1: Y (counts)
  #           Row 2: N (number of words)
  #           Row 3: a (parameter 'a')
  #           Row 4: b (parameter 'b')
  #   - par  : Scalar latent variable 'z' to be evaluated.
  #
  # Output:
  #   - negLP : Negative log-likelihood value for the given parameters and latent variable.
  # ---------------------------------------------------------------------------
  
  neg_logratio_binom <- function(data, par) {
    Y <- data[1, ]
    N <- data[2, ]
    a <- data[3, ]
    b <- data[4, ]
    z <- par
    
    # Calculate success probabilities using the probit link function
    success.prob <- pnorm(a * z + b)
    
    # Calculate negative log-likelihood for binomial cases
    negLP <- -sum(dbinom(Y, N, success.prob, log = TRUE))
    
    return(negLP)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: gamma_multiplier_binom
  # Purpose:
  #   Calculates the gamma multiplier and optimal latent variable 'zmax' for the
  #   rejection sampler. This function determines the acceptance rate in the
  #   rejection sampling process.
  #
  # Inputs:
  #   - Y   : Vector of count data for a single individual.
  #   - N   : Vector of the number of words corresponding to Y.
  #   - a   : Vector of parameter 'a' corresponding to Y.
  #   - b   : Vector of parameter 'b' corresponding to Y.
  #
  # Output:
  #   - gamval : A list containing 'zmax' and 'gammax' values for rejection sampling.
  # ---------------------------------------------------------------------------
  
  gamma_multiplier_binom <- function(Y, N, a, b) {
    if (sum(N - Y) < 0.01) {
      gammax <- 1
      zmax <- Inf
    } else {
      zgam <- optimize(neg_logratio_binom, 
                       data = rbind(Y, N, a, b), 
                       interval = c(-7, 7))
      zmax <- zgam$minimum
      gammax <- exp(-zgam$objective)
    }
    
    gamval <- list(zmax = zmax, gammax = gammax)
    return(gamval)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: imputation_code_binom
  # Purpose:
  #   Performs the imputation step of the MCEM algorithm by sampling latent
  #   variables using rejection sampling based on the current parameter estimates.
  #
  # Inputs:
  #   - Y         : Matrix of count data.
  #   - logT10    : Matrix of log-scale speed data.
  #   - N         : Vector of number of words in each passage.
  #   - a, b      : Vectors of current parameter estimates.
  #   - alpha, beta: Vectors of current parameter estimates.
  #   - sigma_tau : Current estimate of sigma_tau (standard deviation).
  #   - rho       : Current estimate of intra-class correlation.
  #   - K         : Number of imputations to perform.
  #
  # Output:
  #   - imputes : A list containing imputed values for theta and tau, and the
  #              optimal z values.
  # ---------------------------------------------------------------------------
  
  imputation_code_binom <- function(Y, logT10, N, a, b, alpha, beta, sigma_tau, rho, K) {
    # Initialization and parameter setup
    n <- dim(Y)[1]
    zmax <- matrix(rep(0, n), nrow = n)
    gammax <- matrix(rep(0, n), nrow = n)
    mu1 <- matrix(rep(0, n), nrow = n)
    A <- matrix(rep(0, n), nrow = n)
    sigma1_2 <- matrix(rep(0, n), nrow = n)
    sigma2_2 <- matrix(rep(0, n), nrow = n)
    
    for (i in 1:n) {
      index <- which(!is.na(Y[i, ]))
      gamval <- gamma_multiplier_binom(Y[i, index], N[index], a[index], b[index])
      gammax[i] <- gamval$gammax
      zmax[i] <- gamval$zmax
      A[i] <- sum(alpha[index]^2)
      mu1[i] <- -rho * sigma_tau / (1 + A[i] * sigma_tau^2) * sum(alpha[index]^2 * (logT10[i, index] - beta[index]))
      sigma1_2[i] <- 1 / (1 + (A[i] * rho^2 * sigma_tau^2) / (1 + A[i] * (sigma_tau^2 - rho^2 * sigma_tau^2)))
      sigma2_2[i] <- (sigma_tau^2 - rho^2 * sigma_tau^2) / ((sigma_tau^2 - rho^2 * sigma_tau^2) * A[i] + 1)
    }
    
    Z1 <- matrix(rep(0, n * K), nrow = n, ncol
                 