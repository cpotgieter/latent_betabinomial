# =============================================================================
# File: mcem_binom.R
# Author: Cornelis J Potgieter
#
# ---------------------------------------------------------------------------
# Description:
#   This script implements the `run.mcem.binom` function, which executes the
#   Monte Carlo Expectation-Maximization (MCEM) algorithm for parameter estimation
#   in the Binomial variation of the count-time model. The function iteratively
#   imputes latent variables and updates model parameters to approximate the
#   maximum likelihood estimates.
#
#   - `run.mcem.binom`: Main function executing the MCEM algorithm for the
#     Binomial model.
#
# Inputs:
#     - Y         : Matrix of count data (n x J), where n is the number of observations
#                  and J is the number of items/passages.
#     - logT10    : Matrix of log-scale speed data (n x J), corresponding to the counts in Y.
#     - N         : Vector of item lengths (length J), representing the number of words in each passage.
#     - J         : Integer representing the number of passages/items.
#     - k.in      : Integer or vector controlling the number of imputations in the
#                  E-step for each MCEM run (default is 5).
#     - reps.in   : Integer or vector specifying the number of iterations for each
#                  MCEM run (default is 2).
#     - ests.in   : List containing initial estimates for model parameters. If not
#                  provided, Method of Moments (MOM) estimators are used.
#
# Outputs:
#   - MCEM.ests : A list containing matrices of estimated parameters across
#                 MCEM iterations:
#                   - a       : Estimated parameter 'a' for each passage.
#                   - b       : Estimated parameter 'b' for each passage.
#                   - alpha   : Estimated precision parameter 'alpha' for each passage.
#                   - beta    : Estimated parameter 'beta' for each passage.
#                   - vartau  : Estimated variance parameter 'vartau' across iterations.
#                   - rho     : Estimated intra-class correlation 'rho' across iterations.
#
# Dependencies:
#   - `mvtnorm` package
#
# Usage:
#   1. Load the necessary packages:
#        library(mvtnorm)
#
#   2. Source this script to make the `run.mcem.binom` function available:
#        source("path/to/mcem_binom.R")
#
#   3. (Optional) Compute initial parameter estimates using MOM functions:
#        ests.in <- mom.binom(Y, logT10, N, J)
#
#   4. Run the MCEM algorithm:
#        MCEM_results <- run.mcem.binom(Y, logT10, N, J, k.in = 5, reps.in = 2, ests.in)
#
#   5. Access the estimated parameters:
#        estimated_a <- MCEM_results$a
#        estimated_b <- MCEM_results$b
#        # And so on for other parameters
#
# Notes:
#   - Ensure that the input matrices `Y` and `logT10` have the same dimensions.
#   - Set missing data to NA before passing matrices to the functions.
#   - The function includes internal helper functions for negative log-likelihood
#     computation, rejection sampling, imputation, and optimization steps.
#   - For large datasets, consider the computational time due to optimization and
#     imputation steps.
#   - The function uses `dbinom` for the Binomial distribution.
#
# =============================================================================

# -------------------------------
# Function Definitions
# -------------------------------

# Main function implementing the MCEM algorithm for the Binomial model
# Parameters: (a, b, alpha, beta, sigma_tau, rho)
# Y         : Matrix of count data
# logT10    : Matrix of log-scale speed data
# N         : Vector with number of words in each passage
# J         : Number of passages
# k.in      : Number of Monte Carlo imputations
# reps.in   : Number of MCEM iterations
# ests.in   : List with initial values for model parameters
#
# MCEM specification examples:
#   - k.in = 1 and reps.in = 50
#   - k.in = c(1, 10) and reps.in = c(50, 5)

run.mcem.binom <- function(Y, logT10, N, J, k.in = 5, reps.in = 2, ests.in) {
  
  # ---------------------------------------------------------------------------
  # Internal Function: neg_logratio_binom
  # Purpose:
  #   Computes the negative log-likelihood for the Binomial model. This function
  #   is utilized in the rejection sampling step of the MCEM sampler to evaluate
  #   the suitability of sampled latent variables.
  #
  # Inputs:
  #   - data : A matrix where each row corresponds to different parameters:
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
    # Extracting data elements
    Y <- data[1, ]
    N <- data[2, ]
    a <- data[3, ]
    b <- data[4, ]
    z <- par
    
    # Calculating success probabilities using the probit link function
    success.prob <- pnorm(a * z + b)
    
    # Calculating negative log-likelihood for the Binomial model
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
    # Check for a special case where the count data is perfect (no variability)
    if (sum(N - Y) < 0.01) {
      gammax <- 1
      zmax <- Inf
    } else {
      # Optimizing the negative log-ratio over the interval [-7, 7]
      zgam <- optimize(neg_logratio_binom, 
                       data = rbind(Y, N, a, b), 
                       interval = c(-7, 7))
      zmax <- zgam$minimum
      gammax <- exp(-zgam$objective)  # Calculating gamma multiplier
    }
    
    # Returning the calculated gamma multiplier and optimal z
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
  imputation_code_binom <- function(Y, logT10, N, a, b, alpha, 
                                    beta, sigma_tau, rho, K) {
    # Initialization and parameter setup
    n <- dim(Y)[1]
    zmax <- matrix(rep(0, n), nrow = n)
    gammax <- matrix(rep(0, n), nrow = n)
    mu1 <- matrix(rep(0, n), nrow = n)
    A <- matrix(rep(0, n), nrow = n)
    sigma1_2 <- matrix(rep(0, n), nrow = n)
    sigma2_2 <- matrix(rep(0, n), nrow = n)
    
    # Loop over each individual to set up rejection sampling parameters
    for (i in 1:n) {
      index <- which(!is.na(Y[i, ]))  # Identify non-missing entries
      gamval <- gamma_multiplier_binom(Y[i, index], N[index], 
                                       a[index], b[index])
      gammax[i] <- gamval$gammax
      zmax[i] <- gamval$zmax
      A[i] <- sum(alpha[index]^2)
      mu1[i] <- -rho * sigma_tau / (1 + A[i] * sigma_tau^2) * 
        sum(alpha[index]^2 * (logT10[i, index] - beta[index]))
      sigma1_2[i] <- 1 / (1 + (A[i] * rho^2 * sigma_tau^2) / 
                            (1 + A[i] * (sigma_tau^2 - rho^2 * sigma_tau^2)))
      sigma2_2[i] <- (sigma_tau^2 - rho^2 * sigma_tau^2) / 
        ((sigma_tau^2 - rho^2 * sigma_tau^2) * A[i] + 1)
    }
    
    # Imputation step setup
    Z1 <- matrix(rep(0, n * K), nrow = n, ncol = K)
    mu2 <- matrix(rep(0, n * K), nrow = n, ncol = K)
    tau <- matrix(rep(0, n * K), nrow = n, ncol = K)
    
    # Nested loop for actual imputation using rejection sampling
    for (i in 1:n) {
      index <- which(!is.na(Y[i, ]))
      for (k in 1:K) {
        check <- 0
        # Rejection sampling loop
        while (check < 1) {
          # Sample latent variable z from a normal distribution
          Z1[i, k] <- mu1[i] + sqrt(sigma1_2[i]) * rnorm(1)
          U <- runif(1)  # Uniform random number for acceptance
          # Calculate the likelihood ratio
          LP <- exp(-neg_logratio_binom(rbind(Y[i, index], N[index], 
                                              a[index], b[index]),
                                        Z1[i, k])) / gammax[i]
          if (is.nan(LP)) {LP <- 0}  # Handle NaN values
          if (U < LP) {  # Accept or reject the sampled z
            check <- 2
          }
        }
        # Calculate mu2 and tau based on the accepted z
        mu2[i, k] <- -((sigma_tau^2 - rho^2 * sigma_tau^2) * 
                         (sum(alpha[index]^2 * (logT10[i, index] - beta[index]))) - 
                         rho * sigma_tau * Z1[i, k]) / 
          ((sigma_tau^2 - rho^2 * sigma_tau^2) * A[i] + 1)
        tau[i, k] <- mu2[i, k] + sqrt(sigma2_2[i]) * rnorm(1)  # Sample tau
      }
    }
    
    # Returning the imputed values for theta and tau
    theta <- Z1
    imputes <- list(theta = theta, tau = tau, z.opt = zmax)
    return(imputes)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: mcem_count_bin
  # Purpose:
  #   Calculates the negative log-likelihood for the Binomial case during the M-step
  #   of the MCEM algorithm. This function is used to optimize parameters 'a' and 'b'.
  #
  # Inputs:
  #   - par     : Vector containing parameters 'a' and 'b' to be optimized.
  #   - Y       : Vector of observed counts for a single passage.
  #   - N       : Scalar number of words in the passage.
  #   - theta   : Vector of imputed latent variables for the passage.
  #   - K       : Number of imputations.
  #
  # Output:
  #   - nllh    : Negative log-likelihood value for the Binomial model.
  # ---------------------------------------------------------------------------
  mcem_count_bin <- function(par, Y, N, theta, K) {
    a <- par[1]
    b <- par[2]
    # Calculate success probabilities using the probit link
    success.prob <- pnorm(a * theta + b)
    # Replicate Y across K imputations
    Y.rep <- t(array(rep(Y, K), dim = c(length(Y), K)))
    # Calculate negative log-likelihood averaged over imputations
    nllh <- -sum(apply(dbinom(Y.rep, N, success.prob, log = TRUE), 2, mean))
    return(nllh)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: function_s12s22_to_min
  # Purpose:
  #   Defines the objective function for optimizing sigma_tau and rho based on
  #   the imputed latent variables. This function calculates a combined measure
  #   of covariance to be minimized.
  #
  # Inputs:
  #   - data : A list containing imputed values for theta and tau.
  #   - par  : Scalar value representing 's12', a component related to covariance.
  #
  # Output:
  #   - Fval : Combined objective value to be minimized.
  # ---------------------------------------------------------------------------
  function_s12s22_to_min <- function(data, par) {
    theta <- data$theta
    tau <- data$tau
    s12 <- par
    # Calculate s22 based on the current estimate of s12
    s22 <- s12^2 * (1 + mean(theta^2)) - 
      2 * s12 * mean(theta * tau) + 
      mean(tau^2)
    # Calculate components of the objective function
    F1 <- 0.5 * log(s22 - s12^2)
    F2 <- 0.5 / sqrt(s22 - s12^2) * 
      (s22 * mean(theta^2) - 
         2 * s12 * mean(theta * tau) + 
         mean(tau^2))
    Fval <- F1 + F2  # Combined objective value
    return(Fval)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: MCEM_algorithm_one_iteration
  # Purpose:
  #   Executes a single iteration of the MCEM algorithm, performing the E-step
  #   (imputation of latent variables) and the M-step (optimization of model parameters).
  #
  # Inputs:
  #   - Y         : Matrix of count data.
  #   - logT10    : Matrix of log-scale speed data.
  #   - N         : Vector of number of words in each passage.
  #   - J         : Number of passages.
  #   - ests      : Current estimates of model parameters.
  #   - K         : Number of imputations.
  #
  # Output:
  #   - EM.ests   : Updated estimates of model parameters after the iteration.
  # ---------------------------------------------------------------------------
  MCEM_algorithm_one_iteration <- function(Y, logT10, N, J, ests, K) {
    
    # -------------------------------------------------------------------------
    # E-step: Perform imputation of latent variables using the imputation_code_binom function
    # -------------------------------------------------------------------------
    EMimps <- imputation_code_binom(Y, logT10, N, ests$a, ests$b, 
                                    ests$alpha, ests$beta, 
                                    sqrt(ests$vartau), ests$rho, K)
    
    # Initializing vectors to store updated parameter estimates
    EM.a <- rep(NA, J)
    EM.b <- rep(NA, J)
    EM.alpha <- rep(NA, J)
    EM.beta <- rep(NA, J)
    
    # Extract imputed theta and tau matrices
    theta <- t(EMimps$theta)  # Transpose to align dimensions (J x K)
    tau <- t(EMimps$tau)      # Transpose to align dimensions (J x K)
    
    # -------------------------------------------------------------------------
    # M-step: Optimize parameters 'a', 'b', 'alpha', and 'beta' for each passage
    # -------------------------------------------------------------------------
    for (j in 1:J) {
      index <- which(!is.na(Y[, j]))  # Identify non-missing entries for passage j
      
      # ---------------------------------------------------------------------
      # Step 1: Optimize parameters 'a' and 'b' using the Binomial model
      # ---------------------------------------------------------------------
      bin.sol <- tryCatch({
        optim(
          c(ests$a[j], ests$b[j]), 
          mcem_count_bin, 
          Y = Y[index, j], 
          N = N[j], 
          theta = theta[, index],  # Extract relevant theta values for passage j
          K = K, 
          method = "BFGS"
        )
      }, 
      error = function(e) {
        # If BFGS fails, switch to Nelder-Mead optimization
        optim(
          c(ests$a[j], ests$b[j]), 
          mcem_count_bin, 
          Y = Y[index, j], 
          N = N[j], 
          theta = theta[, index], 
          K = K, 
          method = "Nelder-Mead"
        )
      })
      
      # Store optimized parameters 'a' and 'b'
      EM.a[j] <- bin.sol$par[1]
      EM.b[j] <- bin.sol$par[2]
      
      # ---------------------------------------------------------------------
      # Step 2: Update parameters 'beta' and 'alpha' based on imputations
      # ---------------------------------------------------------------------
      if (K > 1) {
        index <- which(!is.na(logT10[, j]))
        # Update 'beta' as the mean of logT10 plus the mean of tau across imputations
        EM.beta[j] <- mean(logT10[index, j] + apply(tau[, index, drop = FALSE], 2, mean))
        # Create a matrix for logT10 values
        logT10_vector <- t(matrix(logT10[index, j], nrow = 1))
        # Calculate residuals
        result <- sweep(tau[, index, drop = FALSE], 2, logT10_vector, "+") - EM.beta[j]
        # Calculate inverse of alpha squared
        alp2_inv <- mean(apply(result^2, 2, mean))
        # Update 'alpha' based on residuals
        EM.alpha[j] <- 1 / sqrt(alp2_inv)
      } else {
        index <- which(!is.na(logT10[, j]))
        # Update 'beta' as the mean of logT10 plus tau
        EM.beta[j] <- mean(logT10[index, j] + tau[, index])
        # Create a matrix for logT10 values
        logT10_vector <- t(matrix(logT10[index, j], nrow = 1))
        # Calculate residuals
        result <- tau[, index] + logT10_vector - EM.beta[j]
        # Calculate inverse of alpha squared
        alp2_inv <- mean(result^2)
        # Update 'alpha' based on residuals
        EM.alpha[j] <- 1 / sqrt(alp2_inv)
      }
      
    }  # End of optimization loop for passages
    
    # -------------------------------------------------------------------------
    # Step 3: Optimize sigma_tau and rho based on the imputed latent variables
    # -------------------------------------------------------------------------
    s12.find <- optim(ests$rho * sqrt(ests$vartau), 
                      fn = function_s12s22_to_min, 
                      method = "BFGS", 
                      data = EMimps)
    s12.min <- s12.find$par
    # Calculate updated variance parameter 'vartau'
    EM.vartau <- s12.min^2 * (1 + mean(EMimps$theta^2)) - 
      2 * s12.min * mean(EMimps$theta * EMimps$tau) + 
      mean(EMimps$tau^2)
    # Update intra-class correlation 'rho'
    EM.rho <- s12.min / sqrt(EM.vartau)
    
    # -------------------------------------------------------------------------
    # Compiling and Returning the Updated Parameter Estimates
    # -------------------------------------------------------------------------
    EM.ests <- list(
      a = EM.a, 
      b = EM.b, 
      alpha = EM.alpha, 
      beta = EM.beta, 
      vartau = EM.vartau, 
      rho = EM.rho
    )
    return(EM.ests)  # Return the updated estimates
  }
  
  # ---------------------------------------------------------------------------
  # Handling the k.in and reps.in Parameters for Multiple Imputations
  # Purpose:
  #   Determines the total number of imputations based on k.in and reps.in, which
  #   control the number of iterations and repetitions for each MCEM run.
  # ---------------------------------------------------------------------------
  nK <- length(k.in)
  if (missing(ests.in)) {
    # If initial estimates are not provided, compute MOM estimates
    ests.in <- mom.binom(Y, logT10, N, J)
  }
  
  # ---------------------------------------------------------------------------
  # Checking and Adjusting Infinite Alpha Values
  # Purpose:
  #   Ensures that the precision parameter 'alpha' does not contain infinite values
  #   by capping them at a reasonable multiple of the maximum finite alpha value.
  # ---------------------------------------------------------------------------
  alpha.check <- ests.in$alpha
  infIndex0 <- which(!is.infinite(alpha.check))
  infIndex1 <- which(is.infinite(alpha.check))
  if (length(infIndex0) > 0 && length(infIndex1) > 0) {
    alpha.check[infIndex1] <- max(alpha.check[infIndex0]) * 10  # Cap infinite values
    ests.in$alpha <- alpha.check  # Update 'alpha' estimates
  }
  
  # ---------------------------------------------------------------------------
  # Setting Up Initial Imputation Run
  # Purpose:
  #   Initializes the total number of iterations (K) for the MCEM runs based on
  #   the provided k.in and reps.in parameters.
  # ---------------------------------------------------------------------------
  n <- dim(Y)[1]
  total.K <- rep(k.in[1], reps.in[1])  # Initialize with the first k.in and reps.in
  
  # Handling multiple imputation configurations if k.in is a vector
  if (nK > 1) {
    for (jj in 2:nK) {
      total.K <- c(total.K, rep(k.in[jj], reps.in[jj]))
    }
  }
  JJ <- length(total.K)  # Total number of MCEM runs
  
  # ---------------------------------------------------------------------------
  # Initializing Storage Matrices for Parameter Estimates
  # Purpose:
  #   Creates matrices to store parameter estimates ('a', 'b', 'alpha', 'beta',
  #   'vartau', 'rho') across all MCEM iterations.
  # ---------------------------------------------------------------------------
  a.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  alpha.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  b.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  beta.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  rho.store <- matrix(rep(0, (JJ + 1)), nrow = (JJ + 1))
  vartau.store <- matrix(rep(0, (JJ + 1)), nrow = (JJ + 1))
  
  # ---------------------------------------------------------------------------
  # Main Loop for the MCEM Algorithm Iterations
  # Purpose:
  #   Executes the MCEM algorithm for each specified run, performing imputation
  #   and parameter optimization steps, and storing the results.
  # ---------------------------------------------------------------------------
  for (jj in 1:JJ) {
    # Store current parameter estimates before updating
    a.store[jj, ] <- ests.in$a
    b.store[jj, ] <- ests.in$b
    alpha.store[jj, ] <- ests.in$alpha
    beta.store[jj, ] <- ests.in$beta
    vartau.store[jj] <- ests.in$vartau
    rho.store[jj] <- ests.in$rho
    
    # Perform one iteration of the MCEM algorithm
    EM.iter <- MCEM_algorithm_one_iteration(Y, logT10, 
                                            N, J, ests.in, 
                                            total.K[jj])
    ests.in <- EM.iter  # Update parameter estimates
  }
  
  # ---------------------------------------------------------------------------
  # Finalizing Parameter Estimates
  # Purpose:
  #   After completing all MCEM iterations, store the final parameter estimates.
  # ---------------------------------------------------------------------------
  ests.in <- EM.iter
  a.store[JJ + 1, ] <- ests.in$a
  b.store[JJ + 1, ] <- ests.in$b
  alpha.store[JJ + 1, ] <- ests.in$alpha
  beta.store[JJ + 1, ] <- ests.in$beta
  vartau.store[JJ + 1] <- ests.in$vartau
  rho.store[JJ + 1] <- ests.in$rho
  # Final parameter estimates can be accessed from the last row of the storage matrices.
  
  # ---------------------------------------------------------------------------
  # Compiling and Returning the Final Estimates
  # ---------------------------------------------------------------------------
  MCEM.ests <- list(
    a = a.store,         # Matrix of 'a' estimates across iterations
    b = b.store,         # Matrix of 'b' estimates across iterations
    alpha = alpha.store, # Matrix of 'alpha' estimates across iterations
    beta = beta.store,   # Matrix of 'beta' estimates across iterations
    numwords.p = N,      # Vector of number of words per passage
    vartau = vartau.store,# Vector of 'vartau' estimates across iterations
    rho = rho.store      # Vector of 'rho' estimates across iterations
  )
  
  return(MCEM.ests)  # Return the list of estimated parameters
}

# =============================================================================
# End of mcem_binom.R
# =============================================================================