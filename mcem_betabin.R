# =============================================================================
# File: mcem_betabin.R
# Author: Cornelis J Potgieter
#
# ---------------------------------------------------------------------------
# Description:
#   This script implements the `run.mcem.betabin` function, which implements
#   the MCEM algorithm for parameter estimation in the Beta-Binomial 
#   variation of the count-time model. The function iteratively
#   imputes latent variables and updates model parameters to approximate the
#   maximum likelihood estimates.
#
#   - `run.mcem.betabin`: Main function executing the MCEM algorithm.
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
#                   - nu      : Estimated dispersion parameter 'nu' for each passage.
#                   - alpha   : Estimated precision parameter 'alpha' for each passage.
#                   - beta    : Estimated parameter 'beta' for each passage.
#                   - numwords.p : Vector of number of words per passage (input N).
#                   - vartau  : Estimated variance parameter 'vartau' across iterations.
#                   - rho     : Estimated intra-class correlation 'rho' across iterations.
#
# Dependencies:
#   - `mvtnorm` package
#
# Usage:
#   1. Load the `mvtnorm` and `extraDistr` packages
#
#   2. Source this script to make the `run.mcem.betabin` function available:
#        source("path/to/mcem_betabin.R")
#
#   3. (Optional) Compute initial parameter estimates using MOM functions:
#        ests.in <- mom.betabin(Y, logT10, N, J)
#
#   4. Run the MCEM algorithm:
#        MCEM_results <- run.mcem.betabin(Y, logT10, N, J, k.in = 5, reps.in = 2, ests.in)
#
#   6. Access the estimated parameters:
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
#   - The function uses `dbbinom` for Beta-Binomial distribution in `extraDistr`
#
# =============================================================================

# -------------------------------
# Function Definitions
# -------------------------------

# Main function implementing the MCEM algorithm
# Parameters: (a, b, nu, alpha, beta, sigma_tau, rho)
# Y: Matrix of count data
# logT10: Matrix of log-scale speed data
# N: Vector with number of words in each passage
# J: Number of passages
# k.in: Number of Monte Carlo imputations
# reps.in: Number of MCEM iterations
# Can have k.in = 1 and reps.in = 50 OR k.in = c(1,10) and reps.in = (50,5)
# ests.in: List with initial values for model parameters
#
# MCEM specificaion examples:
#   - k.in = 1 and reps.in = 50
#   - k.in = c(1, 10) and reps.in = c(50, 5)

run.mcem.betabin <- function(Y, logT10, N, J, k.in = 5, reps.in = 2, ests.in) {
  
  # ---------------------------------------------------------------------------
  # Internal Function: neg_logratio_betabin
  # Purpose:
  #   Computes the negative log-likelihood for the Beta-Binomial model. This function
  #   is utilized in the rejection sampling step of the MCEM sampler to evaluate
  #   the suitability of sampled latent variables.
  #
  # Inputs:
  #   - data : A matrix where each row corresponds to a different parameter:
  #           Row 1: Y (counts)
  #           Row 2: N (number of words)
  #           Row 3: a (parameter 'a')
  #           Row 4: b (parameter 'b')
  #           Row 5: nu (dispersion parameter 'nu')
  #   - par  : Scalar latent variable 'z' to be evaluated.
  #
  # Output:
  #   - negLP : Negative log-likelihood value for the given parameters and latent variable.
  # ---------------------------------------------------------------------------
  neg_logratio_betabin <- function(data, par) {
    # Extracting data elements
    Y <- data[1, ]
    N <- data[2, ]
    a <- data[3, ]
    b <- data[4, ]
    nu <- data[5, ]
    z <- par
    
    # Calculating success probabilities using the probit link function
    success.prob <- pnorm(a * z + b)
    
    # Calculating Beta-Binomial parameters alpha.bb and beta.bb
    alpha.bb <- success.prob * (N - nu) / (nu - 1)
    beta.bb <- alpha.bb * (1 - success.prob) / success.prob
    
    # Handling binomial and beta-binomial cases separately based on nu
    index1 <- which(nu < 1.0001)    # Binomial case where nu ≈ 1
    index2 <- which(nu >= 1.0001)   # Beta-Binomial case where nu > 1
    
    # Calculating negative log-likelihood for binomial cases
    negLP1 <- -sum(dbinom(Y[index1], N[index1], success.prob[index1], log = TRUE))
    
    # Calculating negative log-likelihood for beta-binomial cases
    negLP2 <- -sum(dbbinom(Y[index2], N[index2], alpha.bb[index2], 
                           beta.bb[index2], log = TRUE))
    
    # Summing up the negative log-likelihoods
    negLP <- negLP1 + negLP2
    
    # Handling infinite negative log-likelihood values by assigning a large constant
    if (is.infinite(negLP)) {
      negLP <- 1e20  # Replace with a suitably large constant to penalize
    }
    
    return(negLP)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: gamma_multiplier_betabin
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
  #   - nu  : Vector of dispersion parameters 'nu' corresponding to Y.
  #
  # Output:
  #   - gamval : A list containing 'zmax' and 'gammax' values for rejection sampling.
  # ---------------------------------------------------------------------------
  gamma_multiplier_betabin <- function(Y, N, a, b, nu) {
    # Check for a special case where the count data is perfect (no variability)
    if (sum(N - Y) < 0.01) {
      gammax <- 1
      zmax <- Inf
    } else {
      # Optimizing the negative log-ratio over the interval [-7, 7]
      suppressWarnings({
        zgam <- optimize(neg_logratio_betabin, 
                         data = rbind(Y, N, a, b, nu), 
                         interval = c(-7, 7))
      })
      zmax <- zgam$minimum
      gammax <- exp(-zgam$objective)  # Calculating gamma multiplier
    }
    
    # Returning the calculated gamma multiplier and optimal z
    gamval <- list(zmax = zmax, gammax = gammax)
    return(gamval)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: imputation_code
  # Purpose:
  #   Performs the imputation step of the MCEM algorithm by sampling latent
  #   variables using rejection sampling based on the current parameter estimates.
  #
  # Inputs:
  #   - Y         : Matrix of count data.
  #   - logT10    : Matrix of log-scale speed data.
  #   - N         : Vector of number of words in each passage.
  #   - a, b, nu  : Vectors of current parameter estimates.
  #   - alpha, beta: Vectors of current parameter estimates.
  #   - sigma_tau : Current estimate of sigma_tau (standard deviation).
  #   - rho       : Current estimate of intra-class correlation.
  #   - K         : Number of imputations to perform.
  #
  # Output:
  #   - imputes : A list containing imputed values for theta and tau, and the
  #              optimal z values.
  # ---------------------------------------------------------------------------
  imputation_code <- function(Y, logT10, N, a, b, nu, alpha, beta, sigma_tau, rho, K) {
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
      gamval <- gamma_multiplier_betabin(Y[i, index], N[index], 
                                         a[index], b[index], nu[index])
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
          LP <- exp(-neg_logratio_betabin(rbind(Y[i, index], N[index], 
                                                a[index], b[index], 
                                                nu[index]), Z1[i, k])) / gammax[i]
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
  #   Calculates the negative log-likelihood for the binomial case during the M-step
  #   of the MCEM algorithm. This function is used to optimize parameters 'a' and 'b'
  #   when the dispersion parameter 'nu' is fixed at 1 (binomial case).
  #
  # Inputs:
  #   - par     : Vector containing parameters 'a' and 'b' to be optimized.
  #   - Y       : Vector of observed counts for a single passage.
  #   - N       : Scalar number of words in the passage.
  #   - theta   : Vector of imputed latent variables for the passage.
  #   - K       : Number of imputations.
  #
  # Output:
  #   - nllh    : Negative log-likelihood value for the binomial model.
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
  # Internal Function: mcem_count_betabin_fixednu
  # Purpose:
  #   Calculates the negative log-likelihood for the beta-binomial case with a fixed
  #   dispersion parameter 'nu'. This function is used to optimize parameters 'a' and 'b'
  #   when 'nu' is held constant.
  #
  # Inputs:
  #   - par     : Vector containing parameters 'a' and 'b' to be optimized.
  #   - nu      : Fixed dispersion parameter 'nu' for the passage.
  #   - Y       : Vector of observed counts for a single passage.
  #   - N       : Scalar number of words in the passage.
  #   - theta   : Vector of imputed latent variables for the passage.
  #   - K       : Number of imputations.
  #
  # Output:
  #   - nllh    : Negative log-likelihood value for the beta-binomial model.
  # ---------------------------------------------------------------------------
  mcem_count_betabin_fixednu <- function(par, nu, Y, N, theta, K) {
    a <- par[1]
    b <- par[2]
    # Calculate success probabilities using the probit link
    success.prob <- pnorm(a * theta + b)
    # Calculate Beta-Binomial parameters alpha and beta
    alpha <- success.prob * (N - nu) / (nu - 1)
    beta <- alpha * (1 - success.prob) / success.prob
    # Replicate Y across K imputations
    Y.rep <- t(array(rep(Y, K), dim = c(length(Y), K)))
    # Calculate log-likelihood for beta-binomial distribution
    log.llh.calcs <- array(dbbinom(Y.rep, N, alpha, beta, log = TRUE), 
                           dim = c(K, length(Y)))
    # Calculate negative log-likelihood averaged over imputations
    nllh <- -sum(apply(log.llh.calcs, 2, mean))
    return(nllh)
  }
  
  # ---------------------------------------------------------------------------
  # Internal Function: mcem_count_betabin_nu
  # Purpose:
  #   Optimizes the dispersion parameter 'nu' for the beta-binomial model by
  #   minimizing the negative log-likelihood. This function is used within the
  #   parameter optimization loop during the M-step.
  #
  # Inputs:
  #   - par     : Scalar value for 'nu' to be optimized.
  #   - a.in    : Current estimate of parameter 'a'.
  #   - b.in    : Current estimate of parameter 'b'.
  #   - Y       : Vector of observed counts for a single passage.
  #   - N       : Scalar number of words in the passage.
  #   - theta   : Vector of imputed latent variables for the passage.
  #   - K       : Number of imputations.
  #
  # Output:
  #   - nu.run$value : Negative log-likelihood value for the optimized 'nu'.
  # ---------------------------------------------------------------------------
  mcem_count_betabin_nu <- function(par, a.in, b.in, Y, N, theta, K) {
    nu <- par
    # Attempt to optimize parameters 'a' and 'b' using BFGS method
    nu.run <- tryCatch(
      {
        optim(
          c(a.in, b.in),
          fn = mcem_count_betabin_fixednu,
          nu = nu,
          Y = Y,
          N = N,
          theta = theta,
          K = K,
          method = "BFGS"
        )
      },
      error = function(e1) {
        # If an error occurs, switch to Nelder-Mead method
        optim(
          c(a.in, b.in),
          fn = mcem_count_betabin_fixednu,
          nu = nu,
          Y = Y,
          N = N,
          theta = theta,
          K = K,
          method = "Nelder-Mead"
        )
      }
    )
    
    return(nu.run$value)  # Return the optimized negative log-likelihood
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
    # E-step: Perform imputation of latent variables using the imputation_code function
    # -------------------------------------------------------------------------
    EMimps <- imputation_code(Y, logT10, N, ests$a, ests$b, 
                              ests$nu, ests$alpha, ests$beta, 
                              sqrt(ests$vartau), ests$rho, K)
    
    # Initializing vectors to store updated parameter estimates
    EM.a <- rep(NA, J)
    EM.b <- rep(NA, J)
    EM.nu <- rep(NA, J)
    EM.alpha <- rep(NA, J)
    EM.beta <- rep(NA, J)
    
    # Extract imputed theta and tau matrices
    theta <- t(EMimps$theta)  # Transpose to align dimensions (J x K)
    tau <- t(EMimps$tau)      # Transpose to align dimensions (J x K)
    
    # -------------------------------------------------------------------------
    # M-step: Optimize parameters 'a', 'b', 'nu', 'alpha', and 'beta' for each passage
    # -------------------------------------------------------------------------
    for (j in 1:J) {
      index <- which(!is.na(Y[, j]))  # Identify non-missing entries for passage j
      
      # ---------------------------------------------------------------------
      # Step 1: Optimize parameters 'a' and 'b' using the binomial model
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
      
      # Extract initial parameter estimates from binomial optimization
      initial_a <- bin.sol$par[1]
      initial_b <- bin.sol$par[2]
      
      # ---------------------------------------------------------------------
      # Step 2: Optimize dispersion parameter 'nu' using the beta-binomial model
      # ---------------------------------------------------------------------
      nu.sol <- tryCatch(
        {
          optimize(
            mcem_count_betabin_nu,
            a.in = initial_a,
            b.in = initial_b,
            Y = Y[index, j],
            N = N[j],
            theta = theta[, index],
            K = K,
            interval = c(1.0001, N[j] - 0.0001)
          )
        },
        error = function(e) {
          # If optimization fails, iteratively adjust 'a' and 'b' until successful
          while(TRUE) {
            U <- runif(1)  # Generate a uniform random number
            initial_a <- U * bin.sol$par[1]
            initial_b <- U * bin.sol$par[2]
            
            # Retry optimization with updated 'a' and 'b'
            nu.sol <- tryCatch(
              {
                optimize(
                  mcem_count_betabin_nu,
                  a.in = initial_a,
                  b.in = initial_b,
                  Y = Y[index, j],
                  N = N[j],
                  theta = theta[, index],
                  K = K,
                  interval = c(1.0001, N[j] - 0.0001)
                )
              },
              error = function(e) {
                NULL  # Return NULL if optimization fails again
              }
            )
            
            # Break the loop if a solution is found
            if(!is.null(nu.sol)) {
              break
            }
          }
          nu.sol
        }
      )
      
      # Store current 'a' and 'b' estimates
      a_store <- initial_a
      b_store <- initial_b
      
      # ---------------------------------------------------------------------
      # Step 3: Check and adjust the negative log-likelihood
      # ---------------------------------------------------------------------
      result <- mcem_count_betabin_fixednu(c(initial_a, initial_b), 
                                           nu = nu.sol$minimum, 
                                           Y = Y[index, j], 
                                           N = N[j], 
                                           theta = theta[, index], 
                                           K = K)
      
      # Repeat optimization until a valid (non-NaN) result is obtained
      while (is.nan(result)) {
        U <- runif(1)  # Generate a new uniform random number
        backup_a <- U * a_store
        backup_b <- U * b_store
        # Update 'a' and 'b' with backup values
        result <- mcem_count_betabin_fixednu(c(backup_a, backup_b), 
                                             nu = nu.sol$minimum, 
                                             Y = Y[index, j], 
                                             N = N[j], 
                                             theta = theta[, index], 
                                             K = K)
        initial_a <- backup_a
        initial_b <- backup_b
      }
      
      # ---------------------------------------------------------------------
      # Step 4: Optimize parameters 'a' and 'b' using the beta-binomial model
      # ---------------------------------------------------------------------
      betabin.sol <- tryCatch(
        {
          optim(
            c(initial_a, initial_b),
            mcem_count_betabin_fixednu,
            nu = nu.sol$minimum,
            Y = Y[index, j],
            N = N[j],
            theta = theta[, index],
            K = K,
            method = "BFGS"
          )
        },
        error = function(e1) {
          # If BFGS fails, switch to Nelder-Mead optimization
          optim(
            c(initial_a, initial_b),
            mcem_count_betabin_fixednu,
            nu = nu.sol$minimum,
            Y = Y[index, j],
            N = N[j],
            theta = theta[, index],
            K = K,
            method = "Nelder-Mead"
          )
        }
      )
      
      # ---------------------------------------------------------------------
      # Step 5: Select the model (binomial or beta-binomial) with lower NLLH
      # ---------------------------------------------------------------------
      if (bin.sol$value < betabin.sol$value) {
        # Binomial model provides a better fit
        EM.a[j] <- bin.sol$par[1]
        EM.b[j] <- bin.sol$par[2]
        EM.nu[j] <- 1  # Set dispersion parameter 'nu' to 1 for binomial
      } else {
        # Beta-Binomial model provides a better fit
        EM.a[j] <- betabin.sol$par[1]
        EM.b[j] <- betabin.sol$par[2]
        EM.nu[j] <- nu.sol$minimum  # Update 'nu' with optimized value
      }
      
      # ---------------------------------------------------------------------
      # Step 6: Update parameters 'beta' and 'alpha' based on imputations
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
    # Step 7: Optimize sigma_tau and rho based on the imputed latent variables
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
    # Step 8: Compile and return the updated parameter estimates
    # -------------------------------------------------------------------------
    EM.ests <- list(
      a = EM.a, 
      b = EM.b, 
      nu = EM.nu,
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
    ests.in <- mom.betabin(Y, logT10, N, J)
    # Alternatively, you can use binomial MOM estimates by uncommenting the next line:
    # ests.in <- mom.binom(Y, logT10, N, J) 
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
  alpha.check[infIndex1] <- max(alpha.check[infIndex0]) * 10  # Cap infinite values
  ests.in$alpha <- alpha.check  # Update 'alpha' estimates
  
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
  #   Creates matrices to store parameter estimates ('a', 'b', 'nu', 'alpha', 'beta',
  #   'vartau', 'rho') across all MCEM iterations.
  # ---------------------------------------------------------------------------
  a.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  alpha.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  b.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  nu.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
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
    nu.store[jj, ] <- ests.in$nu
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
  nu.store[JJ + 1, ] <- ests.in$nu
  alpha.store[JJ + 1, ] <- ests.in$alpha
  beta.store[JJ + 1, ] <- ests.in$beta
  vartau.store[JJ + 1] <- ests.in$vartau
  rho.store[JJ + 1] <- ests.in$rho
  # Final parameter estimates.
  #mean_a = a.store[JJ,]
  #mean_b = b.store[JJ,]
  #mean_alpha = alpha.store[JJ,]
  #mean_beta = beta.store[JJ,]
  
  # ---------------------------------------------------------------------------
  # Compiling and Returning the Final Estimates
  # ---------------------------------------------------------------------------
  MCEM.ests <- list(
    a = a.store,         # Matrix of 'a' estimates across iterations
    b = b.store,         # Matrix of 'b' estimates across iterations
    nu = nu.store,       # Matrix of 'nu' estimates across iterations
    alpha = alpha.store, # Matrix of 'alpha' estimates across iterations
    beta = beta.store,   # Matrix of 'beta' estimates across iterations
    numwords.p = N,      # Vector of number of words per passage
    vartau = vartau.store,# Vector of 'vartau' estimates across iterations
    rho = rho.store      # Vector of 'rho' estimates across iterations
  )
  
  return(MCEM.ests)  # Return the list of estimated parameters
}
