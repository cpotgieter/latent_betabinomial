# =============================================================================
# File: mom_estimation.R
# Author: Cornelis J Potgieter
#
# ---------------------------------------------------------------------------
# Description:
#   This script implements two functions, `mom.binom` and `mom.betabin`, which
#   perform Method of Moments (MOM) estimation for the Binomial and Beta-Binomial
#   count models in the joint count-time model as described in the paper 
#   "Joint Analysis of Dispersed Count-Time Data 
#   Using a Bivariate Latent Factor Model".
#
#   - `mom.binom`: Estimates MOM parameters for the Binomial specification.
#   - `mom.betabin`: Estimates MOM parameters for the Beta-Binomial specification.
#
# Inputs:
#   Both functions require the following inputs:
#     - Y       : Matrix of count data (n x J), where n is the number of observations
#                 and J is the number of items.
#     - logT10  : Matrix of log-time data (n x J), corresponding to the counts in Y.
#     - N       : Vector of item lengths (length J).
#     - J       : Integer representing the number of items.
#
# Outputs:
#   Each function returns a list containing the estimated parameters:
#     - For `mom.binom`:
#         - a        : Estimated parameter 'a' for each item.
#         - b        : Estimated parameter 'b' for each item.
#         - alpha    : Estimated precision parameter 'alpha' for each item.
#         - beta     : Estimated parameter 'beta' for each item.
#         - vartau   : Estimated variance parameter 'vartau'.
#         - rho      : Estimated intra-class correlation 'rho'.
#
#     - For `mom.betabin`:
#         - a          : Estimated parameter 'a' for each item.
#         - b          : Estimated parameter 'b' for each item.
#         - nu         : Estimated parameter 'nu' for each item.
#         - alpha      : Estimated precision parameter 'alpha' for each item.
#         - beta       : Estimated parameter 'beta' for each item.
#         - vartau     : Estimated variance parameter 'vartau'.
#         - rho        : Estimated intra-class correlation 'rho'.
#
# Dependencies:
#   - `mvtnorm` package: Required for multivariate normal distribution functions.
#
# Usage:
#   1. Load the `mvtnorm` package: library(mvtnorm)
#
#   2. Source this script to make the functions available in your R session:
#        source("path/to/mom_estimation.R")
#
#   3. Call the functions with appropriate arguments:
#        binom_estimates <- mom.binom(Y, logT10, N, J)
#        betabin_estimates <- mom.betabin(Y, logT10, N, J)
#
# Notes:
#   - Ensure that the input matrices `Y` and `logT10` have the same dimensions.
#   - Set missing data to NA before passing matrices to the functions.
#   - The functions include internal helper functions for optimization and parameter
#     estimation based on the MOM approach.
#   - For large datasets, consider the computational time due to optimization steps.
#
# =============================================================================

# -------------------------------
# Function Definitions
# -------------------------------

mom.binom <- function(Y, logT10, N, J) {
  
  # ---------------------------------------------------------------------------
  # Function: mom.binom
  # Purpose:
  #   Performs MOM estimation under the Binomial model for count data.
  #   Estimates (a, b, alpha, beta, vartau, rho) based on the observed
  #   counts and log-times.
  #
  # Inputs:
  #   Y     - Matrix of count data (n x J).
  #   logT10- Matrix of log-time data (n x J), corresponding to Y.
  #   N     - Vector of item lengths (length J).
  #   J     - Integer representing the number of items.
  #
  # Outputs:
  #   A list containing:
  #     - a        : Estimated parameter 'a' for each item.
  #     - b        : Estimated parameter 'b' for each item.
  #     - alpha    : Estimated parameter 'alpha' for each item.
  #     - beta     : Estimated parameter 'beta' for each item.
  #     - vartau   : Estimated latent variance parameter 'vartau' (single numeric value).
  #     - rho      : Estimated latent correlation 'rho' (single numeric value).
  # ---------------------------------------------------------------------------
  
  # -------------------------------
  # Internal Function: evaluate_rho
  # Purpose:
  #   Computes the Mean Squared Distance (MD) between the observed R and the
  #   cumulative distribution function of a bivariate normal distribution.
  #
  # Inputs:
  #   data - Numeric vector containing two elements:
  #          - R: Computed statistic from the data.
  #          - Q: Quantile value.
  #   par  - Scalar parameter 'r' to be optimized.
  #
  # Output:
  #   MD  - Mean Squared Distance between observed R and the bivariate normal CDF.
  # -------------------------------
  evaluate_rho <- function(data, par) {
    R <- data[1]
    Q <- data[2]
    r <- par
    sigma <- matrix(c(1, r, r, 1), ncol = 2)  # Correlation matrix with correlation 'r'
    # Compute the cumulative probability up to (Q, Q) in a bivariate normal distribution
    MD <- (R - pmvnorm(upper = c(Q, Q), mean = rep(0, 2), sigma = sigma))^2
    return(MD)
  }
  
  # -------------------------------
  # Internal Function: reading_data_parms_MOM
  # Purpose:
  #   MOM parameter estimates of 'a' and 'b' using observed counts and item lengths.
  #
  # Inputs:
  #   Y - Vector of counts for a single item (length n).
  #   N - Scalar representing the item length.
  #
  # Output:
  #   parms - List containing estimated parameters 'a' and 'b'.
  # -------------------------------
  reading_data_parms_MOM <- function(Y, N) {
    Y.bar <- mean(Y)  # Sample mean of counts
    S2 <- mean((Y - Y.bar)^2)  # Sample variance of counts
    
    # Compute Q based on Y.bar and N to handle cases where Y.bar >= N
    if (Y.bar < N) {
      Q <- qnorm(Y.bar / N)
    } else {
      Q <- qnorm((Y.bar + 0.05) / (N + 0.1)) # "Continuity Correction
    }
    
    # Compute R based on sample variance and mean
    R <- (S2 + Y.bar^2 - Y.bar) / (N * (N - 1))
    
    parms.in <- 0.5  # Initial guess for 'r'
    
    # Optimize 'r' to minimize the Mean Squared Distance (MD)
    rho.min <- optim(
      parms.in,
      fn = evaluate_rho,
      method = "Brent",
      data = c(R, Q),
      lower = 0,
      upper = 1
    )
    
    # Estimate parameter 'a' from the optimized 'r'
    a.out <- sqrt(rho.min$par / (1 - rho.min$par))
    
    # Estimate parameter 'b' based on 'Q' and 'a.out'
    b.out <- Q * sqrt(1 + a.out^2)
    
    # Return the estimated parameters
    parms <- list(a = a.out, b = b.out)
    return(parms)
  }
  
  # -------------------------------
  # Method of Moments (MOM) Estimation
  # -------------------------------
  
  n <- dim(Y)[1]  # Number of observations
  # Covariance matrix of log-time data
  nancov_T <- cov(logT10, 
                  use = "pairwise.complete.obs")  
  
  # Initialize matrices to store parameters 'a' and 'b'
  a.out <- matrix(nrow = 1, ncol = J)
  b.out <- matrix(nrow = 1, ncol = J)
  
  # Estimate parameters 'a' and 'b' for each item using MOM
  for (j in 1:J) {
    ab.out <- reading_data_parms_MOM(na.omit(Y[, j]), N[j])
    a.out[j] <- ab.out$a
    b.out[j] <- ab.out$b
  }
  
  # Estimate 'beta' as the mean of log-time data
  beta.out <- apply(logT10, 2, mean, na.rm = TRUE)  
  
  # Initialize variables for covariance adjustments
  count <- 1
  alpha.inv <- numeric(0)
  s22.comp <- numeric(0)
  
  # Compute inverse precision components and covariance adjustments
  for (j in 1:J) {
    if (j == J) {
      break  # Exit loop if j reaches the number of items
    }
    for (k in (j + 1):J) {
      alpha.inv <- rbind(
        alpha.inv,
        matrix(rep(0, J + 1), nrow = 1, ncol = (J + 1))
      )
      alpha.inv[count, j] <- max(0, nancov_T[j, j] - nancov_T[j, k])
      alpha.inv[count, k] <- max(0, nancov_T[k, k] - nancov_T[j, k])
      s22.comp <- rbind(s22.comp, nancov_T[j, k])
      count <- count + 1
    }
  }
  
  # Estimate 'alpha' based on inverse precision components
  alpha.inv_out <- apply(
    alpha.inv[, 1:J],
    2,
    sum,
    na.rm = TRUE
  ) / (J - apply(is.na(alpha.inv[, 1:J]), 2, sum) - 1)
  alpha.out <- 1 / sqrt(alpha.inv_out)  # Final 'alpha' estimates
  
  # Estimate 'vartau' as the mean of covariance components
  vartau.out <- max(0, mean(s22.comp, na.rm = TRUE))
  
  # Initialize matrix to store s12 calculations
  s12.calc <- matrix(rep(0, J), nrow = J)
  
  # Compute 's12' for each item based on covariance with log-time data
  for (j in 1:J) {
    CV <- cov(Y[, j], logT10[, j], use = "pairwise.complete.obs")  # Covariance between count and log-time
    s12.calc[j] <- -CV / N[j] * sqrt(2 * pi) * sqrt(a.out[j]^2 + 1) / a.out[j] *
      exp(0.5 * b.out[j]^2 / (1 + a.out[j]^2))
    s12.calc[j] <- max(-sqrt(vartau.out), min(sqrt(vartau.out), s12.calc[j]))  # Constrain 's12' within bounds
  }
  
  # Average 's12' across all items
  s12.out <- mean(s12.calc)
  
  # Compute intra-class correlation 'rho' based on 'vartau'
  if (vartau.out > 0) {
    s12.cor <- s12.out / sqrt(vartau.out)
  } else {
    s12.cor <- 0
  }
  
  # Compile all MOM estimates into a list
  MOMests <- list(
    a = as.numeric(a.out), 
    b = as.numeric(b.out),
    alpha = as.numeric(alpha.out), 
    beta = as.numeric(beta.out), 
    vartau = as.numeric(vartau.out), 
    rho = as.numeric(s12.cor)
  )
  
  # Return the MOM estimates
  return(MOMests)
}

mom.betabin <- function(Y, logT10, N, J) {
  
  # ---------------------------------------------------------------------------
  # Function: mom.betabinom
  # Purpose:
  #   Performs MOM estimation under the Beta-Binomial model for count data.
  #   Estimates (a, b, nu, alpha, beta, vartau, rho) based on the observed
  #   counts and log-times.
  #
  # Inputs:
  #   Y     - Matrix of count data (n x J).
  #   logT10- Matrix of log-time data (n x J), corresponding to Y.
  #   N     - Vector of item lengths (length J).
  #   J     - Integer representing the number of items.
  #
  # Outputs:
  #   A list containing:
  #     - a        : Estimated parameter 'a' for each item.
  #     - b        : Estimated parameter 'b' for each item.
  #     - nu       : Estimated parameter 'nu' for each item.
  #     - alpha    : Estimated parameter 'alpha' for each item.
  #     - beta     : Estimated parameter 'beta' for each item.
  #     - vartau   : Estimated latent variance parameter 'vartau' (single numeric value).
  #     - rho      : Estimated latent correlation 'rho' (single numeric value).  
  # ---------------------------------------------------------------------------
  # Internal Function: bivnorm.min
  # Purpose:
  #   Computes the Mean Squared Distance (MD) between the observed mu.biv and the
  #   cumulative distribution function of a bivariate normal distribution with
  #   correlation 'r'.
  #
  # Inputs:
  #   qr      - Scalar quantile value.
  #   mu.biv  - Scalar representing the observed statistic mu.biv.
  #   eta.marg- Vector of two elements representing the marginal quantiles.
  #
  # Output:
  #   crit    - Mean Squared Distance between observed mu.biv and the bivariate normal CDF.
  # -------------------------------
  bivnorm.min <- function(qr, mu.biv, eta.marg) {
    r <- pnorm(qr)  # Convert quantile to probability using standard normal CDF
    bivnorm <- pmvnorm(
      upper = eta.marg,
      mean = rep(0, 2),
      corr = matrix(c(1, r, r, 1), nrow = 2),
      keepAttr = FALSE
    )
    crit <- (mu.biv - bivnorm)^2  # Compute squared distance
    return(crit)
  }
  
  # Vectorize the bivnorm.min function for vectorized operations
  bivnorm.min_vectorized <- Vectorize(bivnorm.min, vectorize.args = "qr")
  
  # -------------------------------
  # Internal Function: reading_data_parms_MOM
  # Purpose:
  #   Estimates parameters 'a', 'b', and 'nu' using the Method of Moments based on
  #   the observed counts and item lengths.
  #
  # Inputs:
  #   Y - Matrix of count data for all items (n x J).
  #   N - Vector of item lengths (length J).
  #   J - Integer representing the number of items.
  #
  # Output:
  #   parms - List containing estimated parameters 'a', 'b', and 'nu'.
  # -------------------------------
  reading_data_parms_MOM <- function(Y, N, J) {
    
    # Estimate success probabilities for each item
    success.hat <- apply(Y, 2, mean, na.rm = TRUE) / N  # Vector of success rates (length J)
    
    eta.hat <- rep(NA, J)  # Initialize vector to store eta estimates for each item
    
    # Compute eta.hat for each item based on success.hat and N
    for (j in 1:J) {
      if (success.hat[j] < 1) {
        eta.hat[j] <- qnorm(success.hat[j])  # Standard quantile if success rate < 1
      } else {
        # Adjusted quantile if success rate >= 1 to prevent infinity
        eta.hat[j] <- qnorm(success.hat[j] * N[j]^2 / (1 + N[j]^2))
      }
    }
    
    # Initialize matrices to store indices and pairwise statistics
    index <- NULL
    mu.pair <- NULL
    n.pair <- NULL
    
    # Iterate over all unique pairs of items to compute pairwise means and sample sizes
    for (i in 1:(J - 1)) {
      for (j in (i + 1):J) {
        n.ij <- sum((1 - is.na(Y[, i])) * (1 - is.na(Y[, j])))  # Sample size for pair (i, j)
        if (n.ij > 0) {  # Only consider pairs with non-zero sample size
          index <- rbind(index, c(i, j))  # Store item indices
          mu.pair <- rbind(mu.pair, mean(Y[, i] * Y[, j], na.rm = TRUE) / (N[i] * N[j]))  # Pairwise mean
          n.pair <- rbind(n.pair, n.ij)  # Store sample size
        }
      }
    }
    
    K <- nrow(n.pair)  # Number of valid pairs
    r.store <- array(NA, dim = c(K, 1))  # Initialize array to store correlation estimates
    X <- array(0, dim = c(K, J))  # Initialize matrix to store design matrix
    
    # Estimate pairwise correlations using optimization
    for (k in 1:nrow(n.pair)) {
      mu.biv <- mu.pair[k]  # Observed pairwise mean
      eta.marg <- c(eta.hat[index[k, 1]], eta.hat[index[k, 2]])  # Marginal quantiles
      x <- seq(-3, 3, length.out = 20)  # Sequence of quantile values for optimization
      
      # Compute MD for each quantile in x
      f <- bivnorm.min_vectorized(x, mu.biv, eta.marg)
      x.start <- x[which(f == min(f))[1]]  # Select quantile with minimum MD
      
      # Optimize to find the quantile that minimizes MD
      optim.out <- optim(
        x.start, 
        fn = function(x) bivnorm.min_vectorized(x, mu.biv, eta.marg),
        method = "BFGS"
      )
      
      r.out <- pnorm(optim.out$par)  # Convert optimized quantile to probability
      r.store[k] <- r.out  # Store estimated correlation
      
      # Update design matrix X with indicators for the current pair
      X[k, index[k, 1]] <- 1
      X[k, index[k, 2]] <- 1
    }
    
    # Compute weighting factors based on optimization results
    rat.mom <- exp(
      lm(
        log(r.store) ~ 0 + X,
        weights = sqrt(r.store * (1 - r.store))
      )$coefficients
    )
    
    # Adjust weights that exceed 1 to prevent overestimation
    index <- which(rat.mom > 1)
    if (length(index) > 0) {
      rat.mom[index] <- 1 - 0.5 * (1 - max(rat.mom[-index]))
    }
    
    a.hat <- as.numeric(rat.mom / sqrt(1 - rat.mom^2))  # Estimate parameter 'a' for each item
    
    # Estimate 'mu2.hat' as the mean of squared counts for each item
    mu2.hat <- apply(Y^2, 2, mean, na.rm = TRUE)
    
    Phi2.hat <- rep(0, J)  # Initialize vector to store Phi2 estimates
    
    # Compute Phi2.hat for each item using the bivariate normal CDF
    for (j in 1:J) {
      r <- a.hat[j]^2 / (1 + a.hat[j]^2)  # Compute correlation 'r' based on 'a.hat'
      Phi2.hat[j] <- pmvnorm(
        upper = rep(eta.hat[j], 2),
        mean = rep(0, 2),
        corr = matrix(c(1, r, r, 1), nrow = 2),
        keepAttr = FALSE
      )
    }
    
    # Estimate 'nu.hat' based on observed means and Phi2 estimates
    nu.hat <- (mu2.hat - N^2 * Phi2.hat) / (N * (pnorm(eta.hat) - Phi2.hat))
    nu.hat <- pmax(1, pmin(nu.hat, N))  # Constrain 'nu.hat' within [1, N]
    
    # Estimate parameter 'b' based on 'eta.hat' and 'a.hat'
    b.hat <- eta.hat * sqrt(1 + a.hat^2)
    
    # Compile estimated parameters into a list
    parms <- list(a = a.hat, b = b.hat, nu = nu.hat)
    return(parms)
  }
  
  # -------------------------------
  # Method of Moments (MOM) Estimation
  # -------------------------------
  
  # Estimate parameters 'a', 'b', and 'nu' using MOM
  count.parms <- reading_data_parms_MOM(Y, N, J)
  a.out <- count.parms$a
  b.out <- count.parms$b
  nu.out <- count.parms$nu
  
  n <- dim(Y)[1]  # Number of observations
  nancov_T <- cov(logT10, use = "pairwise.complete.obs")  # Covariance matrix of log-time data
  
  beta.out <- apply(logT10, 2, mean, na.rm = TRUE)  # Estimate 'beta' as the mean of log-time data
  
  # Initialize variables for covariance adjustments
  count <- 1
  alpha.inv <- numeric(0)
  s22.comp <- numeric(0)
  
  # Compute inverse precision components and covariance adjustments
  for (j in 1:J) {
    if (j == J) {
      break  # Exit loop if j reaches the number of items
    }
    for (k in (j + 1):J) {
      alpha.inv <- rbind(
        alpha.inv,
        matrix(rep(0, J + 1), nrow = 1, ncol = (J + 1))
      )
      alpha.inv[count, j] <- max(0, nancov_T[j, j] - nancov_T[j, k])
      alpha.inv[count, k] <- max(0, nancov_T[k, k] - nancov_T[j, k])
      s22.comp <- rbind(s22.comp, nancov_T[j, k])
      count <- count + 1
    }
  }
  
  # Estimate 'alpha' based on inverse precision components
  alpha.inv_out <- apply(
    alpha.inv[, 1:J],
    2,
    sum,
    na.rm = TRUE
  ) / (J - apply(is.na(alpha.inv[, 1:J]), 2, sum) - 1)
  alpha.out <- 1 / sqrt(alpha.inv_out)  # Final 'alpha' estimates
  
  # Estimate 'vartau' as the mean of covariance components
  vartau.out <- max(0, mean(s22.comp, na.rm = TRUE))
  
  # Initialize matrix to store s12 calculations
  s12.calc <- matrix(rep(0, J), nrow = J)
  
  # Compute 's12' for each item based on covariance with log-time data
  for (j in 1:J) {
    CV <- cov(Y[, j], logT10[, j], use = "pairwise.complete.obs")  # Covariance between count and log-time
    s12.calc[j] <- -CV / N[j] * sqrt(2 * pi) * sqrt(a.out[j]^2 + 1) / a.out[j] *
      exp(0.5 * b.out[j]^2 / (1 + a.out[j]^2))
    s12.calc[j] <- max(-sqrt(vartau.out), min(sqrt(vartau.out), s12.calc[j]))  # Constrain 's12' within bounds
  }
  
  # Average 's12' across all items
  s12.out <- mean(s12.calc)
  
  # Compute intra-class correlation 'rho' based on 'vartau'
  if (vartau.out > 0) {
    s12.cor <- s12.out / sqrt(vartau.out)
  } else {
    s12.cor <- 0
  }
  
  # Compile all MOM estimates into a list
  MOMests <- list(
    a = as.numeric(a.out), 
    b = as.numeric(b.out),
    nu = as.numeric(nu.out),
    alpha = as.numeric(alpha.out), 
    beta = as.numeric(beta.out), 
    vartau = as.numeric(vartau.out), 
    rho = as.numeric(s12.cor)
  )
  
  # Return the MOM estimates
  return(MOMests)
}
