# -----------------------------------------------------------------------------
# Filename: simulation_and_modeling.R
# Description: 
#   This R script contains functions to simulate data and model it using both
#   Binomial and Beta-Binomial count models within a bivariate latent factor 
#   framework. The code can replicate the simulation as described in the paper 
#   "Joint Analysis of Dispersed Count-Time Data 
#   Using a Bivariate Latent Factor Model".
#
# Contents of this file:
#   1. simulate.data(seed.num)
#      - Simulates count and log-time data based on specified parameters.
#      - "seed.num" describes different dispersion scenarios (No Dispersion, 
#        Low & High Dispersion) and varying sample sizes (n = 600, 1200, 2000).
#      - Outputs simulated parameters and data along with filenames for storage.
#
#   2. model.data.binom(data)
#      - Fits the Binomial count model to the simulated data.
#      - Utilizes Method of Moments (MOM) for initial parameter estimates.
#      - Implements Monte Carlo Expectation-Maximization (MCEM) for parameter 
#        estimation with two iteration settings.
#      - Records and returns execution times for each modeling step.
#
#   3. model.data.betabin(data)
#      - Fits the Beta-Binomial count model to the simulated data.
#      - Similar structure to model.data.binom but accounts for potential 
#        overdispersion in count data.
#      - Utilizes MOM for initial estimates and MCEM for parameter estimation.
#      - Records and returns execution times for each modeling step.
#
# Dependencies:
#   - mvtnorm: For multivariate normal distributions.
#   - extraDistr: For the Beta-Binomial distribution.
# -----------------------------------------------------------------------------

simulate.data <- function(seed.num) {
  
  # ---------------------------------------------------------------------------
  # Function: simulate.data
  # Purpose:
  #   Simulates a dataset based on the joint count-time latent factor model.
  #   Generates count data (Y_ij) using a Beta-Binomial or Binomial distribution
  #   and log-scale time data (T_ij) using a Normal distribution.
  #
  # Inputs:
  #   seed.num - Integer seed number for reproducibility and to determine
  #              simulation configurations. The paper used seed.num = 1,...,6000
  #              for 1000 Monte Carlo runs of 6 different model configurations.
  #
  # Outputs:
  #   A list containing:
  #     - parms.sim: List of simulated model parameters.
  #     - data.sim: List of simulated data (Count and logT10 matrices).
  #     - filename: Generated filename for saving the simulation.
  #     - full_path: Full path to the file for saving the simulation.
  # ---------------------------------------------------------------------------
  
  # Set a seed for simulating common model parameters -- fixed across 
  # each one of the six different model configurations in the paper.
  set.seed(22724)
  J <- 50  # Number of items (passages)
  
  # Ensure seed.num is an integer
  seed.num <- as.integer(seed.num)
  
  # Load necessary libraries
  library(mvtnorm)    # For multivariate normal distributions
  library(extraDistr) # For Beta-Binomial distribution
  
  # Determine sample size (n) and dispersion scenario based on seed.num
  # There are 6 different scenarios here:
  #   - n = 600, disp.num = 1 (No Dispersion)
  #   - n = 600, disp.num = 2 (Low Dispersion)
  #   - n = 600, disp.num = 3 (High Dispersion)
  #   - n = 1200, disp.num = 1 (No Dispersion)
  #   - n = 1200, disp.num = 2 (Low Dispersion)
  #   - n = 1200, disp.num = 3 (High Dispersion)
  #   - n = 2000, disp.num = 1 (No Dispersion)
  #   - n = 2000, disp.num = 2 (Low Dispersion)
  #   - n = 2000, disp.num = 3 (High Dispersion)
  if (seed.num <= 3000) {
    n <- 600
    disp.num <- as.integer(ceiling(seed.num / 1000))
    if (disp.num == 1) {
      rho.intra <- rep(0, J)  # No Dispersion
    } else if (disp.num == 2) {
      rho.intra <- runif(J, min = 0.03, max = 0.1)  # Low Dispersion
    } else {
      rho.intra <- runif(J, min = 0.1, max = 0.17)  # High Dispersion
    }
  } else if (seed.num <= 6000) {
    n <- 1200
    disp.num <- as.integer(ceiling((seed.num - 3000) / 1000))
    if (disp.num == 1) {
      rho.intra <- rep(0, J)
    } else if (disp.num == 2) {
      rho.intra <- runif(J, min = 0.03, max = 0.1)
    } else {
      rho.intra <- runif(J, min = 0.1, max = 0.17)
    }
  } else {
    n <- 2000
    disp.num <- as.integer(ceiling((seed.num - 6000) / 1000))
    if (disp.num == 1) {
      rho.intra <- rep(0, J)
    } else if (disp.num == 2) {
      rho.intra <- runif(J, min = 0.03, max = 0.1)
    } else {
      rho.intra <- runif(J, min = 0.1, max = 0.17)
    }
  }
  
  # Generate filenames for saving simulation results
  filename <- paste("sim_n", n, "_disp", 
                    disp.num, "_seed", seed.num,
                    ".Rdata", sep = "")
  
  dirname <- paste("Sim_n", n, sep = "")
  
  full_path <- file.path(dirname, filename)
  
  # Simulate item lengths (N_j) for 50 items
  N <- ceiling(c(runif(30, min = 49, max = 59),  # 30 medium-length items
                 runif(20, min = 74, max = 99)))  # 20 longer items
  
  # Assign items to test-takers based on a balanced incomplete block design
  # Ensures that each test-taker reads a subset of items, introducing missingness
  items.assigned <- c(1:10)
  for (j in 1:9) {
    items.assigned <- rbind(items.assigned, c(1:10) + 5 * j)
  }
  items.assigned[10, 6:10] <- items.assigned[10, 6:10] - 50  # Adjust for overlap in the last block
  
  # Define the number of passages read and their probabilities
  # This generates the Missing-at-Random (MAR) structure
  num.read <- 1:10
  p.num.read <- c(0.17, 0.43, 0.03, 0.03, 0.04, 0.05, 0.09, 0.14, 0.01, 0.01)
  
  # Simulate count model parameters
  a <- runif(J, min = 0.3, max = 0.8)  # Discrimination parameters
  p.mean.success <- runif(J, min = 0.8, max = 0.9)  # Mean success probabilities
  b <- qnorm(p.mean.success) * sqrt(1 + a^2)  # Difficulty parameters based on probit link
  
  # Calculate dispersion parameters based on simulated intra-class correlation
  nu <- rho.intra * (N - 1) + 1
  alpha.inv <- runif(J, min = 0.07, max = 0.28)  # Inverse precision parameters
  alpha <- 1 / alpha.inv  # Precision parameters
  beta <- runif(J, min = 1.7, max = 2.1)  # Speed model intercepts
  
  # Latent variable parameters
  vartau <- 0.15  # Variance of latent speed factor
  rho <- 0.4       # Correlation between latent accuracy and speed
  
  # Initialize matrices to store simulated data
  Count <- array(NA, dim = c(n, J))      # Count data matrix
  logT10 <- array(NA, dim = c(n, J))    # Log-time data matrix
  
  # Define latent variable distribution parameters
  mu <- c(0, 0)  # Means of latent traits
  Sigma <- matrix(c(1, rho * sqrt(vartau),
                    rho * sqrt(vartau), vartau), nrow = 2)  # Covariance matrix
  
  # Set unique seed for data generation so that each dataset is unique
  # (earlier seed 22724 ensured common model parameters)
  set.seed(seed.num)
  
  # Simulate latent variables for all test-takers
  theta <- rmvnorm(n, mean = mu, sigma = Sigma)
  tau <- theta[, 2]  # Latent speed factor
  theta <- theta[, 1]  # Latent accuracy factor
  
  # Generate data for each test-taker
  for (i in 1:n) {
    # Select items assigned to the test-taker based on block design
    items.select <- items.assigned[10 - (i %% 10), ]
    
    # Determine number of passages read by the test-taker
    num.pass <- sample(num.read, 1, prob = p.num.read)
    
    # Randomly select passages to read without replacement
    pass.select <- sample(items.select, num.pass, replace = FALSE)
    
    # Calculate success probabilities using the probit link function
    p.success <- pnorm(a * theta[i] + b)
    
    # Parameters for the Beta-Binomial distribution
    alpha.bbin <- (N - 1) / (nu - 1) * p.success
    beta.bbin <- (N - 1) / (nu - 1) * (1 - p.success)
    
    # Simulate count and log-time data for each selected passage
    for (j in pass.select) {
      if (nu[j] < 1.0001) {
        # Use Binomial distribution if dispersion parameter is approximately 1
        Count[i, j] <- rbinom(1, N[j], p.success[j])
      } else {
        # Use Beta-Binomial distribution to account for overdispersion
        Count[i, j] <- rbbinom(1, N[j], alpha.bbin[j], beta.bbin[j])
      }
      # Simulate log-time data based on the Normal distribution
      logT10[i, j] <- rnorm(1, beta[j] - tau[i], alpha.inv[j])
    }
  }
  
  # Compile simulated parameters and data into lists
  parms.sim <- list(a = a, b = b, 
                    rho.intra = rho.intra, nu = nu, 
                    alpha.inv = alpha.inv, alpha = alpha, 
                    beta = beta, vartau = vartau, rho = rho)
  
  data.sim <- list(Count = Count, 
                   logT10 = logT10, 
                   N = N, n = n, J = J)
  
  # Return the simulated parameters, data, and file information
  return(list(parms.sim = parms.sim, 
              data.sim = data.sim, 
              filename = filename,
              full_path = full_path))
  
}

model.data.binom <- function(data) {
  
  # ---------------------------------------------------------------------------
  # Function: model.data.binom
  # Purpose:
  #   Fits the Binomial count model to the simulated data.
  #   Applies Method of Moments (MOM) for initial parameter estimation and
  #   then uses Monte Carlo Expectation-Maximization (MCEM) for parameter
  #   refinement with two iteration settings.
  #
  # Inputs:
  #   data - List containing simulated data:
  #          - Count: Matrix of count data (n x J)
  #          - logT10: Matrix of log-time data (n x J)
  #          - N: Vector of item lengths (J)
  #          - n: Sample size (integer)
  #          - J: Number of items (integer)
  #
  # Outputs:
  #   A list containing:
  #     - mom.bin: MOM estimates.
  #     - mcem.bin1: MCEM estimates (50 runs with K = 1, 2 runs with K = 100).
  #     - mcem.bin2: MCEM estimates (200 runs with K = 1, 5 runs with K = 100).
  # ---------------------------------------------------------------------------
  
  # Extract data components
  Count <- data$Count
  logT10 <- data$logT10
  N <- data$N
  n <- data$n
  J <- data$J
  
  # -------------------------------
  # Method of Moments (MOM) for Binomial Model
  # -------------------------------
  
  # Start timing for MOM estimation
  start_time_binom <- Sys.time()
  
  # Apply MOM to obtain initial parameter estimates for Binomial model
  mom.bin <- mom.binom(Count, logT10, N, J)
  
  # End timing for MOM estimation
  end_time_binom <- Sys.time()
  
  # Calculate execution time
  execution_time_mombinom <- end_time_binom - start_time_binom
  cat("mom.binom execution time:", execution_time_mombinom, "\n")
  
  # Calculate inverse precision and record execution time
  mom.bin$alpha.inv <- 1 / mom.bin$alpha
  mom.bin$time <- as.numeric(execution_time_mombinom, units = "secs")
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Binomial Model - First Run
  # -------------------------------
  
  # Start timing for first MCEM run
  start_time_bin <- Sys.time()
  
  # Run MCEM with k.in = 1 imputation for each of the reps.in = 50 iterations
  mcem.bin.k1.r50 <- run.mcem.binom(Count, logT10, N, J, k.in = 1, reps.in = 50, mom.bin)
  
  # End timing for first MCEM run
  end_time_bin <- Sys.time()
  
  # Calculate execution time
  execution_time_bin_k1_r50 <- end_time_bin - start_time_bin
  cat("mcem.bin execution time:", execution_time_bin_k1_r50, "\n")
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Binomial Model - First Run Continued
  # -------------------------------
  
  # Start timing for finishing first MCEM run
  start_time_bin <- Sys.time()
  
  # Run MCEM with k.in = 100 imputation for each of the reps.in = 2 iterations
  mcem.bin.k100.r2 <- run.mcem.binom(Count, logT10, N, J, k.in = 100, reps.in = 2, 
                                     list(a = mcem.bin.k1.r50$a[51,],
                                          b = mcem.bin.k1.r50$b[51,],
                                          alpha = mcem.bin.k1.r50$alpha[51,],
                                          beta = mcem.bin.k1.r50$beta[51,],
                                          vartau = mcem.bin.k1.r50$vartau[51],
                                          rho = mcem.bin.k1.r50$rho[51]))
  
  # End timing for second MCEM run
  end_time_bin <- Sys.time()
  
  # Calculate execution time
  execution_time_bin_k100_r2 <- end_time_bin - start_time_bin
  cat("mcem.bin execution time:", execution_time_bin_k100_r2, "\n")
  
  # Average parameter estimates from the last two iterations to obtain mcem.bin1
  mcem.bin1 <- list(a = colMeans(mcem.bin.k100.r2$a[2:3,]),
                    b = colMeans(mcem.bin.k100.r2$b[2:3,]),
                    alpha = colMeans(mcem.bin.k100.r2$alpha[2:3,]),
                    beta = colMeans(mcem.bin.k100.r2$beta[2:3,]),
                    alpha.inv = colMeans(1/mcem.bin.k100.r2$alpha[2:3,]),
                    vartau = mean(mcem.bin.k100.r2$vartau[2:3]),
                    rho = mean(mcem.bin.k100.r2$rho[2:3]),
                    time = as.numeric(execution_time_bin_k1_r50 + 
                                        execution_time_bin_k100_r2, units = "secs"))
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Binomial Model - Second Run (Continued from First)
  # -------------------------------
  
  # Start timing for continued MCEM run 
  start_time_bin <- Sys.time()
  
  # Run MCEM with k.in = 1 imputation for a further reps.in = 150 iterations
  mcem.bin.k1.r150 <- run.mcem.binom(Count, logT10, N, J, k.in = 1, reps.in = 150, 
                                     list(a = mcem.bin.k1.r50$a[51,],
                                          b = mcem.bin.k1.r50$b[51,],
                                          alpha = mcem.bin.k1.r50$alpha[51,],
                                          beta = mcem.bin.k1.r50$beta[51,],
                                          vartau = mcem.bin.k1.r50$vartau[51],
                                          rho = mcem.bin.k1.r50$rho[51]))
  
  # End timing for third MCEM run
  end_time_bin <- Sys.time()
  
  # Calculate execution time
  execution_time_bin_k1_r150 <- end_time_bin - start_time_bin
  cat("mcem.bin execution time:", execution_time_bin_k1_r150, "\n")
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Binomial Model - Second Run (Part 2)
  # -------------------------------
  
  # Start timing for MCEM run (finish the Second Run)
  start_time_bin <- Sys.time()
  
  # Run MCEM with k.in = 100 imputation for each of the reps.in = 5 iterations
  mcem.bin.k100.r5 <- run.mcem.binom(Count, logT10, N, J, k.in = 100, reps.in = 5, 
                                     list(a = mcem.bin.k1.r150$a[151,],
                                          b = mcem.bin.k1.r150$b[151,],
                                          alpha = mcem.bin.k1.r150$alpha[151,],
                                          beta = mcem.bin.k1.r150$beta[151,],
                                          vartau = mcem.bin.k1.r150$vartau[151],
                                          rho = mcem.bin.k1.r150$rho[151]))
  
  # End timing for fourth MCEM run
  end_time_bin <- Sys.time()
  
  # Calculate execution time
  execution_time_bin_k100_r5 <- end_time_bin - start_time_bin
  cat("mcem.bin execution time:", execution_time_bin_k100_r5, "\n")
  
  # Average parameter estimates from the last five iterations to obtain mcem.bin2
  mcem.bin2 <- list(a = colMeans(mcem.bin.k100.r5$a[2:6,]),
                    b = colMeans(mcem.bin.k100.r5$b[2:6,]),
                    alpha = colMeans(mcem.bin.k100.r5$alpha[2:6,]),
                    beta = colMeans(mcem.bin.k100.r5$beta[2:6,]),
                    alpha.inv = colMeans(1/mcem.bin.k100.r5$alpha[2:6,]),
                    vartau = mean(mcem.bin.k100.r5$vartau[2:6]),
                    rho = mean(mcem.bin.k100.r5$rho[2:6]),
                    time = as.numeric(execution_time_bin_k1_r50 + 
                                        execution_time_bin_k100_r2 +
                                        execution_time_bin_k100_r5, units = "secs"))
  
  # Return the MOM and MCEM estimates for the Binomial model
  return(list(mom.bin = mom.bin,
              mcem.bin1 = mcem.bin1,
              mcem.bin2 = mcem.bin2))
  
}

model.data.betabin <- function(data) {
  
  # ---------------------------------------------------------------------------
  # Function: model.data.betabin
  # Purpose:
  #   Fits the Beta-Binomial count model to the simulated data.
  #   Applies Method of Moments (MOM) for initial parameter estimation and
  #   then uses Monte Carlo Expectation-Maximization (MCEM) for parameter
  #   refinement with two iteration settings.
  #
  # Inputs:
  #   data - List containing simulated data:
  #          - Count: Matrix of count data (n x J)
  #          - logT10: Matrix of log-time data (n x J)
  #          - N: Vector of item lengths (J)
  #          - n: Sample size (integer)
  #          - J: Number of items (integer)
  #
  # Outputs:
  #   A list containing:
  #     - mom.betabin: MOM estimates for the Beta-Binomial model.
  #     - mcem.betabin1: MCEM estimates (50 runs with K = 1, 2 runs with K = 100).
  #     - mcem.betabin2: MCEM estimates (200 runs with K = 1, 5 runs with K = 100).
  # ---------------------------------------------------------------------------
  
  # Extract data components
  Count <- data$Count
  logT10 <- data$logT10
  N <- data$N
  n <- data$n
  J <- data$J
  
  # -------------------------------
  # Method of Moments (MOM) for Beta-Binomial Model
  # -------------------------------
  
  # Start timing for MOM estimation
  start_time_betabin <- Sys.time()
  
  # Apply MOM to obtain initial parameter estimates for Beta-Binomial model
  mom.betabin <- mom.betabin(Count, logT10, N, J)
  
  # End timing for MOM estimation
  end_time_betabin <- Sys.time()
  
  # Calculate execution time
  execution_time_mombetabin <- end_time_betabin - start_time_betabin
  cat("mom.betabin execution time:", execution_time_mombetabin, "\n")
  
  # Calculate inverse precision, intra-class correlation, and record execution time
  mom.betabin$alpha.inv <- 1 / mom.betabin$alpha
  mom.betabin$rho.intra <- (mom.betabin$nu - 1) / (N - 1)
  mom.betabin$time <- as.numeric(execution_time_mombetabin, units = "secs")
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Beta-Binomial Model - First Run
  # -------------------------------
  
  # Start timing for first MCEM run
  start_time_betabin <- Sys.time()
  
  # Run MCEM with k.in = 1 imputation for each of the reps.in = 50 iterations
  mcem.betabin.k1.r50 <- run.mcem.betabin(
    Count, logT10, N, J,
    k.in = 1,
    reps.in = 50,
    mom.betabin
  )
  
  # End timing for first MCEM run
  end_time_betabin <- Sys.time()
  
  # Calculate execution time
  execution_time_betabin_k1_r50 <- end_time_betabin - start_time_betabin
  cat("mcem.betabin execution time:", execution_time_betabin_k1_r50, "\n")
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Beta-Binomial Model - Second Run
  # -------------------------------
  
  # Start timing for second MCEM run
  start_time_betabin <- Sys.time()
  
  # Run MCEM with k.in = 100 imputations for each of the reps.in = 2 iterations
  mcem.betabin.k100.r2 <- run.mcem.betabin(
    Count, logT10, N, J,
    k.in = 100,
    reps.in = 2,
    list(
      a = mcem.betabin.k1.r50$a[51, ],
      b = mcem.betabin.k1.r50$b[51, ],
      nu = mcem.betabin.k1.r50$nu[51, ],
      alpha = mcem.betabin.k1.r50$alpha[51, ],
      beta = mcem.betabin.k1.r50$beta[51, ],
      vartau = mcem.betabin.k1.r50$vartau[51],
      rho = mcem.betabin.k1.r50$rho[51]
    )
  )
  
  # End timing for second MCEM run
  end_time_betabin <- Sys.time()
  
  # Calculate execution time
  execution_time_betabin_k100_r2 <- end_time_betabin - start_time_betabin
  cat("mcem.betabin execution time:", execution_time_betabin_k100_r2, "\n")
  
  # Average parameter estimates from the last two iterations to obtain mcem.betabin1
  mcem.betabin1 <- list(
    a = colMeans(mcem.betabin.k100.r2$a[2:3, ]),
    b = colMeans(mcem.betabin.k100.r2$b[2:3, ]),
    nu = colMeans(mcem.betabin.k100.r2$nu[2:3, ]),
    rho.intra = colMeans(
      (mcem.betabin.k100.r2$nu[2:3, ] - 1) / 
        (t(array(rep(N, 2), dim = c(J, 2))) - 1)
    ),
    alpha = colMeans(mcem.betabin.k100.r2$alpha[2:3, ]),
    beta = colMeans(mcem.betabin.k100.r2$beta[2:3, ]),
    alpha.inv = colMeans(1 / mcem.betabin.k100.r2$alpha[2:3, ]),
    vartau = mean(mcem.betabin.k100.r2$vartau[2:3]),
    rho = mean(mcem.betabin.k100.r2$rho[2:3]),
    time = as.numeric(
      execution_time_betabin_k1_r50 + 
        execution_time_betabin_k100_r2,
      units = "secs"
    )
  )
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Beta-Binomial Model - Third Run
  # -------------------------------
  
  # Start timing for third MCEM run
  start_time_betabin <- Sys.time()
  
  # Run MCEM with k.in = 1 imputations for each of the reps.in = 150 iterations
  mcem.betabin.k1.r150 <- run.mcem.betabin(
    Count, logT10, N, J,
    k.in = 1,
    reps.in = 150,
    list(
      a = mcem.betabin.k1.r50$a[51, ],
      b = mcem.betabin.k1.r50$b[51, ],
      nu = mcem.betabin.k1.r50$nu[51, ],
      alpha = mcem.betabin.k1.r50$alpha[51, ],
      beta = mcem.betabin.k1.r50$beta[51, ],
      vartau = mcem.betabin.k1.r50$vartau[51],
      rho = mcem.betabin.k1.r50$rho[51]
    )
  )
  
  # End timing for third MCEM run
  end_time_betabin <- Sys.time()
  
  # Calculate execution time
  execution_time_betabin_k1_r150 <- end_time_betabin - start_time_betabin
  cat("mcem.betabin execution time:", execution_time_betabin_k1_r150, "\n")
  
  # -------------------------------
  # Monte Carlo EM (MCEM) for Beta-Binomial Model - Fourth Run
  # -------------------------------
  
  # Start timing for fourth MCEM run
  start_time_betabin <- Sys.time()
  
  # Run MCEM with k.in = 100 imputations for each of the reps.in = 5 iterations
  mcem.betabin.k100.r5 <- run.mcem.betabin(
    Count, logT10, N, J,
    k.in = 100,
    reps.in = 5,
    list(
      a = mcem.betabin.k1.r150$a[151, ],
      b = mcem.betabin.k1.r150$b[151, ],
      nu = mcem.betabin.k1.r150$nu[151, ],
      alpha = mcem.betabin.k1.r150$alpha[151, ],
      beta = mcem.betabin.k1.r150$beta[151, ],
      vartau = mcem.betabin.k1.r150$vartau[151],
      rho = mcem.betabin.k1.r150$rho[151]
    )
  )
  
  # End timing for fourth MCEM run
  end_time_betabin <- Sys.time()
  
  # Calculate execution time
  execution_time_betabin_k100_r5 <- end_time_betabin - start_time_betabin
  cat("mcem.betabin execution time:", execution_time_betabin_k100_r5, "\n")
  
  # Average parameter estimates from the last five iterations to obtain mcem.betabin2
  mcem.betabin2 <- list(
    a = colMeans(mcem.betabin.k100.r5$a[2:6, ]),
    b = colMeans(mcem.betabin.k100.r5$b[2:6, ]),
    nu = colMeans(mcem.betabin.k100.r5$nu[2:6, ]),
    rho.intra = colMeans(
      (mcem.betabin.k100.r5$nu[2:6, ] - 1) / 
        (t(array(rep(N, 2), dim = c(J, 5))) - 1)
    ),
    alpha = colMeans(mcem.betabin.k100.r5$alpha[2:6, ]),
    beta = colMeans(mcem.betabin.k100.r5$beta[2:6, ]),
    alpha.inv = colMeans(1 / mcem.betabin.k100.r5$alpha[2:6, ]),
    vartau = mean(mcem.betabin.k100.r5$vartau[2:6]),
    rho = mean(mcem.betabin.k100.r5$rho[2:6]),
    time = as.numeric(
      execution_time_betabin_k1_r50 + 
        execution_time_betabin_k1_r150 +
        execution_time_betabin_k100_r5,
      units = "secs"
    )
  )
  
  # Return the MOM and MCEM estimates for the Beta-Binomial model
  return(list(
    mom.betabin = mom.betabin,
    mcem.betabin1 = mcem.betabin1,
    mcem.betabin2 = mcem.betabin2
  ))
  
}