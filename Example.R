# ====================================================
# MOM and MCEM Code Demo
# ====================================================

# ====================================================
# Setup Instructions
# ====================================================

# 1. Download Required Files:
#    - Obtain the following files from the GitHub repository:
#        * mom_estimation.R      : Contains functions for Method of Moments estimations.
#        * mcem_binom.R          : Contains functions to perform MCEM for the Binomial model.
#        * mcem_betabin.R        : Contains functions to perform MCEM for the Beta-Binomial model.
#        * simulated_data_illustration.RData : Contains simulated data for demonstration.

# 2. Set Working Directory:
#    - Navigate to the folder where the downloaded files are located.
#      Example:
#      setwd("path/to/your/folder")

# 3. Install Required R Packages:
#    - Ensure that the following packages are installed for statistical computations:
#        * mvtnorm      : For multivariate normal distributions.
#        * extraDistr   : For additional statistical distributions.
#    - Install them using:
#      install.packages(c("mvtnorm", "extraDistr"))

# ====================================================
# Load Simulated Data
# ====================================================

# The RData file contains two lists:
#   1. parms:
#      - Description: True parameter values used to simulate the data.
#      - Example Usage: parms$a retrieves the accuracy discrimination parameters.

#   2. data:
#      - Description: Simulated dataset with the following elements:
#          * Count  : Number of tasks correct per item.
#          * logT10 : Standardized time per 10 tasks per item.
#          * N      : Number of tasks per item.
#          * n      : Sample size (number of observations).
#          * J      : Number of items/passages.

# Load the simulated data into the R environment
load("simulated_data_illustration.RData")

# ====================================================
# Load Required Functions
# ====================================================

# Source the external R scripts containing necessary functions for MOM and MCEM estimations
source("mom_estimation.R")  # Contains mom.binom and mom.betabin functions
source("mcem_binom.R")      # Contains run.mcem.binom function
source("mcem_betabin.R")    # Contains run.mcem.betabin function

# Load the required R packages for statistical computations
library(mvtnorm)      # For multivariate normal distributions
library(extraDistr)   # For additional statistical distributions

# ====================================================
# Illustrating the Binomial Model
# ====================================================

# ----------------------------
# Step 1: Method of Moments (MOM) Estimation for Binomial Model
# ----------------------------

# Calculate the binomial MOM estimators using the mom.binom function
# Inputs:
#   - data$Count : Matrix of counts of correct tasks (n x J)
#   - data$logT10: Matrix of standardized time per 10 tasks (n x J)
#   - data$N     : Vector of number of tasks per item (length J)
#   - data$J     : Number of items/passages
binom_ests_mom <- mom.binom(data$Count, data$logT10, data$N, data$J)

# ----------------------------
# Step 2: Monte Carlo Expectation-Maximization (MCEM) Estimation for Binomial Model
# ----------------------------

# Define the number of imputations and repetitions for MCEM
# k.in: Number of imputations per round (vector)
# reps.in: Number of rounds for each k.in value (vector)
k.in <- c(1, 100)   # First run: 1 imputation, Second run: 100 imputations
reps.in <- c(20, 3) # First run: 20 repetitions, Second run: 3 repetitions

# Perform MCEM estimation for the Binomial model using run.mcem.binom function
# Inputs:
#   - data$Count, data$logT10, data$N, data$J: Simulated data
#   - k.in, reps.in: Number of imputations and repetitions
#   - binom_ests_mom: Initial MOM estimates
binom_mcem_run <- run.mcem.binom(data$Count, data$logT10, data$N, data$J, 
                                 k.in = k.in, reps.in = reps.in, 
                                 binom_ests_mom)

# ----------------------------
# Step 3: Extract and Average MCEM Estimates for Binomial Model
# ----------------------------

# Initialize a list to store averaged MCEM estimates for the Binomial model
binom_ests_mcem <- list()

# Calculate the mean of parameter 'a' across specified MCEM iterations (rows 22 to 24)
binom_ests_mcem$a <- apply(binom_mcem_run$a[22:24, ], 2, mean)

# Calculate the mean of parameter 'b' similarly
binom_ests_mcem$b <- apply(binom_mcem_run$b[22:24, ], 2, mean)

# Initialize 'rho.intra' as a vector of zeros (assumed value implied by binomial model)
binom_ests_mcem$rho.intra <- rep(0, data$J)

# Calculate the mean of the inverse of 'alpha' parameter
binom_ests_mcem$alpha.inv <- apply(1 / binom_mcem_run$alpha[22:24, ], 2, mean)

# Calculate the mean of parameter 'beta'
binom_ests_mcem$beta <- apply(binom_mcem_run$beta[22:24, ], 2, mean)

# Calculate the mean of 'vartau' parameter across specified iterations
binom_ests_mcem$vartau <- mean(binom_mcem_run$vartau[22:24])  

# Calculate the mean of 'rho' (intra-class correlation) parameter
binom_ests_mcem$rho <- mean(binom_mcem_run$rho[22:24])        

# ----------------------------
# Step 4: Method of Moments (MOM) Estimation for Beta-Binomial Model
# ----------------------------

# Calculate the beta-binomial MOM estimators using the mom.betabin function
# Inputs are similar to the binomial MOM estimation
betabin_ests_mom <- mom.betabin(data$Count, data$logT10, data$N, data$J)

# ----------------------------
# Step 5: Monte Carlo Expectation-Maximization (MCEM) Estimation for Beta-Binomial Model
# ----------------------------

# Perform MCEM estimation for the Beta-Binomial model using run.mcem.betabin function
# Inputs:
#   - data$Count, data$logT10, data$N, data$J: Simulated data
#   - k.in, reps.in: Number of imputations and repetitions
#   - betabin_ests_mom: Initial MOM estimates for Beta-Binomial model
betabin_mcem_run <- run.mcem.betabin(data$Count, data$logT10, data$N, data$J, 
                                     k.in = k.in, reps.in = reps.in, 
                                     betabin_ests_mom)

# ----------------------------
# Step 6: Extract and Average MCEM Estimates for Beta-Binomial Model
# ----------------------------

# Initialize a list to store averaged MCEM estimates for the Beta-Binomial model
betabin_ests_mcem <- list()

# Calculate the mean of parameter 'a' across specified MCEM iterations (rows 22 to 24)
betabin_ests_mcem$a <- apply(betabin_mcem_run$a[22:24, ], 2, mean)

# Calculate the mean of parameter 'b' similarly
betabin_ests_mcem$b <- apply(betabin_mcem_run$b[22:24, ], 2, mean)

# Calculate the intra-class correlation 'rho.intra' based on the 'nu' parameter
betabin_ests_mcem$rho.intra <- (apply(betabin_mcem_run$nu[22:24, ], 2, mean) - 1) / (data$N - 1)

# Calculate the mean of the inverse of 'alpha' parameter
betabin_ests_mcem$alpha.inv <- apply(1 / betabin_mcem_run$alpha[22:24, ], 2, mean)

# Calculate the mean of parameter 'beta'
betabin_ests_mcem$beta <- apply(betabin_mcem_run$beta[22:24, ], 2, mean)

# Calculate the mean of 'vartau' parameter across specified iterations
betabin_ests_mcem$vartau <- mean(betabin_mcem_run$vartau[22:24])  

# Calculate the mean of 'rho' (intra-class correlation) parameter
betabin_ests_mcem$rho <- mean(betabin_mcem_run$rho[22:24])        

# ====================================================
# Calculate Root Mean Square Deviation (RMSD) for All Parameters
# ====================================================

# Initialize a list to store RMSD values for each parameter
root_mean_sq_dev <- list()

# ----------------------------
# RMSD for Parameter 'a'
# ----------------------------

# Calculate RMSD between true 'a' and Binomial MCEM estimates
# and between true 'a' and Beta-Binomial MCEM estimates
root_mean_sq_dev$a <- c(
  sqrt(mean((parms$a - binom_ests_mcem$a)^2)),  # RMSD for Binomial model
  sqrt(mean((parms$a - betabin_ests_mcem$a)^2)) # RMSD for Beta-Binomial model
)

# ----------------------------
# RMSD for Parameter 'b'
# ----------------------------

root_mean_sq_dev$b <- c(
  sqrt(mean((parms$b - binom_ests_mcem$b)^2)),  # RMSD for Binomial model
  sqrt(mean((parms$b - betabin_ests_mcem$b)^2)) # RMSD for Beta-Binomial model
)

# ----------------------------
# RMSD for Parameter 'rho.intra' (Dispersion -- Transformed 'nu')
# ----------------------------

# Calculate the true 'rho.intra' based on true 'nu' parameters
parms$rho.intra <- (parms$nu - 1) / (data$N - 1)

root_mean_sq_dev$rho.intra <- c(
  sqrt(mean((parms$rho.intra - binom_ests_mcem$rho.intra)^2)),  # RMSD for Binomial model
  sqrt(mean((parms$rho.intra - betabin_ests_mcem$rho.intra)^2)) # RMSD for Beta-Binomial model
)

# ----------------------------
# RMSD for Parameter 'alpha.inv'
# ----------------------------

root_mean_sq_dev$alpha.inv <- c(
  sqrt(mean((parms$alpha.inv - binom_ests_mcem$alpha.inv)^2)),  # RMSD for Binomial model
  sqrt(mean((parms$alpha.inv - betabin_ests_mcem$alpha.inv)^2)) # RMSD for Beta-Binomial model
)

# ----------------------------
# RMSD for Parameter 'beta'
# ----------------------------

root_mean_sq_dev$beta <- c(
  sqrt(mean((parms$beta - binom_ests_mcem$beta)^2)),  # RMSD for Binomial model
  sqrt(mean((parms$beta - betabin_ests_mcem$beta)^2)) # RMSD for Beta-Binomial model
)

# ----------------------------
# RMSD for Parameter 'vartau'
# ----------------------------

root_mean_sq_dev$vartau <- c(
  sqrt(mean((parms$vartau - binom_ests_mcem$vartau)^2)),  # RMSD for Binomial model
  sqrt(mean((parms$vartau - betabin_ests_mcem$vartau)^2)) # RMSD for Beta-Binomial model
)

# ----------------------------
# RMSD for Parameter 'rho'
# ----------------------------

root_mean_sq_dev$rho <- c(
  sqrt(mean((parms$rho - binom_ests_mcem$rho)^2)),  # RMSD for Binomial model
  sqrt(mean((parms$rho - betabin_ests_mcem$rho)^2)) # RMSD for Beta-Binomial model
)

# ================================
# Final Output: RMSD Summary
# ================================

# Display the RMSD for all parameters to assess estimation accuracy
# This provides a measure of how closely the estimated parameters match the true parameters
# Note -- this is based on a single simulated sample so is not any type of "proof" of method superiority
root_mean_sq_dev
