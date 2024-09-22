# MCEM and MOM Estimation for Binomial and Beta-Binomial Models

This repository contains R scripts implementing the Method of Moments (MOM) and Monte Carlo Expectation-Maximization (MCEM) algorithms for parameter estimation in both the binomial and beta-binomial models. The codes are designed to be easy to use for those who wish to replicate the estimation procedures used in our manuscript.

## Files Overview

1. **`mom_estimation.R`**:  
   This script contains the functions for Method of Moments (MOM) estimators, providing initial starting values for / an alternative to the MCEM algorithm.
   
2. **`mcem_binom.R`**:  
   This script provides the MCEM algorithm for estimating parameters in the binomial count-time model. The script iteratively imputes latent variables and updates model parameters to approximate the maximum likelihood estimates for the binomial model.

3. **`mcem_betabin.R`**:  
   This script provides the MCEM algorithm for estimating parameters in the beta-binomial count-time model. It includes a similar structure to `mcem_binom.R`, but is tailored for the beta-binomial distribution, allowing for over-dispersion in the data.

4. **`simulation_and_modeling.R`**:  
   This script is only relevant for those who wish to recreate the simulations described in our manuscript. It provides the complete code necessary for simulating data under various conditions and fitting the models to the simulated data.

5. **`example.R`**:  
   This script provides a guided example of how to use the MCEM and MOM estimation functions for both the binomial and beta-binomial models. It walks through the process of loading data, running the algorithms, and interpreting the output.
