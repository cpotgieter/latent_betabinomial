# Main function implementing the MCEM algorithm
# Parameters: (a,b,alpha,beta,sigma_tau,rho)
# Y: Matrix of count data
# logT10: Matrix of log-scale speed data
# N: Vector with number of words in each passage
# J: Number of passages
# k.in, reps.in: Parameters controlling the number of imputations
# ests.in: Initial values for model parameters

run.mcem.binom <- function(Y, logT10, N, J, k.in=5, reps.in=2, ests.in) {
  
  # Function to compute the beta-binomial negative log-likelihood
  # Used in the rejection sampling algorithm for the MCEM sampler
  neg_logratio_binom <- function(data, par) {
    # Extracting data elements
    Y <- data[1,]
    N <- data[2,]
    a <- data[3,]
    b <- data[4,]
    z <- par
    
    # Calculating probabilities and distribution parameters
    success.prob <- pnorm(a * z + b)
    
    # Calculating negative log-likelihood for both cases
    negLP <- -sum(dbinom(Y, N, success.prob, log = TRUE))

    return(negLP)
  }
  
  # Function to calculate the gamma multiplier for the rejection sampler
  gamma_multiplier_binom <- function(Y, N, a, b) {
    # Check for special case with perfect score condition
    if (sum(N - Y) < 0.01) {
      gammax <- 1
      zmax <- Inf
    } else {
      # Optimizing the negative log-ratio
      zgam <- optimize(neg_logratio_binom, 
                       data = rbind(Y, N, a, b), 
                       interval = c(-7, 7))
      zmax <- zgam$minimum
      gammax <- exp(-zgam$objective)
    }
    
    # Returning the calculated values
    gamval <- list(zmax = zmax, gammax = gammax)
    return(gamval)
  }
  
  # Function for the imputation step in the MCEM algorithm
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
    
    # Loop over each individual for rejection sampler setup
    for (i in 1:n) {
      index <- which(!is.na(Y[i,]))
      gamval <- gamma_multiplier_binom(Y[i, index], N[index], 
                                         a[index], b[index])
      gammax[i] <- gamval$gammax
      zmax[i] <- gamval$zmax
      A[i] <- sum(alpha[index]^2)
      mu1[i] <- -rho * sigma_tau / (1 + A[i] * sigma_tau^2) * sum(alpha[index]^2 * (logT10[i, index] - beta[index]))
      sigma1_2[i] <- 1 / (1 + (A[i] * rho^2 * sigma_tau^2) / (1 + A[i] * (sigma_tau^2 - rho^2 * sigma_tau^2)))
      sigma2_2[i] <- (sigma_tau^2 - rho^2 * sigma_tau^2) / ((sigma_tau^2 - rho^2 * sigma_tau^2) * A[i] + 1)
    }
    
    # Imputation step setup
    Z1 <- matrix(rep(0, n * K), nrow = n, ncol = K)
    mu2 <- matrix(rep(0, n * K), nrow = n, ncol = K)
    tau <- matrix(rep(0, n * K), nrow = n, ncol = K)
    
    # Nested loop for actual imputation
    for (i in 1:n) {
      index <- which(!is.na(Y[i,]))
      for (k in 1:K) {
        check <- 0
        while (check < 1) {
          Z1[i, k] <- mu1[i] + sqrt(sigma1_2[i]) * rnorm(1)
          U <- runif(1)
          LP <- exp(-neg_logratio_binom(rbind(Y[i, index], N[index], 
                                                a[index], b[index]),
                                        Z1[i, k])) / gammax[i]
          if (U < LP) {
            check <- 2
          }
        }
        mu2[i, k] <- -((sigma_tau^2 - rho^2 * sigma_tau^2) * (sum(alpha[index]^2 * (logT10[i, index] - beta[index]))) - rho * sigma_tau * Z1[i, k]) / ((sigma_tau^2 - rho^2 * sigma_tau^2) * A[i] + 1)
        tau[i, k] <- mu2[i, k] + sqrt(sigma2_2[i]) * rnorm(1)
      }
    }
    
    # Returning the imputed values
    theta <- Z1
    imputes <- list(theta = theta, tau = tau, z.opt = zmax)
  }
  
  # MCEM function for binomial case
  mcem_count_bin <- function(par, Y, N, theta, K) {
    a <- par[1]
    b <- par[2]
    success.prob <- pnorm(a * theta + b)
    Y.rep <- t(array(rep(Y, K), dim = c(length(Y), K)))
    nllh <- -sum(apply(dbinom(Y.rep, N, success.prob, log = TRUE), 2, mean))
    return(nllh)
  }
  
  # Function to optimize sigma_tau and rho from the imputed latent variables
  function_s12s22_to_min <- function(data, par) {
    theta <- data$theta
    tau <- data$tau
    s12 <- par
    s22 <- s12^2 * (1 + mean(theta^2)) - 2 * s12 * mean(theta * tau) + mean(tau^2)
    F1 <- 0.5 * log(s22 - s12^2)
    F2 <- 0.5 / sqrt(s22 - s12^2) * (s22 * mean(theta^2) - 2 * s12 * mean(theta * tau) + mean(tau^2))
    Fval <- F1 + F2
    return(Fval)
  }
  
  # Function for a single iteration of the MCEM algorithm.
  MCEM_algorithm_one_iteration <- function(Y, logT10, N, J, ests, K) {
    
    # Imputation step 
    EMimps <- imputation_code_binom(Y, logT10, N, ests$a, ests$b, 
                                    ests$alpha, ests$beta, 
                                    sqrt(ests$vartau), ests$rho, K)
    
    # Initializing parameters for optimization
    EM.a <- rep(NA, J)
    EM.b <- rep(NA, J)
    EM.alpha <- rep(NA, J)
    EM.beta <- rep(NA, J)
    theta <- t(EMimps$theta)
    tau <- t(EMimps$tau)
    
    # Optimization loop for each passage
    for (j in 1:J) {
      index <- which(!is.na(Y[, j]))
      bin.sol <- optim(c(ests$a[j], ests$b[j]), 
                       mcem_count_bin, 
                       Y = Y[index, j], N = N[j], 
                       theta = theta[, index], 
                       K = K, method = "BFGS")
      
      EM.a[j] <- bin.sol$par[1]
      EM.b[j] <- bin.sol$par[2]
      
      if (K > 1) {
        index <- which(!is.na(logT10[, j]))
        EM.beta[j] <- mean(logT10[index,j] + apply(tau[,index],2,mean))
        logT10_vector <- t(matrix(logT10[index, j], nrow = 1))
        result <- sweep(tau[, index], 2, logT10_vector, "+") - EM.beta[j]
        alp2_inv <- mean(apply(result^2,2,mean))
        EM.alpha[j] <- 1 / sqrt(alp2_inv)
      } else {
        index <- which(!is.na(logT10[, j]))
        EM.beta[j] <- mean(logT10[index,j] + tau[,index])
        logT10_vector <- t(matrix(logT10[index, j], nrow = 1))
        result <- tau[, index] + logT10_vector - EM.beta[j]
        alp2_inv <- mean(result^2)
        EM.alpha[j] <- 1 / sqrt(alp2_inv)
      }
      
    }
    
    # Optimization for sigma_tau and rho
    s12.find <- optim(ests$rho * sqrt(ests$vartau), 
                      fn = function_s12s22_to_min, 
                      method = "BFGS", 
                      data = EMimps)
    s12.min <- s12.find$par
    EM.vartau <- s12.min^2 * (1 + mean(EMimps$theta^2)) - 2 * s12.min * mean(EMimps$theta * EMimps$tau) + mean(EMimps$tau^2)
    EM.rho <- s12.min / sqrt(EM.vartau)
    
    # Returning the estimated parameters
    EM.ests <- list(a = EM.a, b = EM.b, alpha = EM.alpha, 
                    beta = EM.beta, vartau = EM.vartau, rho = EM.rho)
    return(EM.ests)
  }
  
  # Handling the k.in and reps.in parameters for multiple imputations
  # i.e. how many imputations for how many reps
  nK <- length(k.in)
  if (missing(ests.in)) {
    # mom functions defined in different file
    ests.in <- mom.binom(Y,logT10,N,J) 
  }
  
  # Checking and adjusting infinite alpha values.
  alpha.check <- ests.in$alpha
  infIndex0 <- which(!is.infinite(alpha.check))
  infIndex1 <- which(is.infinite(alpha.check))
  alpha.check[infIndex1] <- max(alpha.check[infIndex0]) * 10
  ests.in$alpha <- alpha.check
  
  # Setting up initial imputation run
  n <- dim(Y)[1]
  total.K <- rep(k.in[1], reps.in[1])
  
  # Handling multiple imputation configurations (k.in a vector)
  if (nK > 1) {
    for (jj in 2:nK) {
      total.K <- c(total.K, rep(k.in[jj], reps.in[jj]))
    }
  }
  JJ <- length(total.K)
  
  # Initializing storage matrices for parameter estimates.
  a.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  alpha.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  b.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  beta.store <- matrix(rep(0, (JJ + 1) * J), nrow = (JJ + 1))
  rho.store <- matrix(rep(0, (JJ + 1)), nrow = (JJ + 1))
  vartau.store <- matrix(rep(0, (JJ + 1)), nrow = (JJ + 1))
  
  # Main loop for the MCEM algorithm iterations.
  for (jj in 1:JJ) {
    a.store[jj,] <- ests.in$a
    b.store[jj,] <- ests.in$b
    alpha.store[jj,] <- ests.in$alpha
    beta.store[jj,] <- ests.in$beta
    vartau.store[jj] <- ests.in$vartau
    rho.store[jj] <- ests.in$rho
    EM.iter <- MCEM_algorithm_one_iteration(Y, logT10, 
                                            N, J, ests.in, 
                                            total.K[jj])
    ests.in <- EM.iter
  }
  
  ests.in <- EM.iter
  a.store[JJ+1,] <- ests.in$a
  b.store[JJ+1,] <- ests.in$b
  alpha.store[JJ+1,] <- ests.in$alpha
  beta.store[JJ+1,] <- ests.in$beta
  vartau.store[JJ+1] <- ests.in$vartau
  rho.store[JJ+1] <- ests.in$rho
  # Final parameter estimates.
  #mean_a = a.store[JJ,]
  #mean_b = b.store[JJ,]
  #mean_alpha = alpha.store[JJ,]
  #mean_beta = beta.store[JJ,]
  
  # Returning the final estimates.
  MCEM.ests <- list(
    a = a.store, #formerly mean_a,
    b = b.store, #mean_b,
    alpha = alpha.store, #mean_alpha,
    beta = beta.store, #mean_beta,
    numwords.p = N,
    vartau = vartau.store,
    rho = rho.store
  )
  
  return(MCEM.ests)
}
