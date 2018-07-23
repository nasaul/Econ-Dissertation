# Function that generates IV data 
# Experiment from paper:
# Name    : Adaptive GMM Shrinkage Estimation with Consistent Moment Selection
# Author  : Zhipeng Liao

gen_data <- function(
  cl,                                         # cl \in [0, 0.8]
  theta = c(0.5, 0.4, 0.3),                   # theta 1 x 3
  pi1   = 0.2,                                # pi1   1 x 1
  pi2   = c(0.1, 0.1),                        # pi2   1 x 2
  pi3   = c(0.15, 0.15, 0.2, 0.2, 0.25),      # pi3   1 x 5
  N     = 1000L # N 
){
  
  # Data Generation
  sigma <- diag(
    nrow = 30,
    ncol = 30
  )
  
  for(i in 1:8){
    for(j in 1:8){
      sigma[i, j] = 0.2 ^ abs(i-j)
    }
  }
  
  
  sigma[9:10, 9:10] <- c(1, 0.6, 0.6, 1)
  
  help_matrix <- MASS::mvrnorm(
    n = N,
    mu = rep(0, 30),
    Sigma = sigma
  )
  
  # Creation of fake instruments
  
  L <- rep(0, 20)
  for(j in 1:20){
    L[j] <-  cl + (0.8 - cl) * (j - 1) / 19
  }
  
  F_mat <- help_matrix[, 11:30] + help_matrix[, 9] %*% t(L)
  
  # Data generation
  gen_data <- tibble::tibble(
    W1  = help_matrix[, 1],
    W2  = help_matrix[, 2],
    Z1  = help_matrix[, 3],
    Z21 = help_matrix[, 4],
    Z22 = help_matrix[, 5],
    Z23 = help_matrix[, 6],
    Z24 = help_matrix[, 7],
    Z25 = help_matrix[, 8],
    U   = help_matrix[, 9],
    V   = help_matrix[, 10],
    F1  = F_mat[, 1],
    F2  = F_mat[, 2],
    F3  = F_mat[, 3],
    F4  = F_mat[, 4],
    F5  = F_mat[, 5],
    F6  = F_mat[, 6],
    F7  = F_mat[, 7],
    F8  = F_mat[, 8],
    F9  = F_mat[, 9],
    F10 = F_mat[, 10],
    F11 = F_mat[, 11],
    F12 = F_mat[, 12],
    F13 = F_mat[, 13],
    F14 = F_mat[, 14],
    F15 = F_mat[, 15],
    F16 = F_mat[, 16],
    F17 = F_mat[, 17],
    F18 = F_mat[, 18],
    F19 = F_mat[, 19],
    F20 = F_mat[, 20]
  ) %>%  
    dplyr::mutate(
      X = Z1 * pi1 +
        W1  * pi2[1] +
        W2  * pi2[2] + 
        Z21 * pi3[1] +
        Z22 * pi3[2] +
        Z23 * pi3[3] +
        Z24 * pi3[4] +
        Z25 * pi3[5] +
        V,
      Y = X * theta[1] +
        W1 * theta[2] +
        W2 * theta[3] + 
        U
    )
  
  return(gen_data)
}


# Moment Definition -----------------------------------------------------------------------------------------------
# Define moments that you know are equal to 0.
known_conditions <- function(theta, df){
  X <- as.matrix(df)
  eps <- X[,"Y"] - ( X[, c("X", "W1", "W2")] %*% as.matrix(theta) )
  
  moment <- X[, c("Z1", "W1", "W2")] * rep(c(eps),3) / length(eps)
  
  return(moment)
}
# Define moments that you want to test.
unknown_conditions <- function(theta, df){
  X <- as.matrix(df)
  eps <- X[,"Y"] - ( X[, c("X", "W1", "W2")] %*% as.matrix(theta) )
  
  moment <- as.matrix(
    X[, c(paste("Z", 21:25, sep = ""), paste("F", 1:20, sep = ""))] * 
      rep(c(eps),25)
  )
  
  return(moment)
}


# Simulate and Estimate parameters --------------------------------------------------------------------------------

simulation <- function(cl, lambdas, eps, N, pi1){
  df <- gen_data(
    cl  = cl,
    N   = N,
    pi1 = pi1
  )
  model <- cv_gmm_lasso(
    nfolds       = 2,
    known_cond   = known_conditions,
    unknown_cond = unknown_conditions,
    data         = df,
    theta_0      = c(0, 0, 0),
    eps          = eps,
    nsteps       = 1000L,
    lambdas      = lambdas
  )
  # Extract lambda that achieves minimum error.
  min_lambda <- model %>% 
    group_by(lambda) %>%
    summarise(error = mean(error)) %>%
    filter(error == min(error)) %>% 
    pull(lambda)
  # Estimation of real model
  estimation <- gmm_lasso(
    known_cond   = known_conditions,
    unknown_cond = unknown_conditions,
    data         = df,
    theta_0      = c(0, 0, 0),
    eps          = 1e-8,
    nsteps       = 1000L,
    lambda       = min_lambda
  )
  # Comparates if moments are true or false.
  mom_comp <- near(estimation$moment_parameters, 0, tol = 1e-4)
  # Selected moments.
  selected_mom <- names(subset(mom_comp, mom_comp))
  return(
    list(
      selected_mom = selected_mom,
      mom_par      = estimation$moment_parameters
    )
    
  )
}
