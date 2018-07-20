# Optimizer helper ------------------------------------------------------------------------------------------------
# Returns: function with only one parameter (useful for coordinate descent optimization).
optim_help <- function(
  i,               # Index of parameter to be altered.
  model_param,     # Indicates if the update is going to be to theta o beta.
  known_cond,      # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,    # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  theta,           # Model parameters.
  lambda,          # Penalization hyperparameter
  beta,            # Moments parameters.
  first_beta,      # First estimation of moments parameters.
  data             # data.frame to be used with conditions functions.
){
  
  min_func <- function(alpha){
    # Checks which parameter is going to be altered
    if(model_param == TRUE){
      theta[i] <- alpha
    } else {
      beta[i]  <- alpha
    }
    # Error evaluation
    error <- train_error(
      known_cond   = known_cond,
      unknown_cond = unknown_cond,
      theta        = theta,
      lambda       = lambda,
      beta         = beta,
      first_beta   = first_beta,
      data         = data
    )
    return(error)
  }
  return(min_func)
}
# Optimizer -------------------------------------------------------------------------------------------------------
# Returns: optimized parameter.
optim_gmm <- function(
  i,               # Index of parameter to be altered.
  model_param,     # Indicates if the update is going to be to theta o beta.
  min_p,           # Minimum value of parameter.
  max_p,           # Maximum value of parameter
  known_cond,      # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,    # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  theta,           # Model parameters.
  lambda,          # Penalization hyperparameter
  beta,            # Moments parameters.
  first_beta,      # First estimation of moments parameters.
  data             # data.frame to be used with conditions functions.
){
  
  optim_val <- optimize(
    optim_help(
      i            = i,
      model_param  = model_param,
      known_cond   = known_cond ,
      unknown_cond = unknown_cond,
      theta        = theta, 
      lambda       = lambda,
      beta         = beta,
      first_beta   = first_beta,
      data         = data
    ),
    interval = c(min_p, max_p),
    maximum = FALSE
  )$minimum
  
  return(optim_val)
}

# Error evaluation ------------------------------------------------------------------------------------------------
# Two different errors, two different functions.
# Returns: penalized error.
train_error <- function(
  known_cond,       # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  theta,            # Model Parameter.
  lambda,           # Penalization hyperparameter.
  beta,             # Moments parameters.
  first_beta,       # First estimation of moments parameters.
  data              # data.frame to be used with conditions functions.
){
  # Create a matrix of beta (important when using matrix-vector operations).
  N <- nrow(data)
  beta_mat <- matrix(rep(1, N), ncol = 1) %*% beta
  # Use all moments
  moments <- cbind(
    known_cond(theta, data),
    (unknown_cond(theta, data) - beta_mat)
  ) %>% 
    apply(
      MARGIN = 2,
      FUN = mean
    )
  # Penalized error.
  error <- t(moments) %*% diag(length(moments)) %*% moments + lambda * sum(abs(beta / first_beta))
  return(as.numeric(error))
}

# Returns: unpenalized error.
test_error <- function(
  known_cond,       # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  theta,            # Model Parameter.
  beta,             # Moments parameters.
  data              # data.frame to be used with conditions functions.
){
  # Create a matrix of beta (important when using matrix-vector operations).
  N <- nrow(data)
  beta_mat <- matrix(rep(1, N), ncol = 1) %*% beta
  # Use all moments
  moments <- cbind(
    known_cond(theta, data),
    (unknown_cond(theta, data) - beta_mat)
  ) %>% 
    apply(
      MARGIN = 2,
      FUN = mean
    )
  # Unpenalized error.
  error <- t(moments) %*% diag(length(moments)) %*% moments
  return(as.numeric(error))
}

# GMM Lasso -------------------------------------------------------------------------------------------------------
# Returns: two element list.
#        model_parameters:  Model parameters.
#        moment_parameters: Moments parameters.
gmm_lasso <- function(
  known_cond,       # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  data,             # data.frame to be used with conditions functions.
  theta_0,          # Initial guess of model parameter.
  lambda,           # Penalization hyperparameter.
  nsteps = NULL,    # Maximum steps for optimization.
  eps = 1e-8        # Coondition to check minimum.
){
  # Checks if iteration steps is an integer -------------------------------------------------------------------------
  if(!(is.integer(nsteps) | is.null(nsteps))){
    print("nsteps is not an integer, defaults to 1000.")
    nsteps = 1000L
  }
  # First Step Estimation -------------------------------------------------------------------------------------------
  first_estimation <- gmm::gmm(
    g  = known_cond,
    x  = data,
    t0 = theta_0
  )  
  # Save first beta estimation.
  first_beta <- unknown_cond(
    first_estimation$coefficients,
    data
  ) %>% 
    apply(
      MARGIN = 2,
      FUN    = mean
    )
  # Parameter definition --------------------------------------------------------------------------------------------
  # Model parameters
  K         <- length(theta_0) 
  theta     <- first_estimation$coefficients
  min_theta <- rep(-10, K)
  max_theta <- rep( 10, K)
  # Moments conditions parameters
  M         <- length(first_beta)
  beta      <- first_beta
  min_beta  <- rep(-5, M)
  max_beta  <- rep( 5, M)

  # Optimization ----------------------------------------------------------------------------------------------------
  # Initialize stop criterions.
  condition <- 1
  j         <- 0 
  while(condition > eps & j < nsteps){
    # Update steps.
    j <- j + 1
    # Error evaluation, before a coordinate iteration.
    helper <- train_error(
      known_cond   = known_cond,
      unknown_cond = unknown_cond,
      theta        = theta,
      lambda       = lambda, 
      beta         = beta,
      data         = data,
      first_beta   = first_beta
    )
    # Model parameter optimization.
    for(Z in 1:K){
      # Updates parameter.
      theta[Z] <- optim_gmm(
        i            = Z,
        model_param  = TRUE,
        min_p        = min_theta[Z],
        max_p        = max_theta[Z], 
        known_cond   = known_cond,
        unknown_cond = unknown_cond,
        theta        = theta, 
        lambda       = lambda,
        beta         = beta,
        first_beta   = first_beta,
        data         = data
      )
    }
    # Moment conditions parameters optimization.
    for(U in 1:M){
      # Updates paramater.
      beta[U]   <- optim_gmm(
        i            = U,
        model_param  = FALSE,
        min_p        = min_beta[U],
        max_p        = max_beta[U], 
        known_cond   = known_cond,
        unknown_cond = unknown_cond,
        theta        = theta, 
        lambda       = lambda,
        beta         = beta,
        first_beta   = first_beta,
        data         = data
      )
    }
    # Check how much the penalized error changed.
    condition <- abs(
      helper - train_error(
        known_cond = known_cond,
        unknown_cond = unknown_cond,
        theta = theta,
        lambda = lambda, 
        beta = beta,
        data = data,
        first_beta = first_beta
      )
    )
  }
  
  return(
    list(
      model_parameters  = theta,
      moment_parameters = beta
    )
  )
}

# Train-Test Estimation ---------------------------------------------------------------------------------------------
# Returns: tibble with columns: error, model_parameters and moments_parameters.
train_test_info <- function(
  known_cond,       # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  train_data,       # data.frame to be used with conditions functions.
  test_data,        # data.frame to be used with conditions functions.
  theta_0,          # Initial guess of model parameter.
  lambda,           # Penalization hyperparameter.
  nsteps = NULL,    # Maximum steps for optimization.
  eps = 1e-8        # Coondition to check minimum.
){
  # Fit the model.
  parameter_estimation <- gmm_lasso(
    known_cond   = known_cond,      
    unknown_cond = unknown_cond,
    data         = train_data,
    theta_0      = theta_0,
    lambda       = lambda,
    nsteps       = nsteps,
    eps          = eps
  )
  # Estimate error.
  estimated_error <- test_error(
    known_cond   = known_cond,      
    unknown_cond = unknown_cond,
    theta        = parameter_estimation$model_parameters,
    beta         = parameter_estimation$moment_parameters,
    data         = test_data
  )
  # Model information.
  result <- tibble::tibble(
    error             = estimated_error,
    model_parameters  = list(parameter_estimation$model_parameters),
    moment_parameters = list(parameter_estimation$moment_parameters),
    lambda            = lambda
  )
  return(result)
}

# Evaluate multiple lambdas -----------------------------------------------------------------------------------------
# Returns: tibble with columns: error, model_parameters and moments_parameters for multiple lambdas.
train_test_error <- function(
  known_cond,       # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  train_data,       # data.frame to be used with conditions functions.
  test_data,        # data.frame to be used with conditions functions.
  theta_0,          # Initial guess of model parameter.
  lambdas = seq(0, 1, by = 0.01),          # Penalization hyperparameters.
  nsteps = 1000L,    # Maximum steps for optimization.
  eps = 1e-8        # Coondition to check minimum.
){
 results <- purrr::map_df(
   lambdas,
   ~train_test_info(
     known_cond   = known_cond,
     unknown_cond = unknown_cond,
     train_data   = train_data,
     test_data    = test_data,
     theta_0      = theta_0,
     lambda       = .,
     nsteps       = nsteps,
     eps          = eps
   )
 )
 return(results)
}

# Cross-Validation Estimation ---------------------------------------------------------------------------------------
# Returns: tibble with columns: error, model_parameters and moments_parameters for multiple lambdas.
cv_gmm_lasso <- function(
  nfolds,           # Number of folds for cross-validation.
  known_cond,       # Function with arguments (theta, data.frame), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, data.frame), that returns unknown moment conditions.
  data,             # data.frame to be used with conditions functions.
  theta_0,          # Initial guess of model parameter.
  lambdas = seq(0, 1, by = 0.01),          # Penalization hyperparameters.
  nsteps = NULL,    # Maximum steps for optimization.
  eps = 1e-8        # Coondition to check minimum.
){
  # Create folds.
  folded_data <- data %>% 
    dplyr::mutate(
      fold = base::sample(
        1:nfolds,
        size = nrow(data),
        replace = TRUE
      )
    )
  # Cross-validate the folds.
  cv_results <- purrr::map_df(
    1:nfolds,
    ~train_test_error(
      known_cond   = known_cond,
      unknown_cond = unknown_cond,
      train_data   = dplyr::filter(folded_data, fold != .),
      test_data    = dplyr::filter(folded_data, fold == .) ,
      theta_0      = theta_0,
      lambdas      = lambdas,
      nsteps       = nsteps,
      eps          = eps
    )
  )
  return(cv_results)
}
