# GMM with Identity -----------------------------------------------------------------------------------------------
# Function that estimates GMM
gmm <- function(
  moment_cond,      # Function with arguments (theta, matrix), that returns moment conditions
  data,             # Matrix of data
  theta_0,          # Initial guess
  method = "BFGS",  # Default method for optimization
  matrix = F,
  ...
){
  # Objective function to be minimized.
  obj_func <- function(theta_0){
    mom    <- moment_cond(theta = theta_0, X = data, matrix = matrix)
    t(mom) %*% diag(length(mom)) %*% mom
  }
  # Optimization
  opt <- optim(theta_0, obj_func, method = method, ...)
  return(
    list(
      theta = opt$par,
      "Function at minimum" = opt$value
    )
  )
}

# GMM with Identity: Coordinate Descent ---------------------------------------------------------------------------
# Function that estimates GMM
# Returns:
# Optimimum parameters.
# Objective function at optimimum parameters.
gmm_coord <- function(
  moment_cond,      # Function with arguments (theta, matrix), that returns moment conditions.
  data,             # Matrix of data.
  theta_0,          # Initial guess.
  interval,         # Vector contains.
  nsteps = NULL,    # Steps for the for to run.
  eps = 1e-10,      # Coondition to check minimum.
  matrix = F,
  ...
){
  # Stop if nsteps is not an integer.
  if(!(is.integer(nsteps) | is.null(nsteps))){
    print("nsteps is not an integer.")
    stop()
  }
  # Objective function to be minimized.
  obj_func <- function(theta_0, i){
    # For Coordinate Descent we need that the function 
    # only depends of one variable.
    obj_func_without_theta <- function(theta){
      theta_0[i] <- theta
      mom        <- moment_cond(theta_0, data, matrix = matrix)
      func_value <- t(mom) %*% diag(length(mom)) %*% mom
      return(func_value)
    }
    return(obj_func_without_theta)
  }
  # Checks for steps or while approach.
  if(!is.null(nsteps) & is.integer(nsteps)){
    # Steps of optimization.
    for(j in 1:nsteps){
      # Optimization for every parameter.
      for(i in 1:length(theta_0)){
        # Optimize value theta[i]
        theta_0[i] <- optimize(obj_func(theta_0, i), interval = interval, ...)$minimum
      }
    }
  } else {
    # Create condition for the while to start with.
    condition <- obj_func(theta_0,1)(theta_0[1])
    # Condition that checks the changes within 
    # iterations.
    while(condition >= eps){
      # This helps with the difference of objective function
      # between iterations of optimization.
      helper <- obj_func(theta_0,1)(theta_0[1])
      # Optimization for every parameter.
      for(i in 1:length(theta_0)){
        theta_0[i] <- optimize(obj_func(theta_0, i), interval = interval, ...)$minimum
      }
      # Computes difference of objective function.
      condition <- abs(helper - obj_func(theta_0,1)(theta_0[1]))
    }
  }
  # Saves minimum value.
  min_func <- obj_func(theta_0, 1)(theta_0[1]) %>% as.numeric
  # Value to return.
  return(
    list(
      theta = theta_0, # Saves optimimum theta.
      "Function at minimum" = min_func # Saves value of objective function.
    )
  )
}

# Adaptive Lasso GMM ----------------------------------------------------------------------------------------------

gmm_alasso <- function(
  known_cond,       # Function with arguments (theta, matrix), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, matrix), that returns unknown moment conditions.
  data,             # Matrix of data.
  theta_0,          # Initial guess.
  lambda,           # Penalization Parameter.
  nsteps = NULL,    # Steps for the for to run.
  eps = 1e-8,      # Coondition to check minimum.
  ...
){
  # Stop if nsteps is not an integer.
  if(!(is.integer(nsteps) | is.null(nsteps))){
    print("nsteps is not an integer.")
    stop()
  }
  # 1st Step:
  # Let's bind conditions.
  moment_cond <- function(known_cond, unknown_cond){
    function(theta_0, data, matrix = T){
      moments <- cbind(
        known_cond(
          theta = theta_0,
          df = data,
          matrix = matrix
        ),
        unknown_cond(
          theta = theta_0,
          df = data,
          matrix = matrix
        )
      )
      return(moments)
    }
  }
  # Let's get the first step estimation of beta.
  first_gmm <- gmm::gmm(
    g = moment_cond(known_cond, unknown_cond),
    x = data, 
    t0 = theta_0
  )
  # Values of condition in first estimation
  first_condition <- first_gmm$coefficients %>% 
    unknown_cond(df = data, matrix = T) %>% 
    as.matrix
  n_unkown <- first_condition %>% ncol
  # Moment condition
  first_beta <- first_condition %>% 
    apply(MARGIN = 2, mean) %>% 
    as.vector()
  # Maximium and miminum values 
  first_beta_max <- first_condition %>% 
    apply(MARGIN = 2, max) %>% 
    as.vector()
  first_beta_min <- first_condition %>% 
    apply(MARGIN = 2, min) %>% 
    as.vector()
  # First step theta.
  theta <- first_gmm$coefficients
  # Standar error of theta
  theta_se <- first_gmm %>% 
    broom::tidy() %>% 
    dplyr::distinct(std.error) %>% 
    dplyr::pull()
  min_theta <- theta - 8 * theta_se
  max_theta <- theta + 8 * theta_se 
  # Second Step:
  # Stack parameters
  min_beta <- c(min_theta, first_beta_min) 
  beta     <- c(rep(0, length(theta)), first_beta)
  max_beta <- c(max_theta, first_beta_max)
  # Objective function to be minimized.
  obj_func <- function(beta, i){
    # For Coordinate Descent we need that the function 
    # only depends of one variable.
    obj_func_without_theta <- function(theta){
      # Update value
      beta[i] <- theta
      # Computes known and unknown moments
      mom <- rbind(
        known_cond(beta[1:length(theta_0)], data, matrix = F),
        unknown_cond(beta[1:length(theta_0)], data, matrix = F) - beta[(length(theta_0) + 1):length(beta)]
      )
      # Minimization problem
      func_value <- t(mom) %*% diag(length(mom)) %*% mom + # GMM
        # Penalization
        lambda * sum(abs(beta[(length(theta_0) + 1):length(beta)] / first_beta)) 
      return(func_value)
    }
    return(obj_func_without_theta)
  }
  # print("Start Optimization")
  # Checks for steps or while approach.
  if(!is.null(nsteps) & is.integer(nsteps)){
    # Steps of optimization.
    for(j in 1:nsteps){
      # Optimization for every parameter.
      for(i in 1:length(beta)){
        # Optimize value theta[i]
        beta[i] <- optimize(obj_func(beta, i), interval = c(min_beta[i], max_beta[i]),...)$minimum
      }
    }
  } else {
    # Create condition for the while to start with.
    condition <- obj_func(beta, 1)(beta[1]) 
    # Condition that checks the changes within 
    # iterations.
    while(condition >= eps){
      # This helps with the difference of objective function
      # between iterations of optimization.
      helper <- obj_func(beta,1)(beta[1])
      # Optimization for every parameter.
      for(i in 1:length(beta)){
        beta[i] <- optimize(obj_func(beta, i), interval = c(min_beta[i], max_beta[i]))$minimum
      }
      # Computes difference of objective function.
      condition <- abs(helper - obj_func(beta,1)(beta[1])) 
      # print(condition)
    }
  }
  
  gc()
  
  return(
    list(
      parameters  = beta[1:length(theta_0)],
      tested_cond =  beta[(length(theta_0)+1):length(beta)]
    )
  )
}

# Cross Validation --------------------------------------------------------
# Helper function to estimate the error
error_generator <- function(
  gmm_alasso_list, # A gmm_alasso object
  test_data,       # Data in which the error will be evaluated
  known_cond,      # Function with arguments (theta, matrix), that returns known moment conditions.
  unknown_cond     # Function with arguments (theta, matrix), that returns unknown moment conditions.
){
  tested_cond <- gmm_alasso_list$tested_cond
  parameters  <- gmm_alasso_list$parameters
  conditions  <- if_else(near(tested_cond, 0, tol = 1e-3), TRUE, FALSE)
  
  moments <- cbind(
    known_cond(
      theta = parameters,
      df = test_data,
      matrix = TRUE
    ),
    (unknown_cond(
      theta = parameters,
      df = test_data,
      matrix = TRUE
    ) - tested_cond)
  ) %>% 
    apply(
      MARGIN = 2,
      FUN = sum
    )
  
  error <- as.numeric(t(moments) %*% diag(length(moments)) %*% moments)
  
  return(list(error = error, selected = conditions, moments = moments))
}



train_test_gmm_alasso <- function(
  known_cond,       # Function with arguments (theta, matrix), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, matrix), that returns unknown moment conditions.
  train_data,       # Matrix of data.
  test_data,
  theta_0,          # Initial guess.
  lambda = NULL,    # Penalization Parameter.
  nsteps = NULL,    # Steps for the for to run.
  eps = 1e-8,       # Coondition to check minimum.
  ...
){
  # Determines which lambdas are going to be evaluated
  if(is.null(lambda)){
    lambda <- seq(-.5, 1.5, by = .01)
  }
  
  # Returns a list of parameters
  estimation <- purrr::map(
    lambda,
    ~gmm_alasso(
      known_cond = known_cond,
      unknown_cond = unknown_cond,
      data = train_data,
      theta_0 = theta_0,
      lambda = .,
      nsteps = nsteps,
      eps = eps
    )
  )
  
  error_helper <- purrr::map(
    estimation,
    ~error_generator(
      .,
      test_data = test_data,
      known_cond = known_cond,
      unknown_cond = unknown_cond
    )
  ) 
  
  results <- tibble::tibble(
    lambda = lambda
  ) %>%
    dplyr::mutate(
      selected_number  = purrr::map(error_helper, ~sum(.$selected)) %>% purrr::flatten_int(),
      selected_moments = purrr::map(error_helper, ~.$selected),
      error            = purrr::map(error_helper, ~.$error) %>% purrr::flatten_dbl(),
      moment_values    = purrr::map(error_helper, ~.$moments)
    )
  
  return(results)
  
}

cv_gmm_alasso <- function(
  nfolds,
  known_cond,       # Function with arguments (theta, matrix), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, matrix), that returns unknown moment conditions.
  data,             # Matrix of data.
  theta_0,          # Initial guess.
  lambda = NULL,    # Penalization Parameter.
  nsteps = NULL,    # Steps for the for to run.
  eps = 1e-8,       # Coondition to check minimum.
  ...
){
  
  divide <- nrow(data) / nfolds
  data <- data %>% 
    mutate(helper = rep(sample(1:nfolds), divide))
  
  resul <- purrr::map_df(
    seq(1:nfolds),
    ~train_test_gmm_alasso(
      known_cond = known_cond, 
      unknown_cond = unknown_cond, 
      train_data = dplyr::filter(data, helper != .), 
      test_data = dplyr::filter(data, helper == .),
      theta_0 = theta_0,
      lambda = lambda, 
      nsteps = nsteps, 
      eps = eps  
    )
  )
  
  return(resul)
  
}
