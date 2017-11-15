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
    mom    <- moment_cond(theta_0, data, matrix = matrix)
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
  unknown_cond,      # Function with arguments (theta, matrix), that returns unknown moment conditions.
  data,             # Matrix of data.
  theta_0,          # Initial guess.
  lambda,           # Penalization Parameter.
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
# 1st Step:
  # Let's bind conditions.
  moment_cond <- function(known_cond, unknown_cond){
    function(theta_0, data, matrix = matrix){
      moments <- rbind(known_cond(theta_0, data, matrix), unknown_cond(theta_0, data, matrix))
    }
  }
  # Let's get the first step estimation of beta.
  first_gmm <- gmm(
    moment_cond = moment_cond(known_cond, unknown_cond),
    data        = data, 
    theta_0     = theta_0
  )
  # First step beta.
  first_beta <- first_gmm$theta %>% 
    unknown_cond(X = data, matrix = matrix) %>% 
    as.numeric
  # First step theta.
  theta <- first_gmm$theta
# Second Step:
  # Stack parameters
  beta <- c(theta, first_beta)
  # Objective function to be minimized.
  obj_func <- function(beta, i){
    # For Coordinate Descent we need that the function 
    # only depends of one variable.
    obj_func_without_theta <- function(theta){
      # Update value
      beta[i] <- theta
      # Computes known and unknown moments
      mom        <- known_cond(beta[1:length(theta_0)], data, matrix) %>% 
        rbind(., unknown_cond(beta[1:length(theta_0)], data, matrix) - rep(beta[(length(theta_0) + 1):length(beta)], length(first_beta) ))
      # Minimization problem
      func_value <- t(mom) %*% diag(length(mom)) %*% mom + # GMM
        # Penalization
        lambda * sum(abs(beta[(length(theta_0) + 1):length(beta)] / first_beta)) 
      return(func_value)
    }
    return(obj_func_without_theta)
  }
# Checks for steps or while approach.
  if(!is.null(nsteps) & is.integer(nsteps)){
    # Steps of optimization.
    for(j in 1:nsteps){
      # Optimization for every parameter.
      for(i in 1:length(beta)){
        # Optimize value theta[i]
        beta[i] <- optimize(obj_func(beta, i), interval = interval,...)$minimum
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
        beta[i] <- optimize(obj_func(beta, i), interval = interval, ...)$minimum
      }
      # Computes difference of objective function.
      condition <- abs(helper - obj_func(beta,1)(beta[1]))
    }
  }
 return(beta)
}
