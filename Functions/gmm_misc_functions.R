# GMM with Identity -------------------------------------------------------
# Function that estimates GMM
gmm_1 <- function(
  moment_cond,      # Function with arguments (theta, matrix), that returns moment conditions
  data,             # Matrix of data
  theta_0,          # Initial guess
  method = "BFGS",  # Default method for optimization
  matrix = F,
  ...
  ){
  
  obj_func <- function(theta_0){
    mom <- moment_cond(theta_0, data, matrix = matrix)
    t(mom) %*% diag(length(mom)) %*% mom
  }
  opt <- optim(theta_0, obj_func, method = method, ...)
  
  
  return(
    list(
      "Optimization Result" = opt$par,
      "Function at minimum" = opt$value
    )
  )
}

# GMM with Identity: Coordinate Descent ---------------------------------------------------------------------------

# Function that estimates GMM
gmm_coord <- function(
  moment_cond,      # Function with arguments (theta, matrix), that returns moment conditions
  data,             # Matrix of data
  theta_0,          # Initial guess
  interval,         # Vector contains 
  nsteps = NULL,    # nsteps for the for to run
  eps = 1e-10,      # Coondition to check minimum
  matrix = F,
  ...
){
  stopifnot(is.integer(nsteps) | is.null(nsteps))
  obj_func <- function(theta_0, i){
    obj_func_without_theta <- function(theta){
      theta_0[i] <- theta
      mom <- moment_cond(theta_0, data, matrix = matrix)
      t(mom) %*% diag(length(mom)) %*% mom
    }
  }
  
  theta <- c()
  
  if(is.null(nsteps)){
    condition <- obj_func(theta_0,1)(1)
  }
  
  if(!is.null(nsteps) & is.integer(nsteps)){
    # Cycles of optimization
    for(j in 1:nsteps){
      # Optimization for every parameter
      for(i in 1:length(theta_0)){
        result    <- optimize(obj_func(theta_0, i), interval = interval, ...)
        theta[i]  <- result$minimum
      }
    }
  } else {
    while(condition >= eps ){
      # Optimization for every parameter
      for(i in 1:length(theta_0)){
        result    <- optimize(obj_func(theta_0, i), interval = interval, ...)
        theta[i]  <- result$minimum
      }
      # Checks condition
        condition <- abs(condition - obj_func(theta,1)(1))
    }
  }
  
  min_func<- obj_func(theta, 1)(theta[1]) %>% as.numeric
  
  return(
    list(
      "Optimization Result" = theta,
      "Function at minimum" = min_func
    )
  )
}

# Efficient GMM -----------------------------------------------------------
  # Function that uses two-step optimization of GMM
    # 1 Step. GMM with matrix identity
    # 2 Step  GMM with efficient weight matrix
gmm_eff  <- function(
  moment_cond, # Function with arguments (theta, matrix), that returns moment conditions
  data,        # Matrix of data
  theta_0,     # Initial guess
  method = "BFGS",
  ...
){
# Get efficient Weight Matrix
  gmm_i <- gmm_1(moment_cond, data, theta_0, matrix = F, method = method, ...)
  W <- t(gmm_i$`Moment Conditions`) %*% gmm_i$`Moment Conditions` %>% solve
  
# Define function to minimize
  obj_func <- function(theta_0){
    mom <- moment_cond(theta_0, data, matrix = F)
    t(mom) %*% W %*% mom
  }
# Optimization
  opt <- optim(theta_0, obj_func, method = method, ...)
  
  return(
    list(
      "Optimization Result" = opt
    )
  )
}
