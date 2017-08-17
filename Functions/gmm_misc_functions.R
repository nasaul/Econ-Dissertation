# GMM with Identity -------------------------------------------------------
  # Function that stimates GMM
gmm_1 <- function(
  moment_cond, # Function with arguments (theta, matrix), that returns moment conditions
  data,        # Matrix of data
  theta_0,      # Initial guess
  method = "BFGS",
  matrix = F,
  ...
  ){
  
  obj_func <- function(theta_0){
    mom <- moment_cond(theta_0, data, matrix = matrix)
    t(mom) %*% diag(length(mom)) %*% mom
  }
  opt <- optim(theta_0, obj_func, method = method, ...)
  moments <- moment_cond(opt$par,data)
  
  return(
    list(
      "Optimization Result" = opt,
      "Moment Conditions"   = moments
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
