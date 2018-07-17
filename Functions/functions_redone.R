# Optimizer helper ------------------------------------------------------------------------------------------------
optim_help <- function(
  i,
  model_param,
  known_cond,
  unknown_cond,
  theta, 
  lambda,
  beta,
  first_beta,
  data
){
  
  min_func <- function(alpha){
    
    if(model_param == TRUE){
      theta[i] <- alpha
    } else {
      beta[i]  <- alpha
    }
    
    error <- error_eval(
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
optim_gmm <- function(
  i,
  model_param,
  min_p,
  max_p, 
  known_cond,
  unknown_cond,
  theta, 
  lambda,
  beta,
  first_beta,
  data
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
error_eval <- function(
  known_cond,
  unknown_cond,
  theta,
  lambda,
  beta,
  first_beta,
  data
){
  
  moments <- cbind(
    known_cond(theta, data),
    unknown_cond(theta, data) - beta
  ) %>% 
    apply(
      MARGIN = 2,
      FUN = mean
    )
  
  error <- t(moments) %*% diag(length(moments)) %*% moments + lambda * sum(abs(beta / first_beta))
  return(as.numeric(error))
}
# GMM Lasso -------------------------------------------------------------------------------------------------------
gmm_lasso <- function(
  known_cond,       # Function with arguments (theta, matrix), that returns known moment conditions.
  unknown_cond,     # Function with arguments (theta, matrix), that returns unknown moment conditions.
  data,             # Matrix of data.
  theta_0,          # Initial guess.
  lambda,           # Penalization Parameter.
  nsteps = NULL,    # Steps for the for to run.
  eps = 1e-8,       # Coondition to check minimum.
  ...
){

  # Checks if iteration steps is an integer -------------------------------------------------------------------------
  if(!(is.integer(nsteps) | is.null(nsteps))){
    print("nsteps is not an integer.")
    stop()
  }

  # First Step Estimation -------------------------------------------------------------------------------------------
  
  first_estimation <- gmm::gmm(
    g  = known_cond,
    x  = data,
    t0 = theta_0
  )  
  
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
  min_theta <- rep(-1000, K)
  max_theta <- rep( 1000, K)
  # Moments conditions parameters
  M         <- length(first_beta)
  beta      <- first_beta
  min_beta  <- rep(-1000, M)
  max_beta  <- rep( 1000, M)

  # Optimization ----------------------------------------------------------------------------------------------------

  condition <- 1
  
  while(condition > eps){
    
    helper <- error_eval(
      known_cond = known_cond,
      unknown_cond = unknown_cond,
      theta = theta,
      lambda = lambda, 
      beta = beta,
      data = data,
      first_beta = first_beta
    )
    
    # Model parameter optimization
    for(Z in 1:K){
      theta[Z] <- optim_gmm(
        i = Z,
        model_param = TRUE,
        min_p = min_theta[Z],
        max_p = max_theta[Z], 
        known_cond = known_cond,
        unknown_cond = unknown_cond,
        theta = theta, 
        lambda = lambda,
        beta = beta,
        first_beta  = first_beta,
        data = data
      )
    }
    # Moment conditions parameters optimization
    for(U in 1:M){
      beta[U]   <- optim_gmm(
        i = U,
        model_param = FALSE,
        min_p = min_beta[U],
        max_p = max_beta[U], 
        known_cond = known_cond,
        unknown_cond = unknown_cond,
        theta = theta, 
        lambda = lambda,
        beta = beta,
        first_beta  = first_beta,
        data = data
      )
    }
    # Update condition
    condition <- abs(
      helper - error_eval(
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
      parameters  = theta,
      tested_cond = beta
    )
  )
}