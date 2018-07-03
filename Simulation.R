library(MASS)
library(tibble)
library(magrittr)
library(dplyr)

# Simulate Endogeneity ----------------------------------------------------
# Model:
# (xstar, c) ~ normal((0,0), (1, rho, rho, 1)) ->
#   cor(xstar, c) = rho
# x <- xstar + theta[1] * z1 + theta[2] * z2 + theta[3] * z3 + theta[4] * z4 + theta[5] * z5
# y ~ normal( 
#   alpha[1] + alpha[2] * x + alpha[3] c,
#   sigma_sd^2
# )
# z6 ... z10 are confounders

# Parameters
  N <- 500
  rho <- 0.9
  theta_param <- c(
    0.9, # theta1
    0.8, # theta2
    0.7, # theta3 
    0.6, # theta4 
    0.5  # theta5
  )
  alpha <- c(
    1, # Intercept
    2, # Parameter of interest
    3,  # alpha3
    4
  )
  sigma_sd <- 5
# Data generation
set.seed(158684)

df <- mvrnorm(
  n = N,
  mu = c(0,0),
  Sigma =  matrix(c(1,rho,rho,1), 2, 2)
) %>% 
  as_tibble() %>% 
  rename(
    xstar = V1,
    c = V2
  ) %>% 
  mutate(
    z1  = rnorm(N),
    z2  = rnorm(N),
    z3  = rnorm(N),
    z4  = rnorm(N),
    z5  = rnorm(N),
    z6  = rnorm(N),
    z7  = rnorm(N),
    z8  = rnorm(N),
    z9  = rnorm(N),
    z10 = rnorm(N),
    x  = xstar + theta_param[1] * z1 + theta_param[2] * z2 + theta_param[3] * z3 + theta_param[4] * z4 + theta_param[5] * z5,
    x1 = rnorm(N),
    y = alpha[1] + alpha[2] * x + alpha[3] * c + alpha[4] * x1 +rnorm(N, 0 , sigma_sd)
  )

# Moment Definition -------------------------------------------------------
# Define moments that you know are equal to 0.
known_conditions <- function(theta, df, matrix = T){
  X <- as.matrix(df)
  eps <- X[,"y"] - ( X[, c("x", "x1")] %*% as.matrix(theta) )
  if(matrix == T){
    moment <- X[, c("z1", "z2", "x1")] * rep(c(eps),3) / length(eps)
  } else{
    moment <- (t(X[, c("z1", "z2", "x1")] ) %*% eps) / length(eps) 
  }
  return(moment)
}
# Define moments that you want to test.
unknown_conditions <- function(theta, df, matrix = T){
  X <- as.matrix(df)
  eps <- X[,"y"] - ( X[, c("x", "x1")] %*% as.matrix(theta) )
  if(matrix == T){
    moment <- as.matrix(X[, c("z3", "z4", "z5", "z6", "z7")] * rep(c(eps),5))
  } else{
    moment <- (t(X[, c("z3", "z4", "z5", "z6", "z7")] ) %*% eps) / length(eps) 
  }
  return(moment)
}


# Shrinkage GMM -----------------------------------------------------------
source("./Functions/gmm_misc_functions.R")
# Alasso GMM
gmm_alasso(
  known_cond = known_conditions,
  unknown_cond = unknown_conditions,
  data        = df,
  theta_0     = c(0, 0),
  lambda = 100,
  eps = 1e-8
)

algo <- cv_gmm_alasso(
  known_cond = known_conditions,
  unknown_cond = unknown_conditions,
  data        = df,
  theta_0     = c(0, 0),
  eps = 1e-8
)

algo %>% 
  ggplot(data = ., aes(x = lambda, y = error)) + geom_line()
