#library(gmm)
library(magrittr)
library(dplyr)

# Simulate Endogeneity ----------------------------------------------------
# Model:
# Y = theta * X_1 + theta_1 * Z_1 + theta_2 * Z_2 + theta_3 * X_2 eps
# X = lambda_1 * Z_1 + lambda_2 * Z_2 + u
# Parameters
  n        <- 10000
  theta    <- 5
  theta_1  <- 0.5
  theta_2  <- 0.5
  theta_3  <- 0.7
  lambda_1 <- 3
  lambda_2 <- 5
# Data generation
set.seed(158684)

data_matrix <- data_frame(
  z_1 = rnorm(n = n, mean = 0, sd = 20),
  z_2 = rnorm(n = n, mean = 0, sd = 30),
  x_1   = rnorm(n = n,
              mean = lambda_1 * z_1 + lambda_2 * z_2,
              sd = 1),
  x_2 = rnorm(n = n, mean = 0, sd = 10),
  y   = rnorm(n = n, 
              mean = theta * x_1 + theta_1 * z_1 + theta_2 * z_2 + theta_3 * x_2,
              sd = 1)
) %>% 
  as.matrix

# Define moments ----------------------------------------------------------
# Notice that this always be user input,
# for simulation purposes I can do whatever I want.
moment_conditions <- function(theta, X, matrix = T){
  if(matrix == T){
    eps <- X[,"y"] - ( X[, c("x_1", "x_2")] %*% theta )
    moment <- X[, c("z_1","z_2", "x_2")] * rep(c(eps),3) / length(eps)
  } else{
    eps_ <- X[,"y"] - ( X[, c("x_1", "x_2")] %*% theta )
    moment <- t(X[, c("z_1","z_2", "x_2")] ) %*% eps_ / length(eps_)
  }
  return(moment)
}

# GMM R function ----------------------------------------------------------

# normalid_gmm <- gmm(
#   moment_conditions,
#   data_matrix,
#   t0 = c(0,0),
#   wmatrix ="ident",
#   method = "BFGS"
# )
# 
# normalop_gmm <- gmm(
#   moment_conditions,
#   data_matrix,
#   type     = "twoStep",
#   t0       = c(0,0),
#   wmatrix  = "ident",
#   method   = "BFGS"
# )

# GMM  --------------------------------------------------------------------
source("./Functions/gmm_misc_functions.R")

# Normal GMM Estimation with identity matrix
id_gmm <- gmm(
  moment_cond = moment_conditions,
  data        = data_matrix,
  theta_0     = c(0,0),
  method      = "BFGS",
  matrix      = F
) 

# Coordinate Descent Estimation with identity matrix

coord_gmm <- gmm_coord(
  moment_cond = moment_conditions,
  data        = data_matrix,
  theta_0     = c(0,0),
  interval    = c(-500,500),
  matrix      = F
) 

# Efficient GMM Estimation

# ef_gmm <- gmm_eff(
#   moment_cond = moment_conditions,
#   data        = data_matrix,
#   theta_0     = c(0, 0),
#   method = "BFGS"
# ) 

# Shrinkage GMM ---------------------------------------------------------------------------------------------------
# Define moments that you know are equal to 0.
known_conditions <- function(theta, X, matrix = T){
  if(matrix == T){
    eps <- X[,"y"] - ( X[, c("x_1", "x_2")] %*% theta )
    moment <- X[, c("z_1", "x_2")] * rep(c(eps),2) / length(eps)
  } else{
    eps_ <- X[,"y"] - ( X[, c("x_1", "x_2")] %*% theta )
    moment <- (t(X[, c("z_1", "x_2")] ) %*% eps_) / length(eps_) 
  }
  return(moment)
}
# Define moments that you want to test.
unknown_conditions <- function(theta, X, matrix = T){
  if(matrix == T){
    eps <- X[,"y"] - ( X[, c("x_1", "x_2")] %*% theta )
    moment <- X[, c("z_2")] * rep(c(eps),1) 
  } else{
    eps_ <- X[,"y"] - ( X[, c("x_1", "x_2")] %*% theta )
    moment <- (t(X[, c("z_2")] ) %*% eps_) / length(eps_) 
  }
  return(moment)
}
# Alasso GMM
gmm_alasso(
  known_cond = known_conditions,
  unknown_cond = unknown_conditions,
  data        = data_matrix,
  theta_0     = c(0,0),
  lambda = 10000,
  eps = 1e-8
)
