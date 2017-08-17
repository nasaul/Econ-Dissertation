library(gmm)
library(magrittr)
library(dplyr)

# Simulate Endogeneity ----------------------------------------------------
# Model:
# Y = theta * X + theta_1 * Z_1 + theta_2 * Z_2 + eps
# X = lambda_1 * Z_1 + lambda_2 * Z_2 + u

# Parameters
n        <- 10000
theta    <- 5
theta_1  <- .5
theta_2  <- .5
lambda_1 <- 3
lambda_2 <- 5

# Data generation

set.seed(158684)

data_matrix <- data_frame(
  z_1 = rnorm(n = n, mean = 0, sd = 20),
  z_2 = rnorm(n = n, mean = 0, sd = 30),
  x   = rnorm(n = n,
              mean = lambda_1 * z_1 + lambda_2 * z_2,
              sd = 1),
  y   = rnorm(n = n, 
              mean = theta * x + theta_1 * z_1 + theta_2 * z_2,
              sd = 1)
) %>% 
  as.matrix

# Define moments ----------------------------------------------------------

moment_conditions <- function(theta, X, matrix = T){
  if(matrix == T){
    eps <- X[,"y"] - ( X[,"x"] * theta )
    moment <- X[, c("z_1","z_2")] * rep(c(eps),2) 
  } else{
    eps_ <- X[,"y"] - ( X[,"x"] * theta )
    moment <- t(X[, c("z_1","z_2")] ) %*% eps_
  }
  return(moment)
}

# GMM R function ----------------------------------------------------------

normalid_gmm <- gmm(
  moment_conditions,
  data_matrix,
  t0 = 0,
  wmatrix ="ident",
  method = "BFGS"
)

normalop_gmm <- gmm(
  moment_conditions,
  data_matrix,
  type     = "twoStep",
  t0       = 0,
  wmatrix  = "ident",
  method   = "BFGS"
)

# GMM  --------------------------------------------------------------------
source("./gmm_misc_functions.R")

# Normal GMM Estimation
id_gmm <- gmm_1(
  moment_cond = moment_conditions,
  data        = data_matrix,
  theta_0     = 0,
  method      = "BFGS",
  matrix      = F
) 

# Efficient GMM Estimation

ef_gmm <- gmm_eff(
  moment_cond = moment_conditions,
  data        = data_matrix,
  theta_0     = 0,
  method = "BFGS"
) 

normalid_gmm$coefficients
id_gmm$`Optimization Result`$par
normalop_gmm$coefficients
ef_gmm$`Optimization Result`$par
