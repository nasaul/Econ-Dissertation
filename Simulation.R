library(MASS)
library(tibble)
library(magrittr)
library(dplyr)
library(ggplot2)

# Simulate Endogeneity ----------------------------------------------------
source("./Functions/simulation_functions.R")
df <- gen_data(
  cl = 0.2
)

# Moment Definition -------------------------------------------------------
# Define moments that you know are equal to 0.
known_conditions <- function(theta, df, matrix = T){
  X <- as.matrix(df)
  eps <- X[,"Y"] - ( X[, c("X", "W1", "W2")] %*% as.matrix(theta) )
  if(matrix == T){
    moment <- X[, c("Z1", "W1", "W2")] * rep(c(eps),3) / length(eps)
  } else{
    moment <- (t(X[, c("Z1", "W1", "W2")] ) %*% eps) / length(eps) 
  }
  return(moment)
}
# Define moments that you want to test.
unknown_conditions <- function(theta, df, matrix = T){
  X <- as.matrix(df)
  eps <- X[,"Y"] - ( X[, c("X", "W1", "W2")] %*% as.matrix(theta) )
  if(matrix == T){
    moment <- as.matrix(
      X[, c(paste("Z", 21:25, sep = ""), paste("F", 1:20, sep = ""))] * 
        rep(c(eps),5)
      )
  } else{
    moment <- (t(X[,  c(paste("Z", 21:25, sep = ""), paste("F", 1:20, sep = ""))] ) %*% eps) / length(eps) 
  }
  return(moment)
}


# Shrinkage GMM -----------------------------------------------------------
source("./Functions/gmm_misc_functions.R")

# Cross Validation
cv.model <- cv_gmm_alasso(
  nfolds = 5,
  known_cond = known_conditions,
  unknown_cond = unknown_conditions,
  data        = df,
  theta_0     = c(0, 0),
  eps = 1e-4
)

# cv.model %>% 
#   group_by(lambda) %>% 
#   summarise(error = mean(error), error_max = max(error), error_min = min(error)) %>% 
#   ggplot(aes(x = lambda)) +
#   geom_line(aes(y = error)) +
#   geom_line(aes(y = error_max), colour = "blue") +
#   geom_line(aes(y = error_min), colour = "blue")


# Alasso GMM
gmm_alasso(
  known_cond = known_conditions,
  unknown_cond = unknown_conditions,
  data        = df,
  theta_0     = c(0, 0),
  lambda = 0.025,
  eps = 1e-8
)