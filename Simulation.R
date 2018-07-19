library(MASS)
library(tibble)
library(magrittr)
library(dplyr)
library(ggplot2)

# Simulate Endogeneity ----------------------------------------------------
source("./Functions/Simulation/simulation_function.R")
df <- gen_data(
  cl = 0.2
)

# Moment Definition -------------------------------------------------------
# Define moments that you know are equal to 0.
known_conditions <- function(theta, df){
  X <- as.matrix(df)
  eps <- X[,"Y"] - ( X[, c("X", "W1", "W2")] %*% as.matrix(theta) )
  
  moment <- X[, c("Z1", "W1", "W2")] * rep(c(eps),3) / length(eps)

  return(moment)
}
# Define moments that you want to test.
unknown_conditions <- function(theta, df){
  X <- as.matrix(df)
  eps <- X[,"Y"] - ( X[, c("X", "W1", "W2")] %*% as.matrix(theta) )
  
  moment <- as.matrix(
    X[, c(paste("Z", 21:25, sep = ""), paste("F", 1:20, sep = ""))] * 
      rep(c(eps),25)
  )
  
  return(moment)
}


# Shrinkage GMM -----------------------------------------------------------
source("./Functions/penalized_gmm_functions.R")
# Cross Validation

model <- cv_gmm_lasso(
  nfolds       = 2,
  known_cond   = known_conditions,
  unknown_cond = unknown_conditions,
  data         = df,
  theta_0      = c(0, 0, 0),
  eps          = 1e-4,
  nsteps       = 1000L,
  lambdas      = seq(0, 0.05, length.out = 100) 
)


model %>% View
    group_by(lambda) %>%
    summarise( error_max = max(error), error_min = min(error), error = mean(error)) %>%
    ggplot(aes(x = lambda)) +
    geom_pointrange(aes(y = error, ymin = error_min, ymax = error_max))
