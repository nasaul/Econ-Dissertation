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
source("./Functions/penalized_gmm_methods.R")
# Cross Validation

# cv.model <- cv_gmm_alasso(
#   nfolds = 2,
#   known_cond = known_conditions,
#   unknown_cond = unknown_conditions,
#   data        = df,
#   theta_0     = c(0, 0, 0),
#   lambda = seq(0,.25, length.out = 30),
#   eps = 1e-8
# )
# 
# cv.model %>% 
#   group_by(lambda) %>%
#   summarise( error_max = max(error), error_min = min(error), error = mean(error)) %>%
#   ggplot(aes(x = lambda)) +
#   geom_pointrange(aes(y = error, ymin = error_min, ymax = error_max))
#   
# 
# cv.model %>% 
#   select(lambda, moment_values) %>% 
#   mutate(nombres = purrr::map(moment_values, names)) %>% 
#   tidyr::unnest() %>% 
#   group_by(lambda, nombres) %>% 
#   summarise(
#     moment_values_mean = mean(moment_values),
#     moment_values_min = min(moment_values),
#     moment_values_max = max(moment_values)
#   ) %>% 
#   ggplot(aes(x = nombres)) +
#   geom_errorbar(aes(ymin = moment_values_min, ymax = moment_values_max)) +
#   geom_point(aes(y = moment_values_mean)) + 
#   coord_flip() 
# 
# # Alasso GMM
# gmm_alasso(
#   known_cond = known_conditions,
#   unknown_cond = unknown_conditions,
#   data        = df,
#   theta_0     = c(0, 0, 0),
#   lambda = 0.01,
#   eps = 1e-4
# )
# 
# gmm_lasso(
#   known_cond = known_conditions,
#   unknown_cond = unknown_conditions,
#   data        = df,
#   theta_0     = c(0, 0, 0),
#   lambda = 0.01,
#   eps = 1e-4
# )
# 
