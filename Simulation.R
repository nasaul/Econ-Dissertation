library(MASS)
library(tibble)
library(magrittr)
library(dplyr)
library(ggplot2)
library(furrr)

# Simulate Endogeneity ----------------------------------------------------
source("./Functions/Simulation/simulation_function.R")
source("./Functions/penalized_gmm_functions.R")
set.seed(123456)
plan(multiprocess)
sim_1000 <- furrr::future_map(
  rep(0.2, 1000),
  ~simulation(cl = ., lambdas = seq(0, 0.2, length.out = 8), eps = 1e-8, N = 500, pi1 = 0.2)
)
