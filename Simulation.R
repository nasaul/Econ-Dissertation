library(MASS)
library(tibble)
library(magrittr)
library(dplyr)
library(ggplot2)
library(furrr)
# Simulation for thesis -------------------------------------------------------------------------------------------
source("./Functions/Simulation/simulation_function.R")
source("./Functions/penalized_gmm_functions.R")
set.seed(123456)
plan(multiprocess)
# Number of simulations
nsim <- 18
sim_500_2 <- furrr::future_map(
  1:nsim,
  ~simulation(cl = 0.2, lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-4, N = 500, pi1 = 0.2, sim_num = ., nombre = "sim_500_2"),
  .progress = TRUE
)
sim_500_5  <- furrr::future_map(
  1:nsim,
  ~simulation(cl = 0.5, lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-8, N = 500, pi1 = 0.2, sim_num = ., nombre = "sim_500_5"),
  .progress = TRUE
)
