library(MASS)
library(tibble)
library(magrittr)
library(dplyr)
library(ggplot2)
library(furrr)
# Simulation for thesis -------------------------------------------------------------------------------------------
source("./Functions/Simulation/simulation_function.R")
source("./Functions/penalized_gmm_functions.R")
plan(multiprocess)
nsim <- 1000
set.seed(123456)
sim_500_2 <- furrr::future_map(
  1:nsim,
  ~simulation(cl = 0.2, lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-8, N = 500, pi1 = 0.2, sim_num = ., nombre = "sim_500_2"),
  .progress = TRUE
)
set.seed(123456)
sim_500_5  <- furrr::future_map(
  1:nsim,
  ~simulation(cl = 0.5, lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-8, N = 500, pi1 = 0.2, sim_num = ., nombre = "sim_500_5"),
  .progress = TRUE
)
set.seed(123456)
sim_500_2 <- furrr::future_map(
  1:nsim,
  ~simulation(cl = 0.2, lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-8, N = 500, pi1 = 0.8, sim_num = ., nombre = "sim_500_2_8"),
  .progress = TRUE
)
set.seed(123456)
sim_500_5  <- furrr::future_map(
  1:nsim,
  ~simulation(cl = 0.5, lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-8, N = 500, pi1 = 0.8, sim_num = ., nombre = "sim_500_5_8"),
  .progress = TRUE
)