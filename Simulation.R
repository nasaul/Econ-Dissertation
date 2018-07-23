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
sim_500_2 <- furrr::future_map(
  rep(0.2, 1000),
  ~simulation(cl = ., lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-8, N = 500, pi1 = 0.2)
)

saveRDS(sim_500_2, file = "Results/sim_500_2.rds")
rm(sim_500_2)
gc()
sim_1000_2  <- furrr::future_map(
  rep(0.2, 1000),
  ~simulation(cl = ., lambdas = seq(0.0003030303, 0.01, by = 0.0003030303), eps = 1e-8, N = 1000, pi1 = 0.2)
)

saveRDS(sim_1000_2 , file = "Results/sim_1000_2.rds")