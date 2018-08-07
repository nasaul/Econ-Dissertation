# Packages --------------------------------------------------------------------------------------------------------
library(purrr)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
theme_set(theme_minimal())

# Function to obtain selected moments -----------------------------------------------------------------------------
select_moments <- function(x){
  tibble(sim = 1:length(x)) %>% 
    mutate(
      Z21   = map_lgl(x, ~any(grepl("Z21", .$selected_mom))),
      Z22   = map_lgl(x, ~any(grepl("Z22", .$selected_mom))), 
      Z23   = map_lgl(x, ~any(grepl("Z23", .$selected_mom))),
      Z24   = map_lgl(x, ~any(grepl("Z24", .$selected_mom))), 
      Z25   = map_lgl(x, ~any(grepl("Z25", .$selected_mom))),
      mom_f = map_lgl(x, ~any(grepl("F"  , .$selected_mom)))
    ) 
}

# Function to obtain table results --------------------------------------------------------------------------------
prop_table <- function(x){
  x %>% 
    mutate(
      p1 = if_else(mom_f == TRUE, 1, 0),
      p2 = if_else(p1 == 0 & any(Z21, Z22, Z23, Z24, Z25), 1, 0),
      p3 = if_else(p1 == 0 & all(Z21, Z22, Z23, Z24, Z25), 1, 0)
    ) %>% 
    summarise(
      p1 = sum(p1) / nrow(x),
      p2 = sum(p2) / nrow(x),
      p3 = sum(p3) / nrow(x)
    )
}

# Obtain results --------------------------------------------------------------------------------------------------
sim_2_2 <- readRDS(list.files("Results/sim_500_2/",   full.names = TRUE)) %>% 
  select_moments()
sim_5_2 <- map(list.files("Results/sim_500_5/",   full.names = TRUE), readRDS) %>% 
  select_moments()
sim_2_8 <- map(list.files("Results/sim_500_2_8/", full.names = TRUE), readRDS) %>% 
  select_moments()
sim_5_8 <- map(list.files("Results/sim_500_5_8/", full.names = TRUE), readRDS) %>% 
  select_moments()

# Table results ---------------------------------------------------------------------------------------------------
sim_2_2 %>% prop_table()
sim_5_2 %>% prop_table()
sim_2_8 %>% prop_table()
sim_5_8 %>% prop_table()

# Graphic generator -----------------------------------------------------------------------------------------------

dist_info <- function(x){
  x %>% 
    summarise(
      Z21 = sum(Z21) / nrow(x),
      Z22 = sum(Z22) / nrow(x),
      Z23 = sum(Z23) / nrow(x),
      Z24 = sum(Z24) / nrow(x),
      Z25 = sum(Z25) / nrow(x),
      mom_F = sum(mom_f) / nrow(x)
    )
}

dist_plot <- dist_info(sim_2_2) %>% 
  rbind(dist_info(sim_2_8)) %>% 
  rbind(dist_info(sim_5_2)) %>% 
  rbind(dist_info(sim_5_8)) %>% 
  mutate(
    pi = c(0.2, 0.2, 0.8, 0.8),
    cl = c(0.2, 0.5, 0.2, 0.5)
  ) %>% 
  gather(key, value, -pi, -cl) %>% 
  mutate(key = if_else(key == "mom_F", "F", key)) %>% 
  ggplot(aes(x = reorder(key, -value), y = value)) +
  geom_bar(stat = "identity") +
  facet_grid(pi ~ cl, labeller = label_bquote(cols = c[l] == .(cl), rows = pi[1] == .(pi))) +
  labs(
    x = "",
    y = "",
    title = "Distribución de la selección de momentos"
  ) + 
  theme(
    text = element_text(size = 18)
  )

dist_plot
ggsave("Results/plot_dist.pdf", dist_plot, device = "pdf", width = 9)
