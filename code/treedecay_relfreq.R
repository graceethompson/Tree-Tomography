library(tidyverse)
library(dplyr)
library(ggplot2)
library(broom)

# Use presence and absence to find the relative frequency of species per plot with rot 
rot_summary <- BGS_EMS_treedecay %>%
  group_by(plot, species, rot) %>%
  summarise(num_trees = n(), .groups = "drop")

# Ensure 'species' is a factor with levels only present in the data
rot_summary <- rot_summary %>%
  group_by(plot) %>%
  mutate(species = factor(species)) %>%
  ungroup()

BGS_EMS_treedecay <- BGS_EMS_treedecay %>%
  mutate(damaged_prop = percent_damaged / 100)


















