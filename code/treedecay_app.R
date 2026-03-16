library(tidyverse)
library(dplyr)
library(ggplot2)
library(broom)

# Start of using ERT data 
### BASIC STATISTICS ###
# Mean: Average resistivity value within the tree boundary. Higher values indicate drier conditions.
# Median: Middle value when all pixels are sorted. Less sensitive to outliers than mean.
# SD (St Dev): Measure of variability in resistivity values. Higher SD indicates more heterogeneous moisture distribution.

### HETEROGENEITY METRICS ###
# CV (Coefficient of Variation): SD divided by mean. Normalized measure of relative variability (0 = uniform, higher = more variable). Values typically range from 0.1 to 1.0.
# Gini Coefficient: Measure of inequality in the distribution (0 = perfect equality, 1 = maximum inequality). Higher values indicate concentrated moisture in specific areas.
# Entropy: Information-theoretic measure of disorder. Higher entropy indicates more random/dispersed moisture patterns.

### SPATIAL METRICS ###
# CMA (Central Moisture Accumulation): Proportion of the wettest pixels (lowest 30% resistivity) that are located in the central third of the tree. Values > 0.33 indicate moisture concentration in the center.
# Radial Gradient: Difference between edge mean and center mean resistivity. Positive values indicate drier edges (healthy pattern), negative values indicate drier center (possible heartwood or decay).
# Radial Profile: Mean resistivity values in 8 concentric rings from center to edge. Visualizes the radial moisture distribution pattern.


read.csv("data/ERT_application_results.csv")
absolute_ERTs <- read.csv("data/ERT_application_results.csv")

# Combining ERT and SoT data
combined_SoT_ERT <- left_join(treedecayinfo, absolute_ERTs, by = "tree")
EMS_ERTs_SoTs <- subset(combined_SoT_ERT, plot == "EMS")
BGS_ERTs_SoTs <- subset(combined_SoT_ERT, plot == "BGS")


# Adjusting order
combined_SoT_ERT$species <- factor(combined_SoT_ERT$species, 
                                   levels = c("hem", "rm", "bg", "ro"))


## MEAN RESISTIVITY
ggplot(combined_SoT_ERT, 
       aes(species, mean, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab(expression("Mean Resistivity ("*Omega %.% m*")")) +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_mean_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for mean resistivity
ERT_mean_twowayanova <- aov(mean ~ species * site, data = combined_SoT_ERT)
summary(ERT_mean_twowayanova)


# Tukey HSD test for two-way ANOVA for mean resistivity
TukeyHSD(ERT_mean_twowayanova)



# GINI COEFFICIENT
ggplot(combined_SoT_ERT, 
       aes(species, gini, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Gini Coefficient") +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +

  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_gini_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for Gini coefficient 
ERT_gini_twowayanova <- aov(gini ~ species * site, data = combined_SoT_ERT)
summary(ERT_gini_twowayanova)


# Tukey HSD test for two-way ANOVA for Gini coefficient 
TukeyHSD(ERT_gini_twowayanova)



## CENTRAL MOISTURE ACCUMULATION (CMA)
ggplot(combined_SoT_ERT, 
       aes(species, cma, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Central Moisture Accumulation") +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_cma_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for CMA
ERT_cma_twowayanova <- aov(cma ~ species * site, data = combined_SoT_ERT)
summary(ERT_cma_twowayanova)


# Tukey HSD test for two-way ANOVA for CMA
TukeyHSD(ERT_cma_twowayanova)



## RADIAL GRADIENT
ggplot(combined_SoT_ERT, 
       aes(species, radialgradient, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) + 
  scale_y_continuous(limits = c(-601, 600)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Radial Gradient") +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_radialgradient_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for radial gradient
ERT_radialgradient_twowayanova <- aov(radialgradient ~ species * site, data = combined_SoT_ERT)
summary(ERT_radialgradient_twowayanova)


# Tukey HSD test for two-way ANOVA for CMA
TukeyHSD(ERT_radialgradient_twowayanova)



## MEDIAN RESISTIVITY
ggplot(combined_SoT_ERT, 
       aes(species, median, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab(expression("Median Resistivity ("*Omega %.% m*")")) +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_median_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for median resistivity
ERT_median_twowayanova <- aov(median ~ species * site, data = combined_SoT_ERT)
summary(ERT_median_twowayanova)


# Tukey HSD test for two-way ANOVA for mean resistivity
TukeyHSD(ERT_median_twowayanova)


## STANDARD DEVIATION
ggplot(combined_SoT_ERT, 
       aes(species, sd, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab(expression("Standard Deviation ("*Omega %.% m*")")) +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_sd_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for standard deviation 
ERT_sd_twowayanova <- aov(sd ~ species * site, data = combined_SoT_ERT)
summary(ERT_sd_twowayanova)


# Tukey HSD test for two-way ANOVA for standard deviation
TukeyHSD(ERT_sd_twowayanova)


## COEFFICIENT OF VARIATION
ggplot(combined_SoT_ERT, 
       aes(species, cv, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Coefficient of Variation") +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_cv_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for CV
ERT_cv_twowayanova <- aov(cv ~ species * site, data = combined_SoT_ERT)
summary(ERT_cv_twowayanova)


# Tukey HSD test for two-way ANOVA for CV
TukeyHSD(ERT_cv_twowayanova)



## SHANNON ENTROPY
ggplot(combined_SoT_ERT, 
       aes(species, entropy, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.2, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Shannon Entropy") +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank())

ggsave("output/figures/ert_entropy_by_species_site.pdf", width = 8, height = 5)

# Two-way ANOVA for Shannon entropy 
ERT_entropy_twowayanova <- aov(entropy ~ species * site, data = combined_SoT_ERT)
summary(ERT_entropy_twowayanova)


# Tukey HSD test for two-way ANOVA for Shannon entropy
TukeyHSD(ERT_entropy_twowayanova)



# Try to facet wrap all boxplots
library(tidyr)

simplified_ERT <- combined_SoT_ERT %>%
  select(plot, species, tree, dbh, site, cma, cv, mean, radialgradient) 


long_combined_SoT_ERT <- simplified_ERT %>%
  pivot_longer(
    cols = c(mean, cma, radialgradient, cv),
    names_to = "metric",
    values_to = "value"
  )

long_combined_SoT_ERT$species <- droplevels(long_combined_SoT_ERT$species)

ggplot(long_combined_SoT_ERT,
       aes(species, value, color = species, fill = species)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  
  facet_grid(metric ~ site, scales = "free", space = "free_x",
             labeller = labeller(metric = c(
               cma = "CMA",
               cv = "CV",
               mean = "Mean~(Omega %.% m)",
               radialgradient = "Radial~Gradient"
             ), .default = label_parsed)) +
  
  scale_x_discrete(drop = TRUE) +
  ylab(NULL) + 
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  
  scale_fill_manual(values = c(
    hem = scales::alpha("palegreen4", 0.2),
    rm  = scales::alpha("darkorange2", 0.2),
    bg  = scales::alpha("steelblue3", 0.2),
    ro  = scales::alpha("orchid3", 0.2)
  )) +
  
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
    strip.text.x = element_text(face = "bold"),  
    strip.text.y = element_text(face = "bold")
  )

ggsave("output/figures/ert_combined_metrics_by_species_site.pdf", width = 10, height = 10)


