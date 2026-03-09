library(tidyverse)
library(dplyr)
library(ggplot2)
library(broom)

# Read in tree decay dataset with plot, species, tree tag, dbh, crack detected, SoT percent solid wood, and SoT percent damaged, 
read_csv("Tree_ID_info.csv")
treedecayinfo <- read_csv("Tree_ID_info.csv")

# Create table with BGS and EMS distinction
treedecayinfo <- treedecayinfo %>%
  mutate(site = case_when(
    plot %in% c("C1", "C2", "D1", "D2", "E1", "E5", "F1", "G1", "H1") ~ "EMS",
    TRUE ~ "BGS"))

# Don't know why I have this here 
treedecayinfo <- treedecayinfo %>%
  group_by(site) %>% # From this point forward, any grouped operations are performed within each site
  mutate(species = factor(species)) %>% # Converts species into a factor
  ungroup() # Removes the grouping structure so future operations are not still grouped by site

# Reorganizing the order of names for box plots
treedecayinfo$species <- factor(treedecayinfo$species, 
                                    levels = c("hem", "rm", "bg", "ro"))



# Create box plots comparing species variation in SoT percent of solid wood and damaged wood between sites

library(scales)  # for alpha()

# For SoT percent of solid wood:
ggplot(treedecayinfo) + 
  geom_boxplot(
    aes(species, percent_solid_wood, color = species, fill = species),
    alpha = 0.2, # transparency for the fill
    outlier.shape = NA
    ) + 
  scale_y_continuous(limits = c(0, 101)) +
  geom_jitter(
    aes(species, percent_solid_wood, color = species),
    width = 0.2, alpha = 0.5
  ) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Solid Wood (%)") +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) + 
  
  scale_color_manual(
    values = c(
      hem = "palegreen4",
      rm  = "darkorange2",
      bg  = "steelblue3",
      ro  = "orchid3"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      hem = alpha("palegreen4", 0.2),
      rm  = alpha("darkorange2", 0.2),
      bg  = alpha("steelblue3", 0.2),
      ro  = alpha("orchid3", 0.2)
    )
  ) +
  
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(), 
        strip.text = element_text(size = 11), 
        axis.title = element_text(size = 12))


# Two-way ANOVA test for SoT percent of solid wood: 
SoT_solid_twowayanova <- aov(percent_solid_wood ~ species * site, 
                             data = treedecayinfo)
summary(SoT_solid_twowayanova)

# Tukey HSD test for two-way ANOVA SoT percent of solid wood: 
TukeyHSD(SoT_solid_twowayanova)



# For SoT percent of damaged wood:
ggplot(treedecayinfo) + 
  geom_boxplot(
    aes(species, percent_damaged, color = species, fill = species),
    alpha = 0.2,  # transparency for the fill
    outlier.shape = NA
  ) + 
  scale_y_continuous(limits = c(-1, 100)) +
  geom_jitter(
    aes(species, percent_damaged, color = species),
    width = 0.2, alpha = 0.5
  ) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Damaged Wood (%)") +
  xlab("Species") +
  
  scale_x_discrete(labels = c(
    hem = expression(italic("T. canadensis")),
    rm  = expression(italic("A. rubrum")),
    bg  = expression(italic("N. sylvatica")),
    ro  = expression(italic("Q. rubra"))
  )) +
  
  scale_color_manual(
    values = c(
      hem = "palegreen4",
      rm  = "darkorange2",
      bg  = "steelblue3",
      ro  = "orchid3"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      hem = alpha("palegreen4", 0.2),
      rm  = alpha("darkorange2", 0.2),
      bg  = alpha("steelblue3", 0.2),
      ro  = alpha("orchid3", 0.2)
    )
  ) +
  
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(), 
        strip.text = element_text(size = 11), 
        axis.title = element_text(size = 12))


# Two-way ANOVA test for SoT percent of damaged wood: 
SoT_damaged_twowayanova <- aov(percent_damaged ~ species * site, 
                               data = treedecayinfo)
summary(SoT_damaged_twowayanova)


# Basic statistics for SoT percent of solid wood 
treedecayinfo %>%
  group_by(site) %>%
  summarise(
    n = sum(!is.na(percent_solid_wood)),
    mean = mean(percent_solid_wood, na.rm = TRUE),
    sd = sd(percent_solid_wood, na.rm = TRUE),
    se = sd / sqrt(n),
    median = median(percent_solid_wood, na.rm = TRUE),
    min = min(percent_solid_wood, na.rm = TRUE),
    max = max(percent_solid_wood, na.rm = TRUE),
    .groups = "drop"
  )

# Basic statistics for SoT percent of damaged wood 
treedecayinfo %>%
  group_by(site) %>%
  summarise(
    n = sum(!is.na(percent_damaged)),
    mean = mean(percent_damaged, na.rm = TRUE),
    sd = sd(percent_damaged, na.rm = TRUE),
    se = sd / sqrt(n),
    median = median(percent_damaged, na.rm = TRUE),
    min = min(percent_damaged, na.rm = TRUE),
    max = max(percent_damaged, na.rm = TRUE)
  )


# Between sites comparison
# For SoT percent of solid wood: 
ggplot(treedecayinfo) +
  geom_boxplot(aes(
    x = site, y = percent_solid_wood, color = site, fill = site), 
    linewidth = 0.6,
    alpha = 0.2,
    outlier.shape = NA) +
  geom_jitter(aes(site, percent_solid_wood, color = site), 
              width = 0.2, alpha = 0.5) + 
  ylab("Solid Wood (%)") +
  xlab("Site") +
  scale_y_continuous(limits = c(0, 101)) +
  scale_color_manual(
    values = c(
      BGS = "cyan4",
      EMS  = "gold3"
    )
  ) +
  scale_fill_manual(
    values = c(
      BGS = alpha("cyan4", 0.2),
      EMS  = alpha("gold3", 0.2)
    )
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 11), 
        axis.title = element_text(size = 12))


# One-way ANOVA test for SoT percent of solid wood:
SoT_solid_onewayanova <- aov(percent_solid_wood ~ site, data = treedecayinfo)
summary(SoT_solid_onewayanova)

TukeyHSD(SoT_solid_onewayanova)


# For SoT percent of damaged wood:
ggplot(treedecayinfo) +
  geom_boxplot(aes(
    x = site, y = percent_damaged, color = site, fill = site),
    linewidth = 0.6,
    alpha = 0.2,
    outlier.shape = NA) +
  geom_jitter(aes(site, percent_damaged, color = site),
              width = 0.2, alpha = 0.5) + 
  ylab("Damaged Wood (%)") +
  xlab("Site") +
  scale_y_continuous(limits = c(-1, 100)) +
  scale_color_manual(
    values = c(
      BGS = "cyan4",
      EMS  = "gold3"
    )
  ) +
  scale_fill_manual(
    values = c(
      BGS = alpha("cyan4", 0.2),
      EMS  = alpha("gold3", 0.2)
    )
  ) +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 11), 
        axis.title = element_text(size = 12))


# One-way ANOVA test for SoT percent of damaged wood:
SoT_damaged_onewayanova <- aov(percent_damaged ~ site, data = treedecayinfo)
summary(SoT_damaged_onewayanova)

TukeyHSD(SoT_damaged_onewayanova)


# Alternative data visualization to the first boxplots
summary_damaged <- treedecayinfo %>%
  group_by(species, site) %>%
  summarise(
    n = n(),
    mean_damaged = mean(percent_damaged, na.rm = TRUE),
    sd_damaged = sd(percent_damaged, na.rm = TRUE),
    se_damaged = sd_damaged / sqrt(n)
  ) %>%
  ungroup()

ggplot(treedecayinfo, aes(x = species, y = percent_damaged, color = site)) +
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.2, dodge.width = 0.8),alpha = 0.6) +
  geom_point(data = summary_damaged, 
             aes(x = species, y = mean_damaged),
             position = position_dodge(width = 0.8),
             color = "black", size = 3) +
  geom_errorbar(data = summary_damaged,
                aes(x = species, y = mean_damaged,
                    ymin = mean_damaged - se_damaged,
                    ymax = mean_damaged + se_damaged),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  theme_classic() +
  labs(
    x = "Species",
    y = "Damaged Wood (%)",
    color = "Site"
  )


# Visual exploration of adding tree tags to points
ggplot(treedecayinfo, 
       aes(species, percent_solid_wood, color = species)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 1) +
  geom_text(aes(label = tree), vjust = -0.5, color = "black", 
            position = position_dodge(width = 0.8)) +
  facet_wrap(~ site, scales = "free_x") +
  ylab("Solid Wood (%)") +
  xlab("Species") +
  scale_y_continuous(limits = c(40, 105)) +
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  scale_x_discrete(drop = TRUE) +
  theme_minimal()

 
ggplot(treedecayinfo, aes(species, percent_damaged, color = species)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 1) +
  geom_text(aes(label = tree), vjust = -0.5, color = "black", 
            position = position_dodge(width = 0.8)) +
  facet_wrap(~ site, scales = "free_x") + 
  ylab("Damaged Wood (%)") + 
  xlab("Species") + 
  scale_y_continuous(limits = c(-5, 100)) +
  scale_color_manual(values = c(
    hem = "palegreen4",
    rm  = "darkorange2",
    bg  = "steelblue3",
    ro  = "orchid3"
  )) +
  scale_x_discrete(drop = TRUE) + 
  theme_minimal() 






