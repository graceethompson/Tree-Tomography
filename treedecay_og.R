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
  
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 0.5), 
        strip.text = element_text(size = 10), 
        axis.title = element_text(size = 11))



# Two-way ANOVA test for SoT percent of solid wood: 
SoT_solid_twowayanova <- aov(percent_solid_wood ~ species * site, 
                             data = treedecayinfo)
summary(SoT_solid_twowayanova)

# Tukey HSD test for two-way ANOVA SoT percent of solid wood: 
TukeyHSD(SoT_solid_twowayanova)

# Adjusting anova for overlapping species
overlapping_spp <- treedecayinfo %>%
  filter(species %in% c("rm", "hem"))

overlapping_solid_twowayanova <- aov(percent_solid_wood ~ species * site, 
                             data = overlapping_spp)
summary(overlapping_solid_twowayanova)

TukeyHSD(overlapping_solid_twowayanova)

par(mfrow = c(2,2))
plot(overlapping_solid_twowayanova)


install.packages("lmPerm")
library(lmPerm)

overlapping_solid_permanova <- aovp(percent_solid_wood ~ species * site, 
                                     data = overlapping_spp)
summary(overlapping_solid_permanova)

TukeyHSD(overlapping_solid_permanova)



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
  
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 0.5), 
        strip.text = element_text(size = 10), 
        axis.title = element_text(size = 11))


# Two-way ANOVA test for SoT percent of damaged wood: 
SoT_damaged_twowayanova <- aov(percent_damaged ~ species * site, 
                               data = treedecayinfo)
summary(SoT_damaged_twowayanova)


overlapping_dam_permanova <- aovp(percent_damaged ~ species * site, 
                                    data = overlapping_spp)
summary(overlapping_dam_permanova)

species_dam_permanova <- aovp(percent_damaged ~ species, 
                                  data = overlapping_spp)

summary(species_dam_permanova)

site_dam_permanova <- aovp(percent_damaged ~ site, 
                              data = overlapping_spp)

summary(site_dam_permanova)


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
    q1 = quantile(percent_solid_wood, 0.25, na.rm = TRUE),
    q3 = quantile(percent_solid_wood, 0.75, na.rm = TRUE),
    iqr = IQR(percent_damaged, na.rm = TRUE)
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
    max = max(percent_damaged, na.rm = TRUE),
    q1 = quantile(percent_damaged, 0.25, na.rm = TRUE),
    q3 = quantile(percent_damaged, 0.75, na.rm = TRUE),
    iqr = IQR(percent_damaged, na.rm = TRUE)
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
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "grey90"),
        axis.title.x = element_text(size = 11), 
        axis.title = element_text(size = 12), 
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 0.5))


# One-way ANOVA test for SoT percent of solid wood:
SoT_solid_onewayanova <- aov(percent_solid_wood ~ site, data = treedecayinfo)
summary(SoT_solid_onewayanova)

TukeyHSD(SoT_solid_onewayanova)

# Welch's t-test 
t.test(percent_solid_wood ~ site,
       data = treedecayinfo)

# Normality by site
ggplot(treedecayinfo, aes(x = percent_solid_wood, fill = site)) +
  geom_histogram(bins = 10) +
  facet_wrap(~site) # Skewed

# Shapiro-Wilk test for each site
by(treedecayinfo$percent_solid_wood,
   treedecayinfo$site,
   shapiro.test) # Assumptions are violated

var.test(percent_solid_wood ~ site,
         data = treedecayinfo) # True ratio of variances is not equal to 1

# QQ plot
ggplot(treedecayinfo, aes(sample = percent_solid_wood)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ site) # Not normal

# Non-parametric alternative (if normality assumption violated)
wilcox.test(percent_solid_wood ~ site, data = treedecayinfo)



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
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "grey90"),
        axis.title.x = element_text(size = 11), 
        axis.title = element_text(size = 12), 
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 0.5))


# One-way ANOVA test for SoT percent of damaged wood:
SoT_damaged_onewayanova <- aov(percent_damaged ~ site, data = treedecayinfo)
summary(SoT_damaged_onewayanova)

TukeyHSD(SoT_damaged_onewayanova)

# Welch's t-test
t.test(percent_damaged ~ site,
       data = treedecayinfo) 

# Normality by site
ggplot(treedecayinfo, aes(x = percent_damaged, fill = site)) +
  geom_histogram(bins = 10) +
  facet_wrap(~site) # Skewed

# Shapiro-Wilk test for each site
by(treedecayinfo$percent_damaged,
   treedecayinfo$site,
   shapiro.test) # Assumptions are violated

var.test(percent_damaged ~ site,
         data = treedecayinfo) # True ratio of variances is not equal to 1

# QQ plot
ggplot(treedecayinfo, aes(sample = percent_damaged)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ site) # Not normal

# Non-parametric alternative (if normality assumption violated)
wilcox.test(percent_damaged ~ site, data = treedecayinfo)


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



# Trying a beta regression generalized linear model 
install.packages("betareg")
library(betareg)


treedecayinfo <- treedecayinfo %>%
  mutate(damaged_prop = percent_damaged / 100)

summary(treedecayinfo$damaged_prop)
range(treedecayinfo$damaged_prop, na.rm = TRUE)
sum(is.na(treedecayinfo$damaged_prop))

sum(treedecayinfo$damaged_prop == 0.001)
nrow(treedecayinfo)
table(treedecayinfo$species, treedecayinfo$site)

n <- nrow(treedecayinfo)

treedecayinfo$damaged_prop_adj <- 
  (treedecayinfo$damaged_prop * (n - 1) + 0.5) / n


overlapping_betareg <- betareg(
  damaged_prop ~ species * site, data = overlapping_spp)
summary(overlapping_betareg)








