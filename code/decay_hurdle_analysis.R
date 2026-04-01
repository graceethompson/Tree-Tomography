# =============================================================================
# Tree Decay Analysis: Effects of Species and Site on Wood Decay
# =============================================================================
# This script uses a two-part (hurdle) modeling approach because our response
# variable has two distinct biological processes:
#   1. WHETHER decay occurs at all (presence/absence)
#   2. HOW MUCH decay occurs, given that it is present
#
# ~74% of trees show zero decay — standard regression assumes a continuous response, and
# zero-inflated data like this violates that assumption.
# =============================================================================

library(tidyverse)
library(lme4)  
install.packages("betareg")
library(betareg)
library(dplyr)  
library(ggplot2)
install.packages("emmeans")
library(emmeans) 
install.packages("DHARMa")
library(DHARMa)
#install.packages("patchwork")
#library(patchwork)


treedecayinfo <- read.csv("Tree_ID_info.csv")
treedecayinfo$site <- as.factor(ifelse(treedecayinfo$plot == "BGS", "BGS", "EMS"))
treedecayinfo$species <- as.factor(treedecayinfo$species)

# =============================================================================
# PART 1: Does decay occur? (Binary presence/absence)
# =============================================================================
# Generalized Linear Model (GLM) with a binomial family, which is
# the appropriate when your response variable is binary (0 or 1). 
# We include both species and site as predictors because we have an a priori
# expectation that both could influence decay.
# =============================================================================
treedecayinfo$decay_binary <- ifelse(treedecayinfo$percent_damaged > 0, 1, 0)

m1_binary <- glm(decay_binary ~ species + site,
                 family = binomial,
                 data = treedecayinfo)

summary(m1_binary)

# fancy residual tests before we go too far
m1_simres <- simulateResiduals(m1_binary)
plot(m1_simres) # looks good! Ignore the little bit of red, this isn't bad for DHARMa

# --- Model comparison ----------------------------
# To formally test whether species and site each contribute to explaining
# decay presence, we build a series of simpler models and compare them.
# The anova() with test = "Chisq" asks "does adding this predictor significantly improve model fit?"
# sig p-value = meaningful variation.

# Null model: no predictors baseline for comparison
m_null <- glm(decay_binary ~ 1,
              family = binomial,
              data = treedecayinfo)

# Site-only model
m_site <- glm(decay_binary ~ site,
              family = binomial,
              data = treedecayinfo)

# Species-only model
m_species <- glm(decay_binary ~ species,
                 family = binomial,
                 data = treedecayinfo)

# Full model with both predictors
m_full <- glm(decay_binary ~ species + site,
              family = binomial,
              data = treedecayinfo) # Same as m1_binary

# Test 1: Does site improve the model, and does species add anything beyond site?
# each row tests whether adding the next predictor significantly reduces residual deviance (unexplained variation).
anova(m_null, m_site, m_full, test = "Chisq") # YES, site contributes

# Test 2: Does species improve the model, and does site add anything beyond species?
anova(m_null, m_species, m_full, test = "Chisq") # NO, species does not (full does better)


# =============================================================================
# PART 2: How much decay, given that decay is present? (Continuous severity)
# =============================================================================
# beta regression, which is designed for response variables that are
# proportions strictly between 0 and 1. 
# We subset to only trees where decay is present (percent_damaged > 0),
# because the zeros represent a separate biological process already modeled
# in Part 1. Including zeros in a beta regression would also violate its assumptions.

treedecayinfo$percent_damaged <- treedecayinfo$percent_damaged / 100

# =============================================================================

# Fit species-only model
m2_species_only <- betareg(percent_damaged ~ species,
                           data = filter(treedecayinfo, percent_damaged > 0))

# Fit model with species + site to test whether site matters for severity
m2_beta_full <- betareg(percent_damaged ~ species + site,
                        data = filter(treedecayinfo, percent_damaged > 0))

# Compare AIC between the two models, lower AIC = better fit
# A difference of >2 AIC points is generally considered meaningful
AIC(m2_species_only, m2_beta_full) # close enough to two for me, species-only is better

# This is our final severity model:
m2_beta <- betareg(percent_damaged ~ species,
                   data = filter(treedecayinfo, percent_damaged > 0))


summary(m2_beta)

# =============================================================================
# Now for the fun part. pictures.
#   A) Model-predicted decay presence probability by species (from m1_binary)
#   B) Raw decay prevalence by site (to illustrate the site effect from m1)
#   C) Model-predicted decay severity by species (from m2_beta)
# =============================================================================

# --- Panel A: Decay presence probability by species --------------------------

p1_data <- emmeans(m1_binary, ~ species, type = "response") %>%
  as.data.frame()


p1 <- ggplot(p1_data, aes(x = reorder(species, prob), y = prob)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Species",
       y = "Decay presence probability",
       title = "A) Decay presence") +
  theme_classic() # because classy graphs look better

ggplot(p1_data, aes(x = reorder(species, prob), y = prob)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(x = "Species",
       y = "Decay presence probability",
       title = "A) Decay presence") +
  theme_classic()


# --- Panel B: Raw decay prevalence by site -----------------------------------
# descriptive summary — the observed proportion of trees with decay at each site, with standard error bars.

p2_data <- treedecayinfo %>%
  group_by(site) %>%
  summarise(prop = mean(decay_binary),
            se = sqrt(prop * (1 - prop) / n()))

p2 <- ggplot(p2_data, aes(x = site, y = prop)) +
  geom_col(fill = "grey70", width = 0.5) +
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.15) +
  labs(x = "Site",
       y = "Proportion with decay",
       title = "B) Decay prevalence by site") +
  theme_classic()

ggplot(p2_data, aes(x = site, y = prop)) +
  geom_col(fill = "grey70", width = 0.5) +
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.15) +
  labs(x = "Site",
       y = "Proportion with decay",
       title = "B) Decay prevalence by site") +
  theme_classic()

# --- Panel C: Decay severity by species --------------------------------------
# manually back-transform using plogis() (the inverse logit function) to get values on the proportion scale (0-1).

p3_data <- emmeans(m2_beta, ~ species, type = "response") %>%
  as.data.frame()

p3_data$response <- plogis(p3_data$emmean)
p3_data$LCL      <- plogis(p3_data$asymp.LCL)
p3_data$UCL      <- plogis(p3_data$asymp.UCL)


p3 <- ggplot(p3_data, aes(x = reorder(species, response), y = response)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2) +
  labs(x = "Species",
       y = "Decay severity (proportion damaged)",
       title = "C) Decay severity | decay present") +
  theme_classic()

ggplot(p3_data, aes(x = reorder(species, response), y = response)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2) +
  labs(x = "Species",
       y = "Decay severity (proportion damaged)",
       title = "C) Decay severity | decay present") +
  theme_classic()

# --- Combine panels ----------------------------------------------------------
# combine with patchwork because they look pretty this way

p1 + p2 + p3 # Patchwork does not work in my version of R 



