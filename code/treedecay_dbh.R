# DBH Analysis 
# I want to model a linear relationship between DBH and percent decay
# I want to differentiate between species 
# Facet wrap by site 
# Facet wrap by species 

library(tidyverse)
library(dplyr)
library(ggplot2)

read_csv("Tree_ID_info.csv")
treedecayinfo <- read_csv("Tree_ID_info.csv")

# Create table with BGS and EMS distinction
treedecayinfo <- treedecayinfo %>%
  mutate(site = case_when(
    plot %in% c("C1", "C2", "D1", "D2", "E1", "E5", "F1", "G1", "H1") ~ "EMS",
    TRUE ~ "BGS"))

# Linear regression model
ggplot(treedecayinfo, 
       aes(x = dbh, y = percent_damaged, color = species)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~site) + 
  xlab("DBH (cm)") +
  ylab("Damaged Wood (%)") +
  theme_bw()

ggplot(treedecayinfo, 
       aes(x = dbh, y = percent_damaged, color = site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~species) + 
  xlab("DBH (cm)") +
  ylab("Damaged Wood (%)") +
  theme_bw()

ggplot(treedecayinfo, 
       aes(x = dbh, y = percent_damaged, color = species, shape = site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("DBH (cm)") +
  ylab("Damaged Wood (%)") +
  theme_bw()

ggplot(treedecayinfo, 
       aes(x = dbh, y = percent_damaged, color = site, shape = site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("DBH (cm)") +
  ylab("Damaged Wood (%)") +
  theme_bw()

ggplot(treedecayinfo, 
       aes(x = dbh, y = percent_damaged, color = species)) +
  geom_point(size = 2) +
  geom_smooth(aes(linetype = site), method = "lm", se = FALSE) +
  xlab("DBH (cm)") +
  ylab("Damaged Wood (%)") +
  theme_bw()


dbh_onewayanova <- aov(dbh ~ species, data = treedecayinfo)
summary(dbh_onewayanova)
TukeyHSD(dbh_onewayanova)

BGS_dbh_onewayanova <- aov(dbh ~ species, data = BGS_treedecayinfo)
summary(BGS_dbh_onewayanova)
TukeyHSD(BGS_dbh_onewayanova)

EMS_dbh_onewayanova <- aov(dbh ~ species, data = EMS_treedecayinfo)
summary(EMS_dbh_onewayanova)
TukeyHSD(EMS_dbh_onewayanova)

ggplot(treedecayinfo) +
  geom_boxplot(aes(
    x = species, y = dbh, color = species, fill = species), 
    linewidth = 0.6,
    alpha = 0.2,
    outlier.shape = NA) +
  geom_jitter(aes(species, dbh, color = species), 
              width = 0.2, alpha = 0.5) + 
  ylab("DBH (cm)") +
  xlab("Species") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 11), 
        axis.title = element_text(size = 12))

treedecayinfo %>%
  group_by(species) %>%
  summarise(
    n        = sum(!is.na(dbh)),
    avg      = mean(dbh, na.rm = TRUE),
    std_dev  = sd(dbh, na.rm = TRUE),
    std_err  = std_dev / sqrt(n),
    med      = median(dbh, na.rm = TRUE),
    min_val  = min(dbh, na.rm = TRUE),
    max_val  = max(dbh, na.rm = TRUE)
  )

#------------------------------------------------------------------------------

#library(tidyverse)
library(lme4)  
install.packages("betareg")
library(betareg)
#library(dplyr)  
#library(ggplot2)
install.packages("emmeans")
library(emmeans) 
install.packages("DHARMa")
library(DHARMa)
#install.packages("patchwork")
#library(patchwork)


#treedecayinfo <- read.csv("Tree_ID_info.csv")
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

m1_binary_dbh <- glm(decay_binary ~ species + site + dbh,
                 family = binomial,
                 data = treedecayinfo)

summary(m1_binary_dbh)

# fancy residual tests before we go too far
m1_simres_dbh <- simulateResiduals(m1_binary_dbh)
plot(m1_simres_dbh) # looks good! Ignore the little bit of red, this isn't bad for DHARMa

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

# DBH-only model
m_dbh <- glm(decay_binary ~ dbh,
                 family = binomial,
                 data = treedecayinfo)

# Full model with both predictors
m_full_dbh <- glm(decay_binary ~ species + site + dbh,
              family = binomial,
              data = treedecayinfo) # Same as m1_binary_dbh

# Test 1: Does site improve the model, and does species add anything beyond site?
# each row tests whether adding the next predictor significantly reduces residual deviance (unexplained variation).
anova(m_null, m_site, m_full_dbh, test = "Chisq") # YES, site contributes

# Test 2: Does species improve the model, and does site add anything beyond species?
anova(m_null, m_species, m_full_dbh, test = "Chisq") # NO, species does not (full does better)

# Test 3: Does dbh improve the model, and does site add anything beyond species?
anova(m_null, m_dbh, m_full_dbh, test = "Chisq") # DBH slightly improved model 


m_site_dbh <- glm(decay_binary ~ site + dbh,
                  family = binomial,
                  data = treedecayinfo)

anova(m_null, m_site, m_site_dbh, test = "Chisq")
AIC(m_null, m_species, m_site, m_site_dbh, m_full_dbh)

summary(m_site)
summary(m_site_dbh)

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

# Fit site-only model
m2_site_only <- betareg(percent_damaged ~ site,
                           data = filter(treedecayinfo, percent_damaged > 0))

# Fit dbh-only model
m2_dbh_only <- betareg(percent_damaged ~ dbh,
                           data = filter(treedecayinfo, percent_damaged > 0))


# Fit model with species + site to test whether site matters for severity
m2_beta_spp_site <- betareg(percent_damaged ~ species + site,
                        data = filter(treedecayinfo, percent_damaged > 0))

# Fit model with species + dbh to test whether dbh matters for severity
m2_beta_spp_dbh <- betareg(percent_damaged ~ species + dbh,
                            data = filter(treedecayinfo, percent_damaged > 0))

# Fit model with dbh + site to test whether site matters for severity
m2_beta_site_dbh <- betareg(percent_damaged ~ dbh + site,
                            data = filter(treedecayinfo, percent_damaged > 0))

# Fit model with species + site to test whether site matters for severity
m2_beta_all <- betareg(percent_damaged ~ species + site + dbh,
                            data = filter(treedecayinfo, percent_damaged > 0))


# Compare AIC between the two models, lower AIC = better fit
# A difference of >2 AIC points is generally considered meaningful
AIC(m2_species_only, m2_beta_all) # close enough to two for me, species-only is better

AIC(m2_species_only, m2_beta_spp_dbh) # species-only is better

AIC(m2_species_only, m2_beta_spp_site) # species-only is better

AIC(m2_site_only, m2_beta_all) # site-only is better

AIC(m2_site_only, m2_beta_spp_site) # no difference

AIC(m2_site_only, m2_beta_site_dbh) # dbh-only is better, but not > 2

AIC(m2_dbh_only, m2_beta_all) # dbh only does better

AIC(m2_dbh_only, m2_beta_spp_dbh) # no difference

AIC(m2_dbh_only, m2_beta_site_dbh) # dbh only does 


# You never directly compared all single-predictor models together
AIC(m2_species_only, m2_site_only, m2_dbh_only)

# You never tested species + site + dbh against species-only directly
# (you only compared full vs. species and full vs. site separately)
AIC(m2_species_only, m2_beta_spp_site, m2_beta_spp_dbh, m2_beta_all)

# You never tested species + site head-to-head against species-only clearly
AIC(m2_species_only, m2_beta_spp_site)

# Missing: site + dbh vs species-only
AIC(m2_species_only, m2_beta_site_dbh)


#m2_null <- betareg(percent_damaged ~ 1,
                   #data = filter(treedecayinfo, percent_damaged > 0))

#AIC(m2_null, m2_species_only) # adding species did not improve model

summary(m2_species_only)

summary(m2_beta_spp_dbh)


