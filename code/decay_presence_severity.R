# =============================================================================
# Decay presence and severity analysis
# Binomial GLM for decay presence, beta regression for severity
# Produces 3-panel figure: (A) presence by species, (B) prevalence by site,
#                          (C) severity given decay present by species
# =============================================================================

library(tidyverse)
library(betareg)
library(patchwork)

# --- 1. DATA PREP -----------------------------------------------------------

df <- read_csv("data/Tree_ID_info.csv", show_col_types = FALSE) %>%
  mutate(
    site = if_else(plot %in% c("C1","C2","D1","D2","E1","E5","F1","G1","H1"),
                   "EMS", "BGS"),
    decay_binary = as.integer(decay == "present"),
    species = factor(species, levels = c("rm", "bg", "ro", "hem"))
  )

cat("Sample sizes by species × site:\n")
print(table(df$species, df$site))
cat("\nDecay counts:\n")
print(table(df$species, df$decay))

# --- 2. PANEL A: Decay presence ~ species (binomial GLM) --------------------

glm_presence <- glm(decay_binary ~ species + site, data = df, family = binomial)
cat("\n=== Binomial GLM: decay presence ~ species + site ===\n")
print(summary(glm_presence))

# Predicted probabilities per species (marginalised over site)
newdata_spp <- expand.grid(
  species = levels(df$species),
  site = c("BGS", "EMS")
)
newdata_spp$pred <- predict(glm_presence, newdata = newdata_spp, type = "response")

# Average over sites for marginal species effect
pred_spp <- newdata_spp %>%
  group_by(species) %>%
  summarise(prob = mean(pred), .groups = "drop")

# Bootstrap CIs for species marginal probabilities
set.seed(42)
n_boot <- 2000
boot_mat <- matrix(NA, nrow = n_boot, ncol = nlevels(df$species))
colnames(boot_mat) <- levels(df$species)

for (i in seq_len(n_boot)) {
  idx <- sample(nrow(df), replace = TRUE)
  boot_df <- df[idx, ]
  fit <- tryCatch(
    glm(decay_binary ~ species + site, data = boot_df, family = binomial),
    warning = function(w) suppressWarnings(
      glm(decay_binary ~ species + site, data = boot_df, family = binomial)
    )
  )
  preds <- predict(fit, newdata = newdata_spp, type = "response")
  boot_preds <- data.frame(species = newdata_spp$species, pred = preds)
  marg <- tapply(boot_preds$pred, boot_preds$species, mean)
  boot_mat[i, ] <- marg[levels(df$species)]
}

pred_spp$lo <- apply(boot_mat, 2, quantile, 0.025)
pred_spp$hi <- apply(boot_mat, 2, quantile, 0.975)

p_A <- ggplot(pred_spp, aes(x = species, y = prob)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6, linewidth = 0.7) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Species", y = "Decay presence probability",
       tag = "A") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    plot.tag = element_text(face = "bold")
  )

# --- 3. PANEL B: Decay prevalence by site -----------------------------------

site_prev <- df %>%
  group_by(site) %>%
  summarise(
    n = n(),
    n_decay = sum(decay_binary),
    prop = n_decay / n,
    se = sqrt(prop * (1 - prop) / n),
    .groups = "drop"
  )

cat("\n=== Decay prevalence by site ===\n")
print(site_prev)

# Fisher exact test for site effect
site_tab <- table(df$site, df$decay_binary)
fisher_result <- fisher.test(site_tab)
cat("\nFisher exact test p =", fisher_result$p.value, "\n")

p_B <- ggplot(site_prev, aes(x = site, y = prop)) +
  geom_col(fill = "grey75", width = 0.6) +
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.15) +
  scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Site", y = "Proportion with decay",
       tag = "B") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    plot.tag = element_text(face = "bold")
  )

# --- 4. PANEL C: Decay severity ~ species (beta regression) -----------------

df_decayed <- df %>%
  filter(decay == "present") %>%
  mutate(
    damaged_prop = percent_damaged / 100,
    # Squeeze away from 0 and 1 for beta regression
    damaged_prop_adj = (damaged_prop * (n() - 1) + 0.5) / n(),
    # Ensure species keeps all levels for prediction
    species = factor(species, levels = levels(df$species))
  )

cat("\n=== Decay severity (percent_damaged) among decayed trees ===\n")
df_decayed %>%
  group_by(species, .drop = FALSE) %>%
  summarise(n = n(), mean_pct = mean(percent_damaged), sd_pct = sd(percent_damaged),
            .groups = "drop") %>%
  print()

betareg_sev <- betareg(damaged_prop_adj ~ species, data = df_decayed)
cat("\n=== Beta regression: severity ~ species ===\n")
print(summary(betareg_sev))

# Predicted severity per species (only species present in decayed data)
spp_in_model <- levels(droplevels(df_decayed$species))
pred_sev <- data.frame(
  species = factor(spp_in_model, levels = levels(df_decayed$species))
)
pred_sev$pred <- predict(betareg_sev, newdata = pred_sev, type = "response")

# Bootstrap CIs
boot_sev <- matrix(NA, nrow = n_boot, ncol = length(spp_in_model))
colnames(boot_sev) <- spp_in_model

for (i in seq_len(n_boot)) {
  idx <- sample(nrow(df_decayed), replace = TRUE)
  boot_df <- df_decayed[idx, ]
  fit <- tryCatch(
    betareg(damaged_prop_adj ~ species, data = boot_df),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    preds <- tryCatch(
      predict(fit, newdata = pred_sev, type = "response"),
      error = function(e) rep(NA, nrow(pred_sev))
    )
    boot_sev[i, ] <- preds
  }
}

pred_sev$lo <- apply(boot_sev, 2, quantile, 0.025, na.rm = TRUE)
pred_sev$hi <- apply(boot_sev, 2, quantile, 0.975, na.rm = TRUE)

p_C <- ggplot(pred_sev, aes(x = species, y = pred)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6, linewidth = 0.7) +
  labs(x = "Species", y = "Decay severity (proportion damaged)",
       tag = "C") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    plot.tag = element_text(face = "bold")
  )

# --- 5. PANEL D: Expected severity = P(decay) × E[severity | decay] ---------

# Only species present in both models
spp_both <- intersect(as.character(pred_spp$species), as.character(pred_sev$species))

pred_combined <- data.frame(
  species = factor(spp_both, levels = levels(df$species))
)
pred_combined <- pred_combined %>%
  left_join(pred_spp %>% select(species, prob), by = "species") %>%
  left_join(pred_sev %>% select(species, sev = pred), by = "species") %>%
  mutate(expected_sev = prob * sev)

# Bootstrap CIs for combined metric
boot_combined <- matrix(NA, nrow = n_boot, ncol = length(spp_both))
colnames(boot_combined) <- spp_both

for (i in seq_len(n_boot)) {
  p_presence <- boot_mat[i, spp_both]
  p_severity <- boot_sev[i, spp_both]
  if (!any(is.na(p_presence)) && !any(is.na(p_severity))) {
    boot_combined[i, ] <- p_presence * p_severity
  }
}

pred_combined$lo <- apply(boot_combined, 2, quantile, 0.025, na.rm = TRUE)
pred_combined$hi <- apply(boot_combined, 2, quantile, 0.975, na.rm = TRUE)

cat("\n=== Expected severity = P(decay) × E[severity | decay] ===\n")
print(pred_combined)

p_D <- ggplot(pred_combined, aes(x = species, y = expected_sev)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6, linewidth = 0.7) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Species", y = "Expected damage (proportion)",
       tag = "D") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    plot.tag = element_text(face = "bold")
  )

# --- 6. COMBINE AND SAVE ----------------------------------------------------

fig <- p_A | p_B | p_C | p_D
fig <- fig + plot_layout(widths = c(1, 0.7, 1, 1))

ggsave("output/figures/decay_presence_severity.pdf", fig, width = 15, height = 4.5)
cat("\nSaved: output/figures/decay_presence_severity.pdf\n")

# --- 7. SUPPLEMENTARY: Full model with site interaction ----------------------

cat("\n=== Binomial GLM: decay ~ species * site ===\n")
glm_full <- glm(decay_binary ~ species * site, data = df, family = binomial)
print(summary(glm_full))

# Compare models
cat("\nAIC comparison (additive vs interaction):\n")
cat("  species + site:", AIC(glm_presence), "\n")
cat("  species * site:", AIC(glm_full), "\n")

cat("\n=== ALL DONE ===\n")
