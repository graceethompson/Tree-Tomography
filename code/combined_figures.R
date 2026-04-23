# =============================================================================
# Combined layouts for paper Figures 2, 3, and 4
#   Fig 2: SoT % damaged wood — (A) by site, (B) by species × site
#   Fig 3: ERT metrics — (A,B) mean resistivity; (C,D) coefficient of variation
#   Fig 4: Hurdle model — produced in 3-panel (A/B/C) and 4-panel (A/B/C/D) forms
# Self-contained: reads raw data, refits models, writes PDF + PNG.
# =============================================================================

library(tidyverse)
library(patchwork)
library(betareg)
library(scales)

# --- 1. DATA ----------------------------------------------------------------

df <- read_csv("data/Tree_ID_info.csv", show_col_types = FALSE) %>%
  mutate(
    site = if_else(plot %in% c("C1","C2","D1","D2","E1","E5","F1","G1","H1"),
                   "EMS", "BGS"),
    species = factor(species, levels = c("hem", "rm", "bg", "ro")),
    decay_binary = as.integer(decay == "present")
  )

ert <- read_csv("data/ERT_application_results.csv", show_col_types = FALSE)
combined <- df %>% left_join(ert, by = "tree")

# --- 2. SHARED AESTHETICS ---------------------------------------------------

# Muted, nature-inspired palette — low saturation, easy on the eye:
#   hem (T. canadensis) — sage green (hemlock foliage)
#   rm  (A. rubrum)     — muted terracotta (fall color)
#   bg  (N. sylvatica)  — dusty lavender (fall foliage)
#   ro  (Q. rubra)      — warm taupe (oak bark)
species_colors <- c(
  hem = "#6B8C7A",
  rm  = "#B5725C",
  bg  = "#8B7FAA",
  ro  = "#9E8060"
)
species_labels <- c(
  hem = expression(italic("T. canadensis")),
  rm  = expression(italic("A. rubrum")),
  bg  = expression(italic("N. sylvatica")),
  ro  = expression(italic("Q. rubra"))
)

# Site palette: slate blue for wetland (BGS), warm sand for upland (EMS)
site_colors <- c(BGS = "#6890A4", EMS = "#C4A070")

base_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.tag = element_text(face = "bold")
  )

# --- 3. FIGURE 2: SoT % damaged wood ----------------------------------------

p2A <- ggplot(df, aes(site, percent_damaged, color = site, fill = site)) +
  geom_boxplot(linewidth = 0.6, alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = alpha(site_colors, 0.2)) +
  labs(x = "Site", y = "Damaged Wood (%)", tag = "A") +
  base_theme

p2B <- ggplot(df, aes(species, percent_damaged, color = species, fill = species)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") +
  scale_x_discrete(labels = species_labels) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = alpha(species_colors, 0.2)) +
  labs(x = "Species", y = NULL, tag = "B") +
  base_theme

fig2 <- p2A + p2B + plot_layout(widths = c(1, 2))

ggsave("output/figures/fig2_sot_damaged_combined.pdf", fig2,
       width = 12, height = 5)
ggsave("output/figures/fig2_sot_damaged_combined.png", fig2,
       width = 12, height = 5, dpi = 300)

# --- 4. FIGURE 3: ERT mean resistivity (row 1) + CV (row 2) -----------------

# Row 1: Mean resistivity (Ohm.m)
p3A <- ggplot(combined, aes(site, mean, color = site, fill = site)) +
  geom_boxplot(linewidth = 0.6, alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = alpha(site_colors, 0.2)) +
  labs(x = "Site",
       y = expression("Mean Resistivity ("*Omega %.% m*")"),
       tag = "A") +
  base_theme

p3B <- ggplot(combined, aes(species, mean, color = species, fill = species)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") +
  scale_x_discrete(labels = species_labels) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = alpha(species_colors, 0.2)) +
  labs(x = "Species", y = NULL, tag = "B") +
  base_theme

# Row 2: Coefficient of variation
p3C <- ggplot(combined, aes(site, cv, color = site, fill = site)) +
  geom_boxplot(linewidth = 0.6, alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = alpha(site_colors, 0.2)) +
  labs(x = "Site", y = "Coefficient of Variation", tag = "C") +
  base_theme

p3D <- ggplot(combined, aes(species, cv, color = species, fill = species)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ site, scales = "free_x") +
  scale_x_discrete(labels = species_labels) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = alpha(species_colors, 0.2)) +
  labs(x = "Species", y = NULL, tag = "D") +
  base_theme

fig3 <- wrap_plots(p3A, p3B, p3C, p3D, ncol = 2, widths = c(1, 2))

ggsave("output/figures/fig3_ert_mean_cv_combined.pdf", fig3,
       width = 12, height = 10)
ggsave("output/figures/fig3_ert_mean_cv_combined.png", fig3,
       width = 12, height = 10, dpi = 300)

# --- 5. FIGURE 4: Hurdle model panels --------------------------------------
# For consistency with paper ordering, use species levels c("rm","bg","ro","hem")
# for the hurdle analysis (matches decay_presence_severity.R).

df_h <- df %>%
  mutate(species = factor(as.character(species),
                          levels = c("rm", "bg", "ro", "hem")))

# Panel A: binomial GLM for decay presence
glm_presence <- glm(decay_binary ~ species + site, data = df_h, family = binomial)

newdata_spp <- expand.grid(
  species = levels(df_h$species),
  site = c("BGS", "EMS")
)
newdata_spp$pred <- predict(glm_presence, newdata = newdata_spp, type = "response")

pred_spp <- newdata_spp %>%
  group_by(species) %>%
  summarise(prob = mean(pred), .groups = "drop")

set.seed(42)
n_boot <- 2000
boot_mat <- matrix(NA, nrow = n_boot, ncol = nlevels(df_h$species))
colnames(boot_mat) <- levels(df_h$species)

for (i in seq_len(n_boot)) {
  idx <- sample(nrow(df_h), replace = TRUE)
  boot_df <- df_h[idx, ]
  fit <- suppressWarnings(
    glm(decay_binary ~ species + site, data = boot_df, family = binomial)
  )
  preds <- predict(fit, newdata = newdata_spp, type = "response")
  marg <- tapply(preds, newdata_spp$species, mean)
  boot_mat[i, ] <- marg[levels(df_h$species)]
}

pred_spp$lo <- apply(boot_mat, 2, quantile, 0.025)
pred_spp$hi <- apply(boot_mat, 2, quantile, 0.975)

p4A <- ggplot(pred_spp, aes(x = species, y = prob, color = species)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.7) +
  scale_x_discrete(labels = species_labels) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "Species", y = "Decay presence probability", tag = "A") +
  base_theme

# Panel B: observed prevalence by site
site_prev <- df_h %>%
  group_by(site) %>%
  summarise(
    n = n(),
    n_decay = sum(decay_binary),
    prop = n_decay / n,
    se = sqrt(prop * (1 - prop) / n),
    .groups = "drop"
  )

p4B <- ggplot(site_prev, aes(x = site, y = prop, fill = site)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.15) +
  scale_fill_manual(values = site_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Site", y = "Proportion with decay", tag = "B") +
  base_theme

# Panel C: beta regression severity among decayed trees
# Uses raw proportions with alphabetical factor ordering (bg = N. sylvatica as reference)
# to match the model reported in the manuscript text (beta = 1.70 for rm, 1.53 for ro).
df_decayed <- df_h %>%
  filter(decay == "present") %>%
  mutate(
    damaged_prop = percent_damaged / 100,
    species = factor(as.character(species))  # alphabetical: bg, hem, rm, ro
  )

betareg_sev <- betareg(damaged_prop ~ species, data = df_decayed)

spp_in_model <- levels(df_decayed$species)
pred_sev <- data.frame(
  species = factor(spp_in_model, levels = spp_in_model)
)
pred_sev$pred <- predict(betareg_sev, newdata = pred_sev, type = "response")

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

p4C <- ggplot(pred_sev, aes(x = species, y = pred, color = species)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.7) +
  scale_x_discrete(labels = species_labels) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "Species",
       y = "Decay severity (proportion damaged)",
       tag = "C") +
  base_theme

# Panel D: expected severity = P(decay) × E[severity | decay]
spp_both <- intersect(as.character(pred_spp$species),
                      as.character(pred_sev$species))

pred_combined_sev <- data.frame(
  species = factor(spp_both, levels = levels(df_h$species))
) %>%
  left_join(pred_spp %>% select(species, prob), by = "species") %>%
  left_join(pred_sev %>% select(species, sev = pred), by = "species") %>%
  mutate(expected_sev = prob * sev)

boot_combined_sev <- matrix(NA, nrow = n_boot, ncol = length(spp_both))
colnames(boot_combined_sev) <- spp_both

for (i in seq_len(n_boot)) {
  p_presence <- boot_mat[i, spp_both]
  p_severity <- boot_sev[i, spp_both]
  if (!any(is.na(p_presence)) && !any(is.na(p_severity))) {
    boot_combined_sev[i, ] <- p_presence * p_severity
  }
}

pred_combined_sev$lo <- apply(boot_combined_sev, 2, quantile, 0.025, na.rm = TRUE)
pred_combined_sev$hi <- apply(boot_combined_sev, 2, quantile, 0.975, na.rm = TRUE)

p4D <- ggplot(pred_combined_sev, aes(x = species, y = expected_sev, color = species)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.7) +
  scale_x_discrete(labels = species_labels) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "Species",
       y = "Expected damage (proportion)",
       tag = "D") +
  base_theme

# 3-panel version (paper text: A, B, C)
fig4_3 <- p4A + p4B + p4C + plot_layout(widths = c(1, 0.7, 1))
ggsave("output/figures/fig4_hurdle_3panel.pdf", fig4_3,
       width = 12, height = 4.5)
ggsave("output/figures/fig4_hurdle_3panel.png", fig4_3,
       width = 12, height = 4.5, dpi = 300)

# 4-panel version (A, B, C, D)
fig4_4 <- p4A + p4B + p4C + p4D + plot_layout(widths = c(1, 0.7, 1, 1))
ggsave("output/figures/fig4_hurdle_4panel.pdf", fig4_4,
       width = 15, height = 4.5)
ggsave("output/figures/fig4_hurdle_4panel.png", fig4_4,
       width = 15, height = 4.5, dpi = 300)

cat("\nWrote:\n",
    " output/figures/fig2_sot_damaged_combined.{pdf,png}\n",
    " output/figures/fig3_ert_mean_cv_combined.{pdf,png}\n",
    " output/figures/fig4_hurdle_3panel.{pdf,png}\n",
    " output/figures/fig4_hurdle_4panel.{pdf,png}\n", sep = "")
