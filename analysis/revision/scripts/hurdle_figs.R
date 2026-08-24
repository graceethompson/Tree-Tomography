# Hurdle-model observation and diagnostic figures (Reviewer 1, comments #6 & #7).
#
# The manuscript's hurdle model (code/decay_hurdle_analysis.R) reports only
# model-predicted means and intervals. The reviewer asked for:
#   #6 - the fitted logistic probabilities with the jittered binary observations
#        overlaid by species and site, making the sparse basis of the occurrence
#        model visible;
#   #7 - the beta-regression severity fit with the individual decayed trees
#        overlaid by species (showing per-species decayed n, incl. the single
#        N. sylvatica reference observation), plus residual diagnostics for the
#        severity stage matching the DHARMa checks reported for the binomial stage.
#
# Run from the repo root: Rscript analysis/revision/scripts/hurdle_figs.R

source("analysis/revision/scripts/revision_common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(betareg)
  library(emmeans)
  library(patchwork)
})

df <- load_main()
df$decay_binary <- as.integer(df$percent_damaged > 0)
df$species <- factor(df$species)           # levels bg, hem, rm, ro (bg = N.sylvatica reference)
df$site <- factor(df$site)
df$sp <- factor(df$sp, levels = unname(SPECIES_MAP[c("bg", "hem", "rm", "ro")]))

# ---- occurrence stage: manuscript model ----
m1 <- glm(decay_binary ~ species + site, family = binomial, data = df)

# ---- severity stage: manuscript model (decayed trees only, proportion scale) ----
sev <- df[df$percent_damaged > 0, ]
sev$prop <- sev$percent_damaged / 100
m2 <- betareg(prop ~ species, data = sev)

cat("Decayed n per species (severity-stage basis):\n")
print(table(droplevels(sev$sp)))
cat("\nOccurrence model (binomial GLM):\n"); print(summary(m1)$coefficients)
cat("\nSeverity model (beta regression):\n"); print(summary(m2)$coefficients$mean)

# =====================  Panel A: occurrence  =====================
# Fitted probability (+95% CI) for every species x site cell, with the raw
# binary observations jittered around 0/1.
cells <- as.data.frame(emmeans(m1, ~ species + site, type = "response"))
cells$sp <- factor(unname(SPECIES_MAP[as.character(cells$species)]), levels = levels(df$sp))
# drop cells with no data (specialists occur at one site only)
obs_n <- aggregate(decay_binary ~ species + site, df, length)
cells <- merge(cells, setNames(obs_n, c("species", "site", "n_obs")), by = c("species", "site"))

set.seed(42)
pa <- ggplot() +
  geom_jitter(data = df,
              aes(x = sp, y = decay_binary, colour = site),
              position = position_jitterdodge(jitter.width = 0.18, jitter.height = 0.04,
                                              dodge.width = 0.6),
              size = 1.6, alpha = 0.55) +
  geom_errorbar(data = cells,
                aes(x = sp, ymin = asymp.LCL, ymax = asymp.UCL, group = site),
                position = position_dodge(width = 0.6), width = 0.18, linewidth = 0.5) +
  geom_point(data = cells,
             aes(x = sp, y = prob, fill = site, group = site),
             position = position_dodge(width = 0.6), size = 3, shape = 21, colour = "black") +
  geom_text(data = cells,
            aes(x = sp, y = -0.13, label = paste0("n=", n_obs), group = site),
            position = position_dodge(width = 0.6), size = 2.7, colour = "grey30") +
  scale_colour_manual(values = c(BGS = "#1f77b4", EMS = "#d9822b")) +
  scale_fill_manual(values = c(BGS = "#1f77b4", EMS = "#d9822b")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.18, 1.08)) +
  labs(x = NULL, y = "decay present (0/1) and fitted probability",
       title = "A",
       colour = "site", fill = "site") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(face = "italic"))

# =====================  Panel B: severity  =====================
# NOTE: emmeans already returns betareg predictions on the RESPONSE (proportion)
# scale — do NOT apply plogis() again (code/decay_hurdle_analysis.R does, which
# double-transforms; its Panel C values ~0.5 are plogis(fitted proportion)).
emm2 <- as.data.frame(emmeans(m2, ~ species))
emm2$response <- emm2$emmean
emm2$LCL <- pmax(emm2$asymp.LCL, 0)   # delta-method CI can dip below 0
emm2$UCL <- emm2$asymp.UCL
emm2$sp <- factor(unname(SPECIES_MAP[as.character(emm2$species)]), levels = levels(df$sp))
nsev <- as.data.frame(table(sp = droplevels(sev$sp)))
emm2 <- merge(emm2, nsev, by = "sp", all.x = TRUE)
emm2$Freq[is.na(emm2$Freq)] <- 0

set.seed(7)
pb <- ggplot() +
  geom_jitter(data = sev, aes(x = sp, y = prop, colour = site),
              width = 0.12, height = 0, size = 2.2, alpha = 0.8) +
  geom_errorbar(data = emm2, aes(x = sp, ymin = LCL, ymax = UCL),
                width = 0.18, linewidth = 0.5) +
  geom_point(data = emm2, aes(x = sp, y = response),
             size = 3.2, shape = 21, fill = "grey20", colour = "black") +
  geom_text(data = emm2, aes(x = sp, y = -0.035, label = paste0("n=", Freq)),
            size = 2.9, colour = "grey30") +
  scale_colour_manual(values = c(BGS = "#1f77b4", EMS = "#d9822b")) +
  labs(x = NULL, y = "proportion of section damaged (decayed trees only)",
       title = "B",
       colour = "site") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(face = "italic"))

# =====================  Panels C & D: severity diagnostics  =====================
diag_df <- data.frame(
  fitted = fitted(m2),
  resid = residuals(m2, type = "sweighted2"),
  leverage = hatvalues(m2)
)
# the single N. sylvatica tree is its own reference level: leverage = 1, so its
# standardized residual is undefined (NaN) — the instability the reviewer flags
n_undef <- sum(!is.finite(diag_df$resid))
diag_df <- diag_df[is.finite(diag_df$resid), ]
pc <- ggplot(diag_df, aes(fitted, resid)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "fitted proportion damaged", y = "standardized weighted residual (type 2)",
       title = "C") +
  theme_classic(base_size = 10)

qq <- qqnorm(diag_df$resid, plot.it = FALSE)
pd <- ggplot(data.frame(x = qq$x, y = qq$y), aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "theoretical quantiles", y = "sample quantiles",
       title = "D") +
  theme_classic(base_size = 10)

fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1.35, 1))
ggsave(file.path(OUT_DIR, "CJFR-hurdle-fig.png"), fig, width = 12.5, height = 9, dpi = 170)

# DHARMa simulated residuals for the occurrence stage (as in the manuscript),
# saved alongside for completeness.
suppressPackageStartupMessages(library(DHARMa))
png(file.path(OUT_DIR, "CJFR-hurdle-dharma.png"), width = 1700, height = 850, res = 150)
plot(simulateResiduals(m1, seed = 1))
dev.off()

cat("\nsaved", file.path(OUT_DIR, "CJFR-hurdle-fig.png"), "and CJFR-hurdle-dharma.png\n")
