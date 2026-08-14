# R port of reanalysis.py — revision reanalysis (Reviewers #3, #6, #8).
# Run from repo root: Rscript analysis/revision/scripts/reanalysis.R
source("analysis/revision/scripts/revision_common.R")

sd_pop <- function(x) sqrt(mean((x - mean(x))^2))  # numpy default ddof=0
hr <- strrep("=", 70)

df <- load_main()
df$decay_bin <- as.integer(df$percent_damaged > 0)

cat(hr, "\nN by species x site\n", hr, "\n", sep = "")
print(addmargins(table(df$sp, df$site)))
cat("\nDecay present counts by species x site:\n")
print(ftable(xtabs(~ sp + site + decay_bin, df), col.vars = "decay_bin"))

# ---------- A. STEM SIZE CONFOUND (Reviewer #8) ----------
cat("\n", hr, "\nA. DBH by site & species (Reviewer #8 stem-size confound)\n", hr, "\n", sep = "")
summ <- function(x) c(count = length(x), mean = mean(x), std = sd(x), min = min(x), max = max(x))
print(round(t(sapply(split(df$dbh, df$site), summ)), 6))
print(round(t(sapply(split(df$dbh, df$sp), summ)), 6))

tt <- t.test(df$dbh[df$site == "BGS"], df$dbh[df$site == "EMS"], var.equal = TRUE)
cat(sprintf("\nDBH BGS vs EMS: BGS mean=%.1f, EMS mean=%.1f, t-test p=%.4f\n",
            mean(df$dbh[df$site == "BGS"]), mean(df$dbh[df$site == "EMS"]), tt$p.value))
r <- cor.test(df$dbh, df$percent_damaged)
rs <- suppressWarnings(cor.test(df$dbh, df$percent_damaged, method = "spearman"))
cat(sprintf("dbh vs percent_damaged: Pearson r=%.3f p=%.3f; Spearman rho=%.3f p=%.3f\n",
            r$estimate, r$p.value, rs$estimate, rs$p.value))
d <- df[df$decay_bin == 1, ]
r2 <- cor.test(d$dbh, d$percent_damaged)
cat(sprintf("among decayed (n=%d): dbh vs pct_damaged Pearson r=%.3f p=%.3f\n",
            nrow(d), r2$estimate, r2$p.value))
mn <- round(tapply(df$dbh, df$decay_bin, mean), 1)
cat(sprintf("dbh mean by decay presence: {0: %.1f, 1: %.1f}\n", mn["0"], mn["1"]))
glm_dbh <- glm(decay_bin ~ dbh, df, family = binomial())
sc <- summary(glm_dbh)$coefficients
cat(sprintf("decay~dbh logistic: beta=%.3f p=%.3f\n", sc["dbh", "Estimate"], sc["dbh", "Pr(>|z|)"]))

# ---------- B. ANOVA one-way vs two-way (Reviewer #6) ----------
cat("\n", hr, "\nB. ANOVA on percent_damaged (Reviewer #6)\n", hr, "\n", sep = "")
m1 <- lm(percent_damaged ~ site, df)
cat("One-way (site):\n"); print(round(as.data.frame(anova(m1)), 4))
m2 <- lm(percent_damaged ~ sp + site, df)
cat("\nTwo-way additive (species+site):\n"); print(round(as.data.frame(anova(m2)), 4))
m3 <- lm(percent_damaged ~ sp * site, df)
# NOTE: the sp:site design is rank-deficient (N.sylvatica occurs only at BGS,
# Q.rubra only at EMS), so 2 of the 3 interaction columns are aliased. R's anova
# drops them (df=1, SS=12.4082, p=0.6908); statsmodels' anova_lm keeps df=3 and
# reports SS=15.3824 (p=0.9776), a value that depends on BLAS rounding of the
# degenerate QR directions and is not reproducible across implementations.
# All other rows (sp, site, residuals) are identical. Either way the
# interaction is far from significant.
cat("\nTwo-way with interaction:\n"); print(round(as.data.frame(anova(m3)), 4))

# ---------- C. Binomial GLM presence (site vs species) ----------
cat("\n", hr, "\nC. Decay presence GLM (site vs species)\n", hr, "\n", sep = "")
gnull <- glm(decay_bin ~ 1, df, family = binomial())
gsite <- glm(decay_bin ~ site, df, family = binomial())
gspp  <- glm(decay_bin ~ sp, df, family = binomial())
gfull <- glm(decay_bin ~ sp + site, df, family = binomial())
lrt <- function(a, b) {
  stat <- 2 * (as.numeric(logLik(b)) - as.numeric(logLik(a)))
  dfd <- attr(logLik(b), "df") - attr(logLik(a), "df")
  c(stat, dfd, pchisq(stat, dfd, lower.tail = FALSE))
}
fmt_lrt <- function(x) sprintf("[%.4f, %d, %.4f]", x[1], x[2], x[3])
cat("Null->Site LRT:", fmt_lrt(lrt(gnull, gsite)), "\n")
cat("Site->Full (add species):", fmt_lrt(lrt(gsite, gfull)), "\n")
cat("Null->Species LRT:", fmt_lrt(lrt(gnull, gspp)), "\n")
cat("Species->Full (add site):", fmt_lrt(lrt(gspp, gfull)), "\n")
scf <- summary(gfull)$coefficients
cat("Full model site coef:", round(scf["siteEMS", "Estimate"], 3),
    "p=", round(scf["siteEMS", "Pr(>|z|)"], 4), "\n")
prev <- round(tapply(df$decay_bin, df$site, mean), 3)
cat(sprintf("Prevalence by site: {'BGS': %.3f, 'EMS': %.3f}\n", prev["BGS"], prev["EMS"]))

# ---------- D. Threshold sensitivity: reconstruct PC1 & classification ----------
cat("\n", hr, "\nD. Reconstruct ERT PC1 (species z-scores) + classification thresholds\n", hr, "\n", sep = "")
m <- load_merged()
cat(sprintf("Merged main trees with ERT: %d\n", nrow(m)))
bp <- build_pc1(m, "species")
pc1 <- bp$pc1
cat(sprintf("PC1 var explained=%.1f%%, PC2=%.1f%%\n", bp$ve[1] * 100, bp$ve[2] * 100))
thr <- mean(pc1)
m$struct_loss <- m$percent_damaged > 1
m$anom <- pc1 > thr
cls_map <- c(I = "I_NoDecay", II = "II_Incipient", III = "III_Active", IV = "IV_Cavity")
m$cls <- unname(cls_map[four_cat(m$anom, m$struct_loss)])
cat("Classification counts (threshold=mean, PC1>mean):\n")
ord <- unique(m$cls)
vc <- table(factor(m$cls, levels = ord))
vc <- vc[order(-vc)]  # stable sort: ties keep first-appearance order (as pandas)
for (nm in names(vc)) cat(sprintf("%-13s%d\n", nm, vc[nm]))
cat("By site:\n")
pt <- round(prop.table(table(m$site, m$cls), 1), 2)
print(pt[, c("III_Active", "II_Incipient", "IV_Cavity", "I_NoDecay")])

cat("\nThreshold sensitivity for 'incipient' assignment among structurally-sound trees:\n")
sd1 <- sd_pop(pc1)  # numpy ndarray.std(), ddof=0
sens <- list(list(-0.5, "mean-0.5SD"), list(0, "mean"), list(0.5, "mean+0.5SD"), list(1, "mean+1SD"))
for (s in sens) {
  t2 <- mean(pc1) + s[[1]] * sd1
  anom <- pc1 > t2
  inc <- sum(!m$struct_loss & anom)
  cat(sprintf("  %s: n_incipient=%d, n_anomalous_total=%d / %d\n", s[[2]], inc, sum(anom), nrow(m)))
}
anom2 <- pc1 > (mean(pc1) + 0.5 * sd1)
cls2 <- unname(cls_map[four_cat(anom2, m$struct_loss)])
flips <- sum(m$cls != cls2)
cat(sprintf("Trees changing class when threshold -> mean+0.5SD: %d of %d\n", flips, nrow(m)))

cat("\n1% SoT threshold: trees by damage band:\n")
bands <- cut(df$percent_damaged, c(-0.1, 0, 1, 5, 100), labels = c("0", "(0-1]", "(1-5]", ">5"))
print(table(bands))

# ---------- E. Validation: PC1 vs mean resistivity vs moisture (Reviewer #3) ----------
cat("\n", hr, "\nE. Hemlock validation: PC1 vs mean resistivity vs moisture (Reviewer #3)\n", hr, "\n", sep = "")
v <- load_validation()
vmet <- c("mean_res", "Median", "SD", "CV", "Gini", "Entropy", "CMA", "RadialGradient")
vv <- v
vv$mean <- vv$mean_res  # so build_pc1 orients PC1 against mean resistivity
bpv <- build_pc1(vv, "pooled", cols = vmet)
vpc1 <- bpv$pc1
rpc <- cor.test(vpc1, v$moisture)
rmean <- cor.test(v$mean_res, v$moisture)
cat(sprintf("n=%d hemlocks\n", nrow(v)))
cat(sprintf("PC1 vs moisture: r=%.3f p=%.3f\n", rpc$estimate, rpc$p.value))
cat(sprintf("mean resistivity vs moisture: r=%.3f p=%.3f\n", rmean$estimate, rmean$p.value))
cat(sprintf("Conductance vs moisture: r=%.3f\n", cor(v$Conductance, v$moisture)))
for (met in c("Median", "CV", "Gini", "Entropy", "CMA", "RadialGradient")) {
  rr <- cor.test(v[[met]], v$moisture)
  cat(sprintf("  %s vs moisture: r=%.3f p=%.3f\n", met, rr$estimate, rr$p.value))
}

# ---------- F. Carbon: 47% vs 26% and mean damage among structural-loss ----------
cat("\n", hr, "\nF. Carbon framing numbers\n", hr, "\n", sep = "")
cat(sprintf("%% trees any anomaly (incipient+active+cavity) = %.0f%%\n", mean(m$cls != "I_NoDecay") * 100))
cat(sprintf("%% trees structural loss (>1%% SoT) = %.0f%%\n", mean(df$percent_damaged > 1) * 100))
sl <- df[df$percent_damaged > 1, ]
cat(sprintf("Among structural-loss trees: mean pct damage overall=%.1f%%\n", mean(sl$percent_damaged)))
bys <- round(tapply(sl$percent_damaged, sl$site, mean), 1)
cat(sprintf("  by site: {'BGS': %.1f, 'EMS': %.1f}\n", bys["BGS"], bys["EMS"]))
ball <- round(tapply(df$percent_damaged, df$site, mean), 2)
cat(sprintf("Mean pct damage all trees by site: {'BGS': %.2f, 'EMS': %.2f}\n", ball["BGS"], ball["EMS"]))
