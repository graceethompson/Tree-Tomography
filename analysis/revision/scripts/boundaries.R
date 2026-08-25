# Boundary-by-boundary stability of the four-category phase diagram.
# R port of boundaries.py. Run from the repo root:
#   Rscript analysis/revision/scripts/boundaries.R

source("analysis/revision/scripts/revision_common.R")

# UTF-8 locale so the png devices render em-dash / sqrt / arrow glyphs
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"))) invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))

set.seed(0)
m <- load_merged()
pca <- build_pc1(m, "species")
pc <- pca$pc1
dmg <- m$percent_damaged
mu <- mean(pc)
sdv <- sd(pc)

classify <- function(struct, anom) four_cat(anom, struct)

cat_counts <- function(cc) {
  t <- table(factor(cc, levels = c("I", "II", "III", "IV")))
  paste(sprintf("%s=%d", names(t), as.integer(t)), collapse = " ")
}

# ---- the two axes have different shape ----
cat("SoT structural axis (percent_damaged): distribution\n")
cat("  values >0:", paste(sort(dmg[dmg > 0]), collapse = ", "), "\n")
bands <- cut(dmg, c(-.1, 0, 1, 3, 5, 10, 100),
             labels = c("0", "(0,1]", "(1,3]", "(3,5]", "(5,10]", ">10"))
bt <- table(bands)
cat("  bands:", paste(sprintf("'%s': %d", names(bt), as.integer(bt)), collapse = ", "), "\n")
cat(sprintf("ERT moisture axis: unimodal continuum, SD=%.2f, all trees have a value\n\n", sdv))

base <- classify(dmg > 1, pc > mu)
cat("Baseline categories:", cat_counts(base), "\n")

# ======== BOUNDARY-BY-BOUNDARY STABILITY ========
cat("\n", strrep("=", 80), " \nBOUNDARY STABILITY: which category edges move when thresholds move?\n ",
    strrep("=", 80), "\n", sep = "")

# --- SoT-axis boundaries: I<->IV (anom=0) and II<->III (anom=1). Vary tau_SoT 1..5% ---
cat("\n[SoT-axis boundaries]  I<->IV and II<->III  (vary structural-loss threshold 1%..5%)\n")
for (tau in c(0.5, 1, 2, 3, 5, 10)) {
  cc <- classify(dmg > tau, pc > mu)
  ch <- sum(cc != base)
  cat(sprintf("  tau_SoT=%4s%% : categories={%s}  trees changed vs baseline=%d\n",
              format(tau), cat_counts(cc), ch))
}
n_soT_band <- sum(dmg > 1 & dmg <= 5)
cat(sprintf("  -> trees in the 1%%-5%% 'ambiguous structural' band: %d  (gap: 0 trees in (0,1])\n",
            n_soT_band))

# --- ERT-axis boundaries: I<->II (struct=0) and III<->IV (struct=1). Vary tau_ERT mean±1SD ---
cat("\n[ERT-axis boundaries]  I<->II and III<->IV  (vary moisture-anomaly threshold mean±1SD)\n")
for (k in c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1)) {
  tau <- mu + k * sdv
  cc <- classify(dmg > 1, pc > tau)
  ch <- sum(cc != base)
  cat(sprintf("  tau_ERT=mean%+.2fSD : categories={%s}  trees changed vs baseline=%d\n",
              k, cat_counts(cc), ch))
}

# Decompose ERT instability by boundary, over mean±0.5SD
lo <- mu - 0.5 * sdv; hi <- mu + 0.5 * sdv
sound <- dmg <= 1; dmgd <- dmg > 1
inband <- (pc > lo) & (pc < hi)
cat(sprintf("\nERT 'ambiguous' band = mean±0.5SD = (%.2f,%.2f)\n", lo, hi))
cat(sprintf("  trees in band overall: %d\n", sum(inband)))
cat(sprintf("    of which structurally SOUND  (I<->II at risk): %d\n", sum(inband & sound)))
cat(sprintf("    of which structurally DAMAGED (III<->IV at risk): %d\n", sum(inband & dmgd)))
cat(sprintf("  structurally-damaged trees total: %d; their PC1 range: %.2f..%.2f\n",
            sum(dmgd), min(pc[dmgd]), max(pc[dmgd])))
cat(sprintf("  cavities (IV) at baseline: %d; active (III): %d\n",
            sum(base == "IV"), sum(base == "III")))

# explicit per-boundary flip counts over full ERT sweep and SoT sweep
flips_between <- function(a, b, base_c, alt_c) {
  sum((base_c == a & alt_c == b) | (base_c == b & alt_c == a))
}
cat("\nPer-boundary flips (baseline vs a perturbed threshold):\n")
altERT <- classify(dmg > 1, pc > (mu + 0.5 * sdv))
altSoT <- classify(dmg > 5, pc > mu)
cat(sprintf("  I<->II  (ERT axis, +0.5SD): %d\n", flips_between("I", "II", base, altERT)))
cat(sprintf("  III<->IV(ERT axis, +0.5SD): %d\n", flips_between("III", "IV", base, altERT)))
cat(sprintf("  II<->III(SoT axis, 1->5%%) : %d\n", flips_between("II", "III", base, altSoT)))
cat(sprintf("  I<->IV  (SoT axis, 1->5%%) : %d\n", flips_between("I", "IV", base, altSoT)))

# ======== FIGURE: 1-D stability of the structural (SoT) threshold ========
# Redesigned 2026-08-24: the moisture-threshold sensitivity is covered by the
# threshold-sweep figure; this figure focuses on the structural axis, whose
# strongly bimodal distribution is what makes the 1% cut stable.
nz <- sort(dmg[dmg > 0])
sx <- nz                                   # linear scale
st <- ave(nz, nz, FUN = seq_along)        # stack duplicates vertically
xt <- c(0, 1, 5, 10, 15, 20, 25, 30, 35)

png(file.path(OUT_DIR, "CJFR-boundary-fig.png"),
    width = 8.2, height = 3.4, units = "in", res = 170)
par(mar = c(4.0, 4.0, 1.2, 0.8), mgp = c(2.4, 0.7, 0))
plot(NA, xlim = c(-1, 36), ylim = c(0.4, 4.6), axes = FALSE,
     xlab = "", ylab = "")
usr <- par("usr")
rect(1, usr[3], 5, usr[4],
     col = adjustcolor("#7f8c8d", 0.15), border = NA)
abline(v = 1, lty = 2, col = "#333333", lwd = 1.5)
abline(v = 5, lty = 3, col = "#666666", lwd = 1.1)
# the 42 trees at exactly zero, shown as one annotated marker
points(0, 1, pch = 22, bg = "#3b6fb0", col = "white", cex = 2.4)
text(0, 1, "42", cex = 0.55, col = "white", font = 2)
text(0, 1.75, "42 trees at\nexactly 0%", cex = 0.66, col = "#26456e", adj = c(0.62, 0))
# the 15 damaged trees at their recorded values
points(sx, st, pch = 21, bg = "#b83232", col = "white", cex = 1.6, lwd = 0.6)
text(11.5, 4.25,
     "threshold range 1-5% (shaded):\nonly 2 trees inside (3% and 4%);\nclassifications otherwise identical",
     cex = 0.66, col = "#333333")
text(1, 0.55, "1% (manuscript)", cex = 0.62, col = "#333333", adj = 1.06)
text(5, 0.55, "5%", cex = 0.62, col = "#666666", adj = -0.45)
axis(1, at = xt, labels = xt, cex.axis = 0.85)
axis(2, at = 1:3, labels = 1:3, cex.axis = 0.8, las = 1)
box()
title(xlab = "SoT % of section damaged (recorded value)", cex.lab = 0.85)
title(ylab = "trees at value (stacked)", cex.lab = 0.85)
invisible(dev.off())
cat("\nsaved figure\n")
