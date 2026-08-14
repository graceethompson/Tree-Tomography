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

# ======== FIGURE: phase diagram with threshold BANDS ========
add_alpha <- function(col, a) adjustcolor(col, alpha.f = a)
ysqrt <- sqrt(pmax(dmg, 0))
yt <- c(0, 1, 5, 10, 20, 30)
mk <- c(A.rubrum = 24, T.canadensis = 21, N.sylvatica = 22, Q.rubra = 23)

png(file.path(OUT_DIR, "CJFR-boundary-fig.png"),
    width = 8.2, height = 6.4, units = "in", res = 170)
par(mar = c(4.2, 4.2, 2.6, 1), mgp = c(2.4, 0.7, 0))
xr <- range(pc); xpad <- 0.05 * diff(xr)
plot(NA, xlim = c(xr[1] - xpad, xr[2] + xpad), ylim = c(-0.2, sqrt(36) + 0.35),
     axes = FALSE, xlab = "", ylab = "")
# bands
usr <- par("usr")
rect(lo, usr[3], hi, usr[4], col = add_alpha("#f39c12", 0.13), border = NA)      # ERT ambiguous band (wide)
rect(usr[1], sqrt(1), usr[2], sqrt(5), col = add_alpha("#7f8c8d", 0.18), border = NA) # SoT ambiguous band
abline(v = mu, col = "#c0392b", lty = 2, lwd = 1.4)
abline(h = sqrt(1), col = "#333333", lty = 2, lwd = 1.4)
for (s in names(SPECIES_COLS)) {
  idx <- m$sp == s
  points(pc[idx], ysqrt[idx] + rnorm(sum(idx), 0, 0.03),
         pch = mk[[s]], bg = SPECIES_COLS[[s]], col = "white", lwd = 0.5, cex = 1.25)
}
axis(1, cex.axis = 0.85)
axis(2, at = sqrt(yt), labels = yt, cex.axis = 0.85, las = 1)
box()
title(xlab = "ERT PC1  (moisture axis — CONTINUUM, no natural break)", cex.lab = 0.85)
title(ylab = "SoT structural loss %  (√ scale — BIMODAL, gap at 0→6%)", cex.lab = 0.85)
# quadrant labels
text(min(pc) + 0.2, sqrt(28), "IV: Cavity", cex = 0.85, col = "#7f0000", font = 2, adj = 0)
text(max(pc) - 1.4, sqrt(28), "III: Active", cex = 0.85, col = "#7f4f00", font = 2, adj = 0)
text(min(pc) + 0.2, 0.05, "I: No decay", cex = 0.85, col = "#26456e", font = 2, adj = 0)
text(max(pc) - 1.6, 0.05, "II: Incipient", cex = 0.85, col = "#4d7f00", font = 2, adj = 0)
text(mu, sqrt(34) + 0.28, "moisture threshold\n(wide band = many trees flip)",
     cex = 0.62, col = "#c0392b")
text(max(pc) - 0.1, sqrt(3), "structural threshold\n(band empty = stable)",
     cex = 0.62, col = "#333333", adj = 1)
legend("right", legend = names(SPECIES_COLS), pch = mk[names(SPECIES_COLS)],
       pt.bg = unname(SPECIES_COLS[names(SPECIES_COLS)]), col = "white",
       pt.cex = 1.2, cex = 0.72, bty = "n")
title(main = "Boundary stability: the structural (SoT) split is crisp; the moisture (ERT) split is fuzzy",
      cex.main = 0.85)
invisible(dev.off())
cat("\nsaved figure\n")
