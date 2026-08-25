# Where does 'incipient' start? Sensitivity of the sound->incipient cut and a
# 6-cell dead-band scheme. R port of incipient.py. Run from the repo root:
#   Rscript analysis/revision/scripts/incipient.R

source("analysis/revision/scripts/revision_common.R")

# UTF-8 locale so the png devices render em-dash / sqrt / arrow glyphs
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"))) invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))

m <- load_merged()
# species-standardized resistivity deviation; ORIENT so higher = WETTER (right)
m$xdev <- xdev_species_median(m)
m$struct <- m$percent_damaged > 1
sound <- !m$struct; dmgd <- m$struct
mb <- m$site == "BGS"; me <- m$site == "EMS"

# ---- WHERE DOES INCIPIENT START? sensitivity of the sound->incipient cut ----
cat("Incipient onset = sound tree confidently WETTER than its species median.\n")
cat("(xdev = species-standardized; wetter = higher; +0.5 SD = confidently wet)\n")
cat(sprintf("%-32s%5s%7s%7s\n", "cut (SD wetter than median)", "nII", "BGS%", "EMS%"))
cuts <- list(c(0, "median (any wetter)"), c(0.25, "+0.25 SD"),
             c(0.5, "+0.5 SD"), c(1.0, "+1.0 SD"))
for (cl in cuts) {
  cut_ <- as.numeric(cl[1]); lab <- cl[2]
  inc <- sound & (m$xdev > cut_)
  b <- sum(inc[mb]) / sum(sound[mb]) * 100
  e <- sum(inc[me]) / sum(sound[me]) * 100
  cat(sprintf("  %-30s%5d%7.0f%7.0f\n", lab, sum(inc), b, e))
}

# ---- internal 'positive control': moisture of structurally-confirmed decay (active) trees ----
# active = damaged & wet (xdev>0). Their xdev shows 'decay-associated moisture' range.
active <- m[dmgd & (m$xdev > 0), ]
cat(sprintf("\nInternal reference — structurally-damaged WET (active) trees, n=%d:\n", nrow(active)))
cat(sprintf("  their xdev (wetness): median=%.2f, range %.2f..%.2f\n",
            median(active$xdev), min(active$xdev), max(active$xdev)))
cat("  -> a sound tree wetter than ~+0.5 SD overlaps the active-decay moisture range (supports +0.5 SD onset)\n")

# ---- 6-cell / dead-band classification ----
classify <- function(cut_) {
  wet <- m$xdev > cut_; dry <- m$xdev < -cut_
  ifelse(sound & wet, "II",
    ifelse(sound & !wet, "I",
      ifelse(dmgd & dry, "IV", "III")))
}
cat_counts <- function(cc) {
  t <- table(factor(cc, levels = c("I", "II", "III", "IV")))
  paste(sprintf("%s=%d", names(t), as.integer(t)), collapse = " ")
}
for (cut_ in c(0.5, 1.0)) {
  cat(sprintf("\n6-cell scheme, dead band ±%s SD:  {%s}\n",
              sprintf("%.1f", cut_), cat_counts(classify(cut_))))
}

# ---- FIGURE: dead-band phase diagram, nested band widths in one panel ----
# Rationale for the widths (PC1 sd = 1.89): the alternative data-driven
# thresholds sit at -0.13 SD (median), +0.53 SD (Otsu) and +1.03 SD (GMM
# crossover) from the mean, so +/-0.5 SD spans the median-to-Otsu range of
# defensible cuts and +/-1 SD additionally covers the GMM crossover.
set.seed(2)
yt <- c(0, 1, 5, 10, 20, 30)
y <- sqrt(pmax(m$percent_damaged, 0))
xr <- range(m$xdev); xpad <- 0.05 * diff(xr)
xlim <- c(xr[1] - xpad, xr[2] + xpad)

png(file.path(OUT_DIR, "CJFR-6cell.png"),
    width = 8.4, height = 5.6, units = "in", res = 170)
par(mar = c(4.4, 4.2, 1.2, 0.8), mgp = c(2.4, 0.7, 0))
cat_ <- classify(0)   # baseline species-median split; bands overlay uncertainty
plot(NA, xlim = xlim, ylim = c(-0.2, sqrt(37) + 0.2), axes = FALSE,
     xlab = "", ylab = "")
usr <- par("usr")
rect(-1.0, usr[3], 1.0, usr[4], col = adjustcolor("gray82", alpha.f = 0.45), border = NA)
rect(-0.5, usr[3], 0.5, usr[4], col = adjustcolor("gray62", alpha.f = 0.45), border = NA)
abline(v = 0, col = "#666666", lty = 1, lwd = 1)
abline(h = sqrt(1), col = "#666666", lty = 1, lwd = 1)
for (cc in names(CAT_COLS)) {
  idx <- cat_ == cc
  points(m$xdev[idx], y[idx] + rnorm(sum(idx), 0, 0.03),
         pch = 21, bg = CAT_COLS[[cc]], col = "white", lwd = 0.5, cex = 1.3)
}
axis(1, cex.axis = 0.8)
axis(2, at = sqrt(yt), labels = yt, cex.axis = 0.8, las = 1)
box()
title(xlab = "moisture anomaly  (species-standardized;  drier \u2190   0 = species median   \u2192 wetter)",
      cex.lab = 0.8)
title(ylab = "SoT structural loss %  (\u221a scale)", cex.lab = 0.85)
xl <- xlim
text(xl[1] + 0.2, sqrt(31), "IV Cavity", cex = 0.75,
     col = CAT_COLS[["IV"]], font = 2, adj = c(0, 1))
text(xl[2] - 0.2, sqrt(31), "III Active", cex = 0.75,
     col = CAT_COLS[["III"]], font = 2, adj = c(1, 1))
text(xl[1] + 0.2, sqrt(2.6), "I No decay", cex = 0.75,
     col = CAT_COLS[["I"]], font = 2, adj = c(0, 0.5))
text(xl[2] - 0.2, sqrt(2.6), "II Incipient", cex = 0.75,
     col = CAT_COLS[["II"]], font = 2, adj = c(1, 0.5))
text(0, sqrt(35.5), "transitional \u00b10.5 SD", cex = 0.62, col = "#444444")
text(0.75, sqrt(33), "\u00b11 SD", cex = 0.62, col = "#777777", adj = 0)
legend("left", legend = names(CAT_COLS), pch = 21, pt.bg = unname(CAT_COLS),
       col = "white", pt.cex = 1.25, cex = 0.75, bty = "n", title = "Category")
invisible(dev.off())
cat("\nsaved 6-cell figure\n")
