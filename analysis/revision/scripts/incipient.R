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

# ---- FIGURE: 6-cell dead-band phase diagram (two dead-band widths) ----
set.seed(2)
yt <- c(0, 1, 5, 10, 20, 30)
y <- sqrt(pmax(m$percent_damaged, 0))
xr <- range(m$xdev); xpad <- 0.05 * diff(xr)
xlim <- c(xr[1] - xpad, xr[2] + xpad)

png(file.path(OUT_DIR, "CJFR-6cell.png"),
    width = 15, height = 6, units = "in", res = 140)
par(mfrow = c(1, 2), mar = c(4.4, 4.2, 2.4, 0.8), oma = c(0, 0, 3.4, 0),
    mgp = c(2.4, 0.7, 0))
first <- TRUE
for (cut_ in c(0.5, 1.0)) {
  cat_ <- classify(cut_)
  plot(NA, xlim = xlim, ylim = c(-0.2, sqrt(37) + 0.2), axes = FALSE,
       xlab = "", ylab = "")
  rect(-cut_, par("usr")[3], cut_, par("usr")[4],
       col = adjustcolor("gray85", alpha.f = 0.6), border = NA)
  abline(v = 0, col = "#666666", lty = 1, lwd = 1)
  abline(h = sqrt(1), col = "#666666", lty = 1, lwd = 1)
  for (cc in names(CAT_COLS)) {
    idx <- cat_ == cc
    points(m$xdev[idx], y[idx] + rnorm(sum(idx), 0, 0.03),
           pch = 21, bg = CAT_COLS[[cc]], col = "white", lwd = 0.5, cex = 1.3)
  }
  axis(1, cex.axis = 0.8)
  if (first) axis(2, at = sqrt(yt), labels = yt, cex.axis = 0.8, las = 1) else
    axis(2, at = sqrt(yt), labels = FALSE)
  box()
  title(xlab = "moisture anomaly  (species-standardized;  drier ←   0 = species median   → wetter)",
        cex.lab = 0.8)
  if (first) title(ylab = "SoT structural loss %  (√ scale)", cex.lab = 0.85)
  # cell labels
  xl <- xlim
  text(xl[1] + 0.2, sqrt(31), "IV Cavity\n(damaged, confidently dry)", cex = 0.7,
       col = CAT_COLS[["IV"]], font = 2, adj = c(0, 1))
  text(xl[2] - 0.2, sqrt(31), "III Active\n(damaged, wet)", cex = 0.7,
       col = CAT_COLS[["III"]], font = 2, adj = c(1, 1))
  text(xl[1] + 0.2, 0.05, "I No decay\n(sound, not-wet)", cex = 0.7,
       col = CAT_COLS[["I"]], font = 2, adj = c(0, 0.5))
  text(xl[2] - 0.2, 0.05, "II Incipient\n(sound, confidently wet)", cex = 0.7,
       col = CAT_COLS[["II"]], font = 2, adj = c(1, 0.5))
  text(0, sqrt(34) + 0.25,
       sprintf("dead band ±%.1f SD\n(transition → defaults to row baseline)", cut_),
       cex = 0.65, col = "#555555")
  title(main = sprintf("Dead band ±%.1f SD", cut_), cex.main = 0.9)
  if (first) {
    legend("left", legend = names(CAT_COLS), pch = 21, pt.bg = unname(CAT_COLS),
           col = "white", pt.cex = 1.25, cex = 0.75, bty = "n", title = "Category")
  }
  first <- FALSE
}
mtext(paste0("6-cell scheme: 4 confident corners + a \"normal\" transition column that defaults to the row baseline\n",
             "(sound→No-decay, damaged→Active). Only \"confidently wet sound\" (II) and \"confidently dry damaged\" (IV) are off-baseline."),
      outer = TRUE, cex = 0.85, line = 0.6)
invisible(dev.off())
cat("\nsaved 6-cell figure\n")
