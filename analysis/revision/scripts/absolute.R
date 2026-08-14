# R port of absolute.py (revision sensitivity analysis, cjfr-2026-0202).
# Run from the repo root: Rscript analysis/revision/scripts/absolute.R
source("analysis/revision/scripts/revision_common.R")

skewness <- function(x) { mu <- mean(x); s <- sqrt(mean((x - mu)^2)); mean((x - mu)^3) / s^3 }

m <- load_merged()
m$res <- m$mean               # mean resistivity, Ohm-m (absolute)
m$struct <- m$percent_damaged > 1
ns <- !m$struct
mb <- m$site == "BGS"; me <- m$site == "EMS"

# ---- validation calibration: moisture ~ resistivity (n=12) ----
v <- load_validation()
v$res <- v$mean_res
fit <- lm(moisture ~ res, data = v)
itc <- unname(coef(fit)[1]); slp <- unname(coef(fit)[2])
rv <- cor(v$res, v$moisture)
n <- nrow(v); xbar <- mean(v$res); Sxx <- sum((v$res - xbar)^2)
resid <- v$moisture - (itc + slp * v$res)
rstd <- sqrt(sum(resid^2) / (n - 2))
pred_moist <- function(x) itc + slp * x
pred_se <- function(x) rstd * sqrt(1 + 1 / n + (x - xbar)^2 / Sxx)   # prediction SE
res_for_moist <- function(Mcut) {  # invert: resistivity giving predicted moisture=Mcut, with band
  x <- (Mcut - itc) / slp
  band <- pred_se(x) / abs(slp)
  c(x, x - 1.96 * band, x + 1.96 * band)
}
m$pmoist <- pred_moist(m$res)
r100 <- res_for_moist(100)[1]
cat(sprintf("Calibration: moisture(%%) = %.1f %+.4f*res   r=%.2f, n=%d, resid s=%.1f%%\n",
            itc, slp, rv, n, rstd))
cat(sprintf("Predicted moisture range across 57 trees: %.0f-%.0f%%\n",
            min(m$pmoist), max(m$pmoist)))

cat("\n--- Absolute mean resistivity (Ohm-m) by species & site ---\n")
spp <- sort(unique(m$sp))
cat(sprintf("%-13s%6s%7s%7s%7s\n", "", "count", "mean", "min", "max"))
cat(sprintf("%-13s\n", "sp"))
for (s in spp) {
  d <- m$res[m$sp == s]
  cat(sprintf("%-13s%6d%7.1f%7.1f%7.1f\n", s, length(d),
              round(mean(d)), round(min(d)), round(max(d))))
}
cat(sprintf("{'mean': {'BGS': %.1f, 'EMS': %.1f}, 'median': {'BGS': %.1f, 'EMS': %.1f}}\n",
            round(mean(m$res[mb])), round(mean(m$res[me])),
            round(median(m$res[mb])), round(median(m$res[me]))))

# ---- modality on absolute resistivity ----
modality <- function(x, label) {
  g <- gmm_1d(x)
  verdict <- if (g$bic2 + 2 < g$bic1) "bimodal" else "unimodal"
  cat(sprintf("%s: BIC 1=%.1f 2=%.1f -> %s; skew=%.2f\n",
              label, g$bic1, g$bic2, verdict, skewness(x)))
  invisible(g)
}
modality(m$res, "abs resistivity")
modality(m$pmoist, "pred moisture")

# ---- candidate ABSOLUTE thresholds (anomaly = LOW resistivity / HIGH moisture) ----
res <- m$res
t_otsu <- otsu_threshold(res)
t_gmm <- gmm_1d(res)$crossover
rules <- list(
  list(name = "Sample median resistivity",         t = median(res)),
  list(name = "Sample mean resistivity",           t = mean(res)),
  list(name = "Otsu break (data-driven)",          t = t_otsu),
  list(name = "GMM crossover (data-driven)",       t = t_gmm),
  list(name = "Lit: silver-fir wetwood ~200",      t = 200.0),
  list(name = "Lit: silver-fir low-ER ~166",       t = 166.0),
  list(name = "Validation: moisture>90% (<505)",   t = res_for_moist(90)[1]),
  list(name = "Validation: moisture>100% (<402)",  t = res_for_moist(100)[1]),
  list(name = "Validation: moisture>110% (<299)",  t = res_for_moist(110)[1]))
bar92 <- strrep("=", 92)
cat("\n", bar92, "\n", sep = "")
cat(sprintf("%-34s%7s%9s%7s%7s%9s%9s\n",
            "ABSOLUTE threshold rule", "res<", "~moist>", "nAnom", "nIncip", "BGSinc%", "EMSinc%"))
cat(bar92, "\n", sep = "")
tbl <- list()
for (r in rules) {
  anom <- res < r$t
  inc <- ns & anom
  bgsr <- sum(inc[mb]) / sum(mb) * 100; emsr <- sum(inc[me]) / sum(me) * 100
  mo <- pred_moist(r$t)
  cat(sprintf("%-34s%7.0f%9.0f%7d%7d%9.0f%9.0f\n",
              r$name, r$t, mo, sum(anom), sum(inc), bgsr, emsr))
  tbl[[length(tbl) + 1]] <- data.frame(
    rule = r$name, res_thresh = as.integer(round(r$t)),
    moist_equiv = as.integer(round(mo)), n_anom = sum(anom), n_incip = sum(inc),
    BGS_inc = as.integer(round(bgsr)), EMS_inc = as.integer(round(emsr)))
}
write.csv(do.call(rbind, tbl),
          file.path(OUT_DIR, "CJFR-absolute-thresholds.csv"),
          row.names = FALSE, quote = FALSE)

# stability of anomaly label across absolute rule set
lab <- sapply(rules, function(r) as.integer(res < r$t))   # 57 x 9
lab_min <- apply(lab, 1, min); lab_max <- apply(lab, 1, max)
cat(sprintf("\nAlways-anomalous (all abs rules): %d, never: %d, unstable: %d\n",
            sum(lab_min == 1), sum(lab_max == 0), sum(lab_min != lab_max)))

# does site DIRECTION hold on absolute axis?
grid <- seq(quantile(res, 0.05, type = 7), quantile(res, 0.95, type = 7), length.out = 200)
bgs <- numeric(length(grid)); ems <- numeric(length(grid))
for (i in seq_along(grid)) {
  inc <- ns & (res < grid[i])
  bgs[i] <- sum(inc[mb]) / sum(mb) * 100
  ems[i] <- sum(inc[me]) / sum(me) * 100
}
cat(sprintf("\nSite DIRECTION on absolute axis: BGS>=EMS at %.0f%% of thresholds; BGS>EMS at %.0f%%; crossovers where EMS>BGS at %.0f%%\n",
            mean(bgs >= ems) * 100, mean(bgs > ems) * 100, mean(ems > bgs) * 100))
cat(sprintf("BGS incip range %.0f-%.0f%%, EMS %.0f-%.0f%%\n",
            min(bgs), max(bgs), min(ems), max(ems)))

# ---------------- FIGURE (3 panels) ----------------
# matplotlib figsize (16.5, 4.9) at dpi 170 -> 2805 x 833 px
png(file.path(OUT_DIR, "CJFR-absolute-fig.png"),
    width = 2805, height = 833, res = 170, pointsize = 10)
par(mfrow = c(1, 3), mgp = c(2.4, 0.7, 0))
par(cex = 1)   # undo the automatic mfrow cex shrink so text matches matplotlib
OHM <- intToUtf8(0x03A9); MIDDOT <- intToUtf8(0x00B7)
EMDASH <- intToUtf8(0x2014); ELLIP <- intToUtf8(0x2026)
ohm_m <- paste0(OHM, MIDDOT, "m")
spp_order <- c("A.rubrum", "T.canadensis", "N.sylvatica", "Q.rubra")

# A: absolute resistivity by species (the confound made visible)
set.seed(0)
par(mar = c(4, 6.5, 3.5, 1))
plot(NA, xlim = rev(range(m$res) + c(-25, 25)), ylim = c(-0.5, 3.5),
     xlab = paste0("mean resistivity (", ohm_m, ")  ", EMDASH, " lower = wetter"),
     ylab = "", yaxt = "n", main = "")
for (i in seq_along(spp_order)) {
  s <- spp_order[i]
  d <- m$res[m$sp == s]
  points(d, rnorm(length(d), i - 1, 0.08), col = "white",
         bg = adjustcolor(SPECIES_COLS[[s]], alpha.f = 0.8), pch = 21, cex = 1.15, lwd = 0.5)
  segments(mean(d), i - 1 - 0.28, mean(d), i - 1 + 0.28, col = SPECIES_COLS[[s]], lwd = 2.5)
}
axis(2, at = 0:3, las = 1, cex.axis = 0.95,
     labels = paste0(spp_order, ifelse(spp_order %in% c("N.sylvatica", "Q.rubra"), "*", "")))
title(main = "A.  Absolute axis is species-confounded\n(* = single-site specialist)",
      adj = 0, cex.main = 1.05, font.main = 1)

# B: predicted-moisture distribution by site with candidate thresholds
par(mar = c(4, 4, 3.5, 1))
xs2 <- seq(min(m$pmoist) - 5, max(m$pmoist) + 5, length.out = 300)
kb <- approx(density(m$pmoist[mb])$x, density(m$pmoist[mb])$y, xs2, rule = 2)$y
ke <- approx(density(m$pmoist[me])$x, density(m$pmoist[me])$y, xs2, rule = 2)$y
ymax <- max(kb, ke) * 1.05
plot(NA, xlim = range(xs2), ylim = c(-0.0022, ymax),
     xlab = "predicted moisture content (%)  [from validation calibration]",
     ylab = "density", main = "")
lines(xs2, kb, col = "#1f77b4", lwd = 2)
lines(xs2, ke, col = "#d9822b", lwd = 2)
points(m$pmoist[mb], rep(-0.001, sum(mb)), col = adjustcolor("#1f77b4", 0.6), pch = 16, cex = 0.55)
points(m$pmoist[me], rep(-0.0006, sum(me)), col = adjustcolor("#d9822b", 0.6), pch = 16, cex = 0.55)
for (mc in list(c(90, "#999999"), c(100, "#555555"), c(110, "#111111"))) {
  x0 <- as.numeric(mc[[1]])
  abline(v = x0, col = mc[[2]], lty = 3, lwd = 1.2)
  text(x0, ymax * 0.92, sprintf("%d%%", x0), cex = 0.6)
}
legend("topright", legend = c("BGS", "EMS"), col = c("#1f77b4", "#d9822b"),
       lwd = 2, bty = "n", cex = 0.9)
title(main = "B.  Absolute (predicted moisture) axis by site",
      adj = 0, cex.main = 1.05, font.main = 1)

# C: site direction sweep on absolute resistivity, anchors marked
par(mar = c(4, 4, 3.5, 1))
ymax2 <- max(bgs, ems) * 1.04
plot(NA, xlim = rev(range(grid)), ylim = c(0, ymax2),
     xlab = paste0("anomaly threshold: resistivity below ", ELLIP, " (", ohm_m, ")"),
     ylab = "% of site trees \"incipient\"", main = "")
kp <- bgs >= ems
if (any(kp)) polygon(c(grid[kp], rev(grid[kp])), c(bgs[kp], rev(ems[kp])),
                     col = adjustcolor("#1f77b4", alpha.f = 0.10), border = NA)
ko <- ems > bgs
if (any(ko)) polygon(c(grid[ko], rev(grid[ko])), c(bgs[ko], rev(ems[ko])),
                     col = adjustcolor("#d9822b", alpha.f = 0.18), border = NA)
lines(grid, bgs, col = "#1f77b4", lwd = 2.4)
lines(grid, ems, col = "#d9822b", lwd = 2.4)
anchor <- c(median = median(res), mean = mean(res), Otsu = t_otsu,
            "moist>100%" = r100, "wetwood~200" = 200)
for (nm in names(anchor)) {
  abline(v = anchor[[nm]], col = "#444444", lty = 3, lwd = 1.0)
  text(anchor[[nm]], ymax2 * 0.98, nm, srt = 90, adj = c(1, 0), cex = 0.62, col = "#444444")
}
legend("topleft", legend = c("BGS (wetland)", "EMS (upland)"),
       col = c("#1f77b4", "#d9822b"), lwd = 2.4, bty = "n", cex = 0.9)
title(main = "C.  Site direction across absolute thresholds\n(shaded orange = where upland exceeds wetland)",
      adj = 0, cex.main = 1.05, font.main = 1)
invisible(dev.off())
cat("\nsaved figure + CSV\n")
