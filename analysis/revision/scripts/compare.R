# Does the ERT axis choice change the binning? PC1 vs CMA vs absolute /
# species-normalized resistivity. R port of compare.py. Run from the repo root:
#   Rscript analysis/revision/scripts/compare.R

source("analysis/revision/scripts/revision_common.R")

# UTF-8 locale so the png devices render em-dash / sqrt / arrow glyphs
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"))) invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))

m <- load_merged()
pca <- build_pc1(m, "species")
pc <- pca$pc1
Z <- pca$Z
mets <- ERT_METRICS
m$pc1 <- pc
m$res_z <- Z$mean
m$struct <- m$percent_damaged > 1
mu <- mean(pc)

# ---- PC1 loadings: does PC1 "contain" CMA? ----
cat("PC1 loadings (does PC1 include CMA?):\n")
ld <- pca$loadings
for (met in names(ld)[order(-abs(ld))]) {
  cat(sprintf("   %-14s%+.2f\n", met, ld[[met]]))
}
cat(sprintf("   -> CMA correlation with PC1: r=%+.2f   (CMA loads on PC2, not PC1)\n",
            cor(m$cma, pc)))
cat(sprintf("   -> mean resistivity vs PC1: r=%+.2f\n", cor(m$mean, pc)))

four <- function(anom, struct) four_cat(anom, struct)
# Scheme A: original PC1 (anomaly = wet/heterogeneous = PC1>mean)
A <- four(pc > mu, m$struct)
# Scheme B: CMA (anomaly = central moisture accumulation > 0.33, the paper's own 'even' point)
B <- four(m$cma > 0.33, m$struct)
# Scheme C: absolute mean resistivity (anomaly = wetter than sample median)
C <- four(m$mean < median(m$mean), m$struct)
# Scheme D: species-normalized resistivity only (anomaly = res_z<0)
D <- four(m$res_z < 0, m$struct)

cat_counts <- function(cc) {
  t <- table(factor(cc, levels = c("I", "II", "III", "IV")))
  paste(sprintf("%s=%d", names(t), as.integer(t)), collapse = " ")
}
schemes <- list("A PC1 (original)" = A, "B CMA>0.33" = B,
                "C abs-resistivity" = C, "D species-resistivity" = D)
for (nm in names(schemes)) cat(sprintf("\n%s: {%s}\n", nm, cat_counts(schemes[[nm]])))

agree <- function(x, y) mean(x == y) * 100
cat("\n--- AGREEMENT with original PC1 scheme (A) ---\n")
for (nm in c("B CMA", "C abs-res", "D sp-res")) {
  arr <- list("B CMA" = B, "C abs-res" = C, "D sp-res" = D)[[nm]]
  cat(sprintf("  A vs %s: %.0f%% identical  (%d of 57 trees change category)\n",
              nm, agree(A, arr), sum(A != arr)))
}
cat("  A vs B, restricted to the moisture-axis calls only (I/II among sound; III/IV among damaged):\n")
sound <- !m$struct; dmgd <- m$struct
cat(sprintf("     I/II (sound):  %.0f%% (%d/%d differ)\n",
            agree(A[sound], B[sound]), sum(A[sound] != B[sound]), sum(sound)))
cat(sprintf("     III/IV (damaged): %.0f%% (%d/%d differ)\n",
            agree(A[dmgd], B[dmgd]), sum(A[dmgd] != B[dmgd]), sum(dmgd)))

# confusion A vs B
cat("\nConfusion A (rows) vs B (cols):\n")
lv <- c("I", "II", "III", "IV")
print(table(PC1 = factor(A, levels = lv), CMA = factor(B, levels = lv)))

# ---------------- FIGURE: phase diagram under 3 binnings ----------------
catcol <- c(I = "#3b6fb0", II = "#5aa02c", III = "#d98a1f", IV = "#b83232")
yt <- c(0, 1, 5, 10, 20, 30)
y <- sqrt(pmax(m$percent_damaged, 0))

panel <- function(x, anom_line, xlabel, cats, title, xinvert = FALSE, ylab = FALSE) {
  xr <- range(x); xpad <- 0.05 * diff(xr)
  xlim <- c(xr[1] - xpad, xr[2] + xpad)
  if (xinvert) xlim <- rev(xlim)
  plot(NA, xlim = xlim, ylim = c(-0.2, sqrt(36)), axes = FALSE, xlab = "", ylab = "")
  abline(h = sqrt(1), col = "#888888", lty = 2, lwd = 1)
  abline(v = anom_line, col = "#888888", lty = 2, lwd = 1)
  for (cat_ in names(catcol)) {
    idx <- cats == cat_
    points(x[idx], y[idx] + rnorm(sum(idx), 0, 0.03),
           pch = 21, bg = catcol[[cat_]], col = "white", lwd = 0.5, cex = 1.2)
  }
  axis(1, cex.axis = 0.8)
  axis(2, at = sqrt(yt), labels = yt, cex.axis = 0.8, las = 1)
  box()
  title(xlab = xlabel, cex.lab = 0.8)
  if (ylab) title(ylab = "SoT structural loss %  (√ scale)", cex.lab = 0.85)
  title(main = title, cex.main = 0.85)
}

set.seed(1)
png(file.path(OUT_DIR, "CJFR-binning-compare.png"),
    width = 16, height = 5.2, units = "in", res = 150)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.4, 1), oma = c(0, 0, 2.2, 0), mgp = c(2.4, 0.7, 0))
panel(pc, mu, "ERT PC1 (species-normalized composite)", A, "A. Original — bin by PC1", ylab = TRUE)
legend("topleft", legend = names(catcol), pch = 21, pt.bg = unname(catcol),
       col = "white", pt.cex = 1.2, cex = 0.8, bty = "n", title = "Category")
panel(m$cma, 0.33, "CMA — central moisture accumulation", B, "B. Bin by CMA (>0.33)")
panel(m$mean, median(m$mean), "mean resistivity (Ω·m, absolute)", C,
      "C. Bin by absolute resistivity", xinvert = TRUE)
mtext("Does the ERT axis choice change the binning? Same trees, three anomaly criteria",
      outer = TRUE, cex = 0.9, font = 1, line = 0.6)
invisible(dev.off())

# second fig: PC1 diagram colored by CMA value (continuous) to show CMA is a SEPARATE axis
png(file.path(OUT_DIR, "CJFR-pc1-cma.png"),
    width = 7.6, height = 5.6, units = "in", res = 150)
layout(matrix(1:2, 1, 2), widths = c(6.4, 1.2))
par(mar = c(4.2, 4.2, 3.4, 0.6), mgp = c(2.4, 0.7, 0))
pal <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(256)
vv <- pmin(pmax(m$cma, 0), 0.5) / 0.5
ptcol <- pal[pmax(1, ceiling(vv * 256))]
xr <- range(pc); xpad <- 0.05 * diff(xr)
plot(NA, xlim = c(xr[1] - xpad, xr[2] + xpad), ylim = c(-0.2, sqrt(36)),
     axes = FALSE, xlab = "", ylab = "")
abline(h = sqrt(1), col = "#888888", lty = 2, lwd = 1)
abline(v = mu, col = "#888888", lty = 2, lwd = 1)
points(pc, y + rnorm(length(y), 0, 0.03), pch = 21, bg = ptcol, col = "black",
       lwd = 0.4, cex = 1.4)
axis(1, cex.axis = 0.8)
axis(2, at = sqrt(yt), labels = yt, cex.axis = 0.8, las = 1)
box()
title(xlab = "ERT PC1", cex.lab = 0.9)
title(ylab = "SoT structural loss % (√)", cex.lab = 0.9)
title(main = "PC1 phase diagram colored by CMA:\nCMA does NOT line up with the PC1 axis (r=+0.05) — it is orthogonal (PC2)",
      cex.main = 0.8)
# colorbar
par(mar = c(4.2, 0.8, 3.4, 2.6))
zz <- matrix(seq(0, 0.5, length.out = 256), nrow = 1)
image(x = 1, y = seq(0, 0.5, length.out = 256), z = zz, col = pal,
      axes = FALSE, xlab = "", ylab = "")
axis(4, cex.axis = 0.8, las = 1)
mtext("CMA (blue=wet core, red=dry core)", side = 4, line = 1.7, cex = 0.75)
box()
invisible(dev.off())
cat("\nsaved figures\n")
