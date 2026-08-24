# R port of threshold_fig.py (revision sensitivity analysis, cjfr-2026-0202).
# Run from the repo root: Rscript analysis/revision/scripts/threshold_fig.R
source("analysis/revision/scripts/revision_common.R")

m <- load_merged()
b <- build_pc1(m, "species")
pc <- b$pc1
pc_pooled <- build_pc1(m, "pooled")$pc1
m$struct <- m$percent_damaged > 1

mu <- mean(pc); med <- median(pc); sdv <- sd(pc)
t_otsu <- otsu_threshold(pc)
t_gmm <- gmm_1d(pc)$crossover
rules <- list("mean (published)" = list(t = mu,     col = "#c0392b"),
              "median"           = list(t = med,    col = "#2980b9"),
              "Otsu break"       = list(t = t_otsu, col = "#8e44ad"),
              "GMM crossover"    = list(t = t_gmm,  col = "#27ae60"))

# Panel B sweep: 5th-95th percentile (numpy default = R quantile type 7)
ts <- seq(quantile(pc, 0.05, type = 7), quantile(pc, 0.95, type = 7), length.out = 200)
mask_bgs <- m$site == "BGS"; mask_ems <- m$site == "EMS"; ns <- !m$struct
bgs <- numeric(length(ts)); ems <- numeric(length(ts))
for (i in seq_along(ts)) {
  inc <- ns & (pc > ts[i])
  bgs[i] <- sum(inc[mask_bgs]) / sum(mask_bgs) * 100
  ems[i] <- sum(inc[mask_ems]) / sum(mask_ems) * 100
}

# ---------- FIGURE ----------
png(file.path(OUT_DIR, "CJFR-threshold-fig.png"),
    width = 3282, height = 840, res = 175, pointsize = 10)
par(mfrow = c(1, 3), mar = c(4, 4, 3.5, 1), mgp = c(2.4, 0.7, 0))

# Panel A: density with candidate thresholds
xs2 <- seq(min(pc) - 0.4, max(pc) + 0.4, length.out = 300)
kd <- density(pc)
ky <- approx(kd$x, kd$y, xs2, rule = 2)$y
brk <- seq(min(pc), max(pc), length.out = 16)     # 15 bins over the data range
h <- hist(pc, breaks = brk, plot = FALSE)
plot(NA, xlim = range(xs2), ylim = c(0, 0.42), xaxs = "i", yaxs = "i",
     xlab = "ERT PC1  (higher = lower resistivity / more heterogeneous)",
     ylab = "density", main = "")
plot(h, freq = FALSE, col = "grey86", border = "white", add = TRUE)
lines(xs2, ky, col = "black", lwd = 2)
labs <- character(0); lcols <- character(0)
for (name in names(rules)) {
  r <- rules[[name]]
  if (is.finite(r$t)) {
    abline(v = r$t, col = r$col, lty = 2, lwd = 1.7)
    labs <- c(labs, sprintf("%s = %.2f", name, r$t))
    lcols <- c(lcols, r$col)
  }
}
legend("topright", legend = labs, col = lcols, lty = 2, lwd = 1.7, bty = "n", cex = 0.8)
EMDASH <- intToUtf8(0x2014)
title(main = "A", adj = 0, cex.main = 1, font.main = 2)

# Panel B: threshold sweep, BGS vs EMS incipient %
plot(NA, xlim = range(ts), ylim = range(c(bgs, ems)) + c(0, 4),
     xlab = "anomaly threshold on ERT PC1",
     ylab = "% of site trees classified \"incipient\"", main = "")
keep <- bgs >= ems
polygon(c(ts[keep], rev(ts[keep])), c(bgs[keep], rev(ems[keep])),
        col = adjustcolor("#1f77b4", alpha.f = 0.10), border = NA)
lines(ts, bgs, col = "#1f77b4", lwd = 2.4)
lines(ts, ems, col = "#d9822b", lwd = 2.4)
for (name in names(rules)) {
  r <- rules[[name]]
  if (is.finite(r$t)) abline(v = r$t, col = r$col, lty = 3, lwd = 1.3)
}
legend("topright", legend = c("BGS (wetland)", "EMS (upland)"),
       col = c("#1f77b4", "#d9822b"), lwd = 2.4, bty = "n", cex = 0.9)
title(main = "B", adj = 0, cex.main = 1, font.main = 2)

# Panel C: same distribution under the pooled normalization (preserves
# absolute between-site/species differences); replaces the former
# standalone normalization figure (CJFR-threshold-rules.png, archived).
kd3 <- density(pc_pooled)
xs3 <- seq(min(pc_pooled) - 0.4, max(pc_pooled) + 0.4, length.out = 300)
ky3 <- approx(kd3$x, kd3$y, xs3, rule = 2)$y
brk3 <- seq(min(pc_pooled), max(pc_pooled), length.out = 16)
h3 <- hist(pc_pooled, breaks = brk3, plot = FALSE)
plot(NA, xlim = range(xs3), ylim = c(0, max(h3$density, ky3) * 1.08),
     xaxs = "i", yaxs = "i",
     xlab = "ERT PC1, pooled normalization",
     ylab = "density", main = "")
plot(h3, freq = FALSE, col = "grey86", border = "white", add = TRUE)
lines(xs3, ky3, col = "black", lwd = 2)
abline(v = mean(pc_pooled), col = "#c0392b", lty = 2, lwd = 1.7)
legend("topright", legend = sprintf("mean = %.2f", mean(pc_pooled)),
       col = "#c0392b", lty = 2, lwd = 1.7, bty = "n", cex = 0.8)
title(main = "C", adj = 0, cex.main = 1, font.main = 2)
invisible(dev.off())

cat("saved. BGS never below EMS across sweep:",
    ifelse(all(bgs >= ems), "True", "False"), "\n")
cat(sprintf("BGS incipient range across sweep: %.0f-%.0f%% ; EMS: %.0f-%.0f%%\n",
            min(bgs), max(bgs), min(ems), max(ems)))
