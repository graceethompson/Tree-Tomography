# R port of threshold.py (revision sensitivity analysis, cjfr-2026-0202).
# Run from the repo root: Rscript analysis/revision/scripts/threshold.R
source("analysis/revision/scripts/revision_common.R")

# ---- small formatting helpers (match the Python printout) ----
pydict <- function(keys, vals) {
  paste0("{", paste(sprintf("'%s': %s", keys, vals), collapse = ", "), "}")
}
fmtnum <- function(x) vapply(x, function(v) format(v), character(1))
bar66 <- strrep("=", 66)

# scipy defaults: biased (population) skewness; Fisher, biased kurtosis
skewness <- function(x) { mu <- mean(x); s <- sqrt(mean((x - mu)^2)); mean((x - mu)^3) / s^3 }
kurtosis_fisher <- function(x) { mu <- mean(x); v <- mean((x - mu)^2); mean((x - mu)^4) / v^2 - 3 }

m <- load_merged()

for (scheme in c("species", "pooled")) {
  b <- build_pc1(m, scheme)
  m[[paste0("pc1_", scheme)]] <- b$pc1
  cat(sprintf("\n=== %s normalization: PC1 %.1f%%, PC2 %.1f%% ===\n",
              scheme, b$ve[1] * 100, b$ve[2] * 100))
  cat("  loadings:", pydict(ERT_METRICS, fmtnum(round(b$loadings[ERT_METRICS], 2))), "\n")
}

cat("\ncorr(species-PC1, pooled-PC1) =",
    format(round(cor(m$pc1_species, m$pc1_pooled), 3)), "\n")

# ---- MODALITY: is there a natural break? ----
cat("\n", bar66, " \nMODALITY of the anomaly axis (natural break?)\n ", bar66, "\n", sep = "")
modality <- function(x, label) {
  g <- gmm_1d(x)
  verdict <- if (g$bic2 + 2 < g$bic1) "2-comp (bimodal) favored" else "unimodal favored"
  cat(sprintf("%s: BIC 1-comp=%.1f  2-comp=%.1f  -> %s\n", label, g$bic1, g$bic2, verdict))
  cat(sprintf("   skew=%.2f kurtosis=%.2f\n", skewness(x), kurtosis_fisher(x)))
  g
}
g2_sp <- modality(m$pc1_species, "species-PC1")
g2_pl <- modality(m$pc1_pooled, "pooled-PC1")

# ---- Threshold rules on species-PC1 (primary axis) ----
pc <- m$pc1_species
mu <- mean(pc); med <- median(pc); sdv <- sd(pc)
t_otsu <- otsu_threshold(pc)
t_gmm <- g2_sp$crossover
rules <- list("mean (published)" = mu, "median" = med,
              "GMM 2-comp crossover" = t_gmm, "Otsu break" = t_otsu,
              "mean-0.5SD" = mu - 0.5 * sdv, "mean+0.5SD" = mu + 0.5 * sdv,
              "mean+1SD" = mu + 1 * sdv)
m$struct <- m$percent_damaged > 1
cat("\n", bar66, " \nTHRESHOLD SENSITIVITY on species-PC1\n ", bar66, "\n", sep = "")
cat(sprintf("%-22s%7s%7s%7s%9s%9s  species incipient %%\n",
            "rule", "thr", "nAnom", "nIncip", "BGSinc%", "EMSinc%"))
spp <- sort(unique(m$sp))
for (name in names(rules)) {
  t <- rules[[name]]
  anom <- pc > t
  incip <- (!m$struct) & anom
  bgs <- sum(incip[m$site == "BGS"]) / sum(m$site == "BGS") * 100
  ems <- sum(incip[m$site == "EMS"]) / sum(m$site == "EMS") * 100
  sppc <- vapply(spp, function(s) round(sum(incip[m$sp == s]) / sum(m$sp == s) * 100),
                 numeric(1))
  cat(sprintf("%-22s%7.2f%7d%7d%9.0f%9.0f  %s\n", name, t, sum(anom), sum(incip),
              bgs, ems, pydict(spp, sprintf("%d", sppc))))
}

# how many trees flip anomaly label across the full rule set
labels <- sapply(rules, function(t) as.integer(pc > t))   # 57 x 7
lab_min <- apply(labels, 1, min); lab_max <- apply(labels, 1, max)
cat(sprintf("\nTrees whose anomaly label is NOT stable across all rules: %d/%d\n",
            sum(lab_min != lab_max), nrow(m)))
cat(sprintf("Always-anomalous: %d, never-anomalous: %d\n",
            sum(lab_min == 1), sum(lab_max == 0)))

# ---- Within vs between: site contrast under each normalization at its own mean ----
cat("\n", bar66,
    " \nWITHIN (species) vs BETWEEN (pooled) normalization: site incipient contrast\n ",
    bar66, "\n", sep = "")
for (scheme in c("species", "pooled")) {
  a <- m[[paste0("pc1_", scheme)]]
  t <- mean(a)
  incip <- (!m$struct) & (a > t)
  bgs <- mean(incip[m$site == "BGS"]) * 100
  ems <- mean(incip[m$site == "EMS"]) * 100
  sp_means <- vapply(spp, function(s) round(mean(a[m$sp == s]), 2), numeric(1))
  cat(sprintf("%8s: BGS incip %.0f%% vs EMS %.0f%%  | per-species mean PC1 %s\n",
              scheme, bgs, ems, pydict(spp, fmtnum(sp_means))))
}

# ---- External anchor demonstration via validation regression ----
cat("\n", bar66,
    " \nEXTERNAL ANCHOR demo (validation moisture <-> mean resistivity)\n ",
    bar66, "\n", sep = "")
v <- load_validation()
fit <- lm(moisture ~ mean_res, data = v)
itc <- unname(coef(fit)[1]); slp <- unname(coef(fit)[2])
rv <- cor(v$mean_res, v$moisture)
cat(sprintf("validation: moisture = %.1f + %.3f*mean_res  (r=%.2f)\n", itc, slp, rv))
# invert: resistivity at a moisture anchor
for (Mcut in c(90, 100, 110)) {
  Rcut <- (Mcut - itc) / slp
  n_below <- sum(m$mean < Rcut)
  cat(sprintf("  moisture>%d%% <=> mean_res<%.0f Ohm-m -> %d/%d main trees flagged (species-blind, absolute)\n",
              Mcut, Rcut, n_below, nrow(m)))
}

# ---------- FIGURE ----------
# matplotlib figsize (12, 4.6) at dpi 170 -> 2040 x 782 px
png(file.path(OUT_DIR, "CJFR-threshold-rules.png"),
    width = 2040, height = 782, res = 170, pointsize = 10)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.4, 0.7, 0))
rule_cols <- c("mean (published)" = "#c0392b", "median" = "#2980b9",
               "GMM 2-comp crossover" = "#27ae60", "Otsu break" = "#8e44ad")
panels <- list(
  list(col = "pc1_species", title = "Within-species normalization (published)"),
  list(col = "pc1_pooled",  title = "Pooled normalization (preserves absolute differences)"))
for (p in panels) {
  a <- m[[p$col]]
  brk <- seq(min(a), max(a), length.out = 17)      # 16 bins over the data range
  h <- hist(a, breaks = brk, plot = FALSE)
  kd <- density(a)
  xs <- seq(min(a) - 0.3, max(a) + 0.3, length.out = 300)
  ky <- approx(kd$x, kd$y, xs, rule = 2)$y
  ylim <- c(0, max(h$density, ky) * 1.06)
  plot(NA, xlim = range(xs), ylim = ylim, xaxs = "i", yaxs = "i",
       xlab = "ERT PC1 (higher = lower resistivity / more heterogeneous)",
       ylab = "density", main = "")
  plot(h, freq = FALSE, col = "grey85", border = "white", add = TRUE)
  lines(xs, ky, col = "black", lwd = 1.8)
  if (p$col == "pc1_species") {
    labs <- character(0); lcols <- character(0)
    for (name in names(rule_cols)) {
      t <- rules[[name]]
      if (is.finite(t)) {
        abline(v = t, col = rule_cols[[name]], lty = 2, lwd = 1.6)
        labs <- c(labs, sprintf("%s=%.2f", name, t))
        lcols <- c(lcols, rule_cols[[name]])
      }
    }
    legend("topright", legend = labs, col = lcols, lty = 2, lwd = 1.6,
           bty = "n", cex = 0.7)
  } else {
    abline(v = mean(a), col = "#c0392b", lty = 2, lwd = 1.6)
    legend("topright", legend = sprintf("mean=%.2f", mean(a)), col = "#c0392b",
           lty = 2, lwd = 1.6, bty = "n", cex = 0.7)
  }
  title(main = p$title, cex.main = 1, font.main = 1)
}
invisible(dev.off())
cat("\nsaved figure\n")
