# Statistical tests supporting the revision (R port of statistical_tests.py):
#  1. Fisher exact tests on the incipient site contrast (by threshold; generalists-only)
#  2. Calibration-anchor uncertainty: CI on the mean calibration line vs single-tree PI
# Run from repo root: Rscript analysis/revision/scripts/statistical_tests.R
source("analysis/revision/scripts/revision_common.R")

m <- load_merged()
bp <- build_pc1(m, "species")
pc <- bp$pc1
m$pc1 <- pc
m$struct <- m$percent_damaged > 1
sdv <- sd(pc)          # ddof=1, as in the Python script
sound <- !m$struct

# Sample odds ratio (a*d)/(b*c), matching scipy.stats.fisher_exact's OR;
# p-value from R's fisher.test (identical two-sided p).
fisher_line <- function(a, b, c, d) {
  or <- (a * d) / (b * c)
  p <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value
  list(or = if (is.finite(or)) sprintf("%.2f", or) else "inf", p = p)
}

cat("1. Fisher exact: incipient site contrast by ERT threshold\n")
cases <- list(list(-0.5, "mean-0.5SD"), list(0, "mean (published)"),
              list(0.5, "mean+0.5SD"), list(1, "mean+1SD"))
for (cs in cases) {
  t0 <- mean(pc) + cs[[1]] * sdv
  inc <- sound & (pc > t0)
  a <- sum(inc & m$site == "BGS"); b <- sum(sound & m$site == "BGS" & !inc)
  c <- sum(inc & m$site == "EMS"); d <- sum(sound & m$site == "EMS" & !inc)
  f <- fisher_line(a, b, c, d)
  cat(sprintf("  %-18s: BGS %d/%d vs EMS %d/%d   OR=%s  p=%.3f\n",
              cs[[2]], a, a + b, c, c + d, f$or, f$p))
}

inc0 <- sound & (pc > mean(pc))
a <- sum(inc0 & m$site == "BGS"); n1 <- sum(m$site == "BGS")
c <- sum(inc0 & m$site == "EMS"); n2 <- sum(m$site == "EMS")
f <- fisher_line(a, n1 - a, c, n2 - c)
cat(sprintf("  As %% of ALL trees (31%% vs 11%%): %d/%d vs %d/%d  OR=%s  p=%.3f\n",
            a, n1, c, n2, f$or, f$p))

g <- m[m$sp %in% c("A.rubrum", "T.canadensis"), ]
incg <- (!g$struct) & (g$pc1 > mean(pc))
a <- sum(incg & g$site == "BGS"); n1 <- sum(g$site == "BGS")
c <- sum(incg & g$site == "EMS"); n2 <- sum(g$site == "EMS")
f <- fisher_line(a, n1 - a, c, n2 - c)
cat(sprintf("  GENERALISTS ONLY: BGS %d/%d (%.0f%%) vs EMS %d/%d (%.0f%%)  OR=%s  p=%.3f\n",
            a, n1, a / n1 * 100, c, n2, c / n2 * 100, f$or, f$p))

cat("\n2. Absolute anchor uncertainty (hemlock calibration, n=12)\n")
v <- load_validation()
res <- v$mean_res
fit <- lm(moisture ~ res, data = data.frame(res = res, moisture = v$moisture))
intercept <- coef(fit)[1]; slope <- coef(fit)[2]
n <- nrow(v); xbar <- mean(res); Sxx <- sum((res - xbar)^2)
rstd <- sqrt(sum(residuals(fit)^2) / (n - 2))
tcrit <- qt(0.975, n - 2)
for (M in c(90, 100, 110)) {
  x <- (M - intercept) / slope
  ci <- rstd * sqrt(1 / n + (x - xbar)^2 / Sxx) / abs(slope) * tcrit
  pi <- rstd * sqrt(1 + 1 / n + (x - xbar)^2 / Sxx) / abs(slope) * tcrit
  cat(sprintf("  moisture>%d%% -> res<%.0f Ohm-m; 95%% CI (stand-level anchor): %.0f..%.0f; 95%% PI (single tree): %.0f..%.0f\n",
              M, x, x - ci, x + ci, x - pi, x + pi))
}
