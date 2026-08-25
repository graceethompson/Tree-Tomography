# Vertical variation in decay columns from the multi-height validation scans
# (Reviewer 1, comment #10).
#
# The 12 validation hemlocks were imaged below (Lower), at (DBH), and above
# (Upper) breast height (images/hemlock_SoT, images/hemlock_ERT). The reviewer
# asks that these be used to estimate how much a single BH measurement plane
# under- or over-states decay, since butt rots are widest near the ground.
#
# Data: data/hemlock/SOT_results.csv (per-image SoT % damaged) and
# data/hemlock/ERT_results_2026-03-05.csv (per-image ERT metrics; replicate
# analyses of the same image are averaged).
#
# Run from the repo root: Rscript analysis/revision/scripts/vertical_variation.R

source("analysis/revision/scripts/revision_common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- SoT: % damaged by height ----
sot <- read.csv("data/hemlock/SOT_results.csv", fileEncoding = "UTF-8-BOM")
names(sot) <- c("filename", "pct_damaged")
sot$tree <- sub("^(HF_[0-9]+)_.*$", "\\1", sot$filename)
sot$height <- sub("^HF_[0-9]+_(.*)\\.jpg$", "\\1", sot$filename)

# HF_031205 originally listed "Lower" twice (2 then 0) with no Upper row;
# the source file was corrected 2026-08-24 after inspecting the tomograms
# (the Upper scan shows the small central defect = 2%; the Lower scan is
# clean = 0%). Guard against any re-introduced duplicates:
stopifnot(!anyDuplicated(sot[c("tree", "height")]))
sot$height <- factor(sot$height, levels = c("Lower", "DBH", "Upper"))

wide <- reshape(sot[c("tree", "height", "pct_damaged")],
                idvar = "tree", timevar = "height", direction = "wide")
names(wide) <- sub("pct_damaged\\.", "", names(wide))
multi <- wide[complete.cases(wide[c("Lower", "DBH", "Upper")]), ]

cat("=== SoT % damaged by height (validation hemlocks) ===\n")
print(wide, row.names = FALSE)
cat("\nTrees with all three heights: n =", nrow(multi), "\n")

cat("\nPer-tree vertical profile (Lower -> DBH -> Upper):\n")
multi$range <- apply(multi[c("Lower", "DBH", "Upper")], 1, function(x) max(x) - min(x))
multi$max_h <- c("Lower", "DBH", "Upper")[apply(multi[c("Lower", "DBH", "Upper")], 1, which.max)]
print(multi, row.names = FALSE)

cat("\nWithin-tree spread: mean range =", round(mean(multi$range), 1),
    "percentage points; max =", max(multi$range), "\n")
cat("Height with the largest damage: ", paste(multi$max_h, collapse = ", "), "\n")

# Does the DBH plane understate the maximum decay in the stem section?
multi$max3 <- pmax(multi$Lower, multi$DBH, multi$Upper)
cat("\nDBH vs max across heights:\n")
for (i in seq_len(nrow(multi)))
  cat(sprintf("  %s: DBH %d%% vs max %d%% (at %s)\n",
              multi$tree[i], multi$DBH[i], multi$max3[i], multi$max_h[i]))
cat("Mean (max - DBH) =", round(mean(multi$max3 - multi$DBH), 2), "points\n")

# paired tests (small n — descriptive)
wl <- wilcox.test(multi$Lower, multi$DBH, paired = TRUE, exact = FALSE)
wu <- wilcox.test(multi$Upper, multi$DBH, paired = TRUE, exact = FALSE)
cat(sprintf("\nWilcoxon signed-rank (n=%d): Lower vs DBH p=%.2f; Upper vs DBH p=%.2f\n",
            nrow(multi), wl$p.value, wu$p.value))

# structural-loss classification agreement across planes (manuscript threshold >1%)
for (tau in c(1, 5)) {
  cls <- data.frame(tree = multi$tree,
                    DBH = multi$DBH > tau,
                    any_height = multi$max3 > tau)
  n_miss <- sum(cls$any_height & !cls$DBH)
  cat(sprintf("Threshold >%d%%: DBH-only misses %d of %d trees flagged at some height\n",
              tau, n_miss, sum(cls$any_height)))
}

# ---- ERT: metrics by height (replicates averaged) ----
ert <- read.csv("data/hemlock/ERT_results_2026-03-05.csv")
ert$tree <- sub("^(HF_[0-9]+)_.*$", "\\1", ert$Filename)
ert$height <- sub("^HF_[0-9]+_(.*)\\.jpg$", "\\1", ert$Filename)
erta <- aggregate(cbind(Mean, Median, SD, CV, Gini, Entropy, CMA, RadialGradient) ~ tree + height,
                  ert, mean)
erta$height <- factor(erta$height, levels = c("Lower", "DBH", "Upper"))
multi_ert <- names(which(table(erta$tree) == 3))

cat("\n=== ERT mean resistivity (Ohm-m) by height, trees with 3 planes ===\n")
ew <- reshape(erta[erta$tree %in% multi_ert, c("tree", "height", "Mean")],
              idvar = "tree", timevar = "height", direction = "wide")
names(ew) <- sub("Mean\\.", "", names(ew))
ew <- ew[c("tree", "Lower", "DBH", "Upper")]
ew[-1] <- round(ew[-1])
print(ew, row.names = FALSE)

# within-tree vs between-tree variability in mean resistivity
sub <- erta[erta$tree %in% multi_ert, ]
within_cv <- aggregate(Mean ~ tree, sub, function(x) sd(x) / mean(x))
tree_means <- aggregate(Mean ~ tree, sub, mean)
cat(sprintf("\nMean within-tree CV of mean resistivity across heights: %.2f\n",
            mean(within_cv$Mean)))
cat(sprintf("Between-tree CV of tree-mean resistivity: %.2f\n",
            sd(tree_means$Mean) / mean(tree_means$Mean)))

# ---- figure ----
p1 <- ggplot(sot[sot$tree %in% multi$tree, ],
             aes(height, pct_damaged, group = tree, colour = tree)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  annotate("text", x = "Lower", y = 1, label = "1% threshold", vjust = -0.6,
           hjust = 1.1, size = 2.8, colour = "grey30") +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 2.4) +
  labs(x = NULL, y = "SoT % of section damaged",
       title = "A",
       colour = NULL) +
  theme_classic(base_size = 10)

p2 <- ggplot(sub, aes(height, Mean, group = tree, colour = tree)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 2.4) +
  labs(x = NULL, y = expression("ERT mean resistivity" ~ (Omega %.% m)),
       title = "B",
       colour = NULL) +
  theme_classic(base_size = 10)

fig <- p1 | p2
ggsave(file.path(OUT_DIR, "CJFR-vertical-variation.png"), fig,
       width = 11, height = 4.6, dpi = 170)
cat("\nsaved", file.path(OUT_DIR, "CJFR-vertical-variation.png"), "\n")
