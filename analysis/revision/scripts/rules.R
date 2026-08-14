# R port of rules.py — graded moisture bands and III/IV physical rule.
# Run from repo root: Rscript analysis/revision/scripts/rules.R
source("analysis/revision/scripts/revision_common.R")
options(width = 160)

m <- load_merged()
bp <- build_pc1(m, "species")
pc <- bp$pc1
m$pc1 <- pc
m$res_z <- bp$Z$mean
mu <- mean(pc); sdv <- sd(pc)   # ddof=1
otsu <- 1.00  # from prior analysis

m$struct <- m$percent_damaged > 1
# graded moisture bands
band <- function(p) {
  if (p <= mu) "normal" else if (p <= otsu) "possible" else "confident"
}
m$mband <- vapply(m$pc1, band, character(1))

# pandas value_counts: descending count, ties in order of first appearance
value_counts <- function(x) {
  ord <- unique(x)
  cnt <- table(factor(x, levels = ord))
  cnt[order(-cnt)]
}
dict_int <- function(cnt) {
  paste0("{", paste(sprintf("'%s': %d", names(cnt), cnt), collapse = ", "), "}")
}

cat(sprintf("PC1: mean=%.2f, SD=%.2f; Otsu 'confident' break=1.0\n", mu, sdv))
cat("Moisture-anomaly bands (all 57):", dict_int(value_counts(m$mband)), "\n")
cat("Among SOUND trees (I/II split):", dict_int(value_counts(m$mband[!m$struct])), "\n")

# III/IV physical rule among damaged trees: cavity = dry core = high resistivity AND low central moisture
dmg <- m[m$struct, ]
dmg$cavity_phys <- dmg$mean >= 300 & dmg$cma < 0.15
dmg$class_phys <- ifelse(dmg$cavity_phys, "IV cavity", "III active")
dmg$class_composite <- ifelse(dmg$pc1 > mu, "III active", "IV cavity")
dmg <- dmg[order(dmg$mean), ]
cat("\nIII/IV among 15 damaged trees:\n")
tab <- data.frame(tree = dmg$tree, sp = dmg$sp, percent_damaged = dmg$percent_damaged,
                  mean = sprintf("%.2f", dmg$mean), cma = sprintf("%.4f", dmg$cma),
                  pc1 = sprintf("%.6f", dmg$pc1),
                  class_composite = dmg$class_composite, class_phys = dmg$class_phys)
print(tab, row.names = FALSE)
cat("composite cavities:", sum(dmg$class_composite == "IV cavity"),
    "| physical(dry-core) cavities:", sum(dmg$class_phys == "IV cavity"), "\n")

# absolute equivalents
cat("\nAbsolute translation (approx, species-confounded): PC1 mean ~ 369 Ohm-m ~ 103% moisture; Otsu(1.0) ~ 387 Ohm-m ~ 101% moisture\n")
gap <- sort(m$percent_damaged[m$percent_damaged > 0])
cat("Structural gap check: damage values >0: [", paste(gap, collapse = ", "), "]\n", sep = "")
