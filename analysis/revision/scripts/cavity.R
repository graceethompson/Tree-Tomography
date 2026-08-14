# Cavity (IV) vs active (III) among the structurally-damaged trees: composite
# PC1 rule vs a physical crack/dry-core rule. R port of cavity.py. Run from the
# repo root:
#   Rscript analysis/revision/scripts/cavity.R

source("analysis/revision/scripts/revision_common.R")

# UTF-8 locale so the png devices render em-dash / sqrt / arrow glyphs
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"))) invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))

m <- load_merged()
pca <- build_pc1(m, "species")
pc <- pca$pc1
m$pc1 <- pc
mu <- mean(pc)
dm <- m[m$percent_damaged > 1, ]
dm$base <- ifelse(dm$pc1 > mu, "III active", "IV cavity")

# Physically-motivated cavity signals among damaged trees:
# (a) crack detected  (b) dry core: low CMA (few low-resistivity pixels in inner third)
cat("Damaged trees — is there a physical 'dry cavity' signature?\n")
cat(sprintf("CMA (central moisture accumulation; LOW = dry core): III active mean=%.2f, IV cavity mean=%.2f\n",
            mean(dm$cma[dm$base == "III active"]), mean(dm$cma[dm$base == "IV cavity"])))
# alternative dry-core cavity rule
dm$dry_core <- dm$cma < 0.15
dm$crack <- dm$crack_detected == "yes"
cat("\nProposed physical cavity rule = structural loss AND (crack OR dry core CMA<0.15):\n")
dm$cav_phys <- dm$crack | dm$dry_core
dm <- dm[order(dm$pc1), ]
pybool <- function(x) ifelse(x, "True", "False")
cat(sprintf("%5s %12s %16s %6s %6s %9s %9s %10s %9s\n",
            "tree", "sp", "percent_damaged", "crack", "cma", "dry_core", "pc1", "base", "cav_phys"))
for (i in seq_len(nrow(dm))) {
  r <- dm[i, ]
  cat(sprintf("%5d %12s %16d %6s %6.4f %9s %9.6f %10s %9s\n",
              r$tree, r$sp, r$percent_damaged, pybool(r$crack), r$cma,
              pybool(r$dry_core), r$pc1, r$base, pybool(r$cav_phys)))
}
cat(sprintf("\nComposite (PC1<mean) cavities: %d  | physical (crack/dry-core) cavities: %d\n",
            sum(dm$base == "IV cavity"), sum(dm$cav_phys)))
cat("Agreement between the two cavity definitions:",
    sum((dm$base == "IV cavity") == dm$cav_phys), "/", nrow(dm), "\n")

# gap check on PC1 among damaged trees
s <- sort(dm$pc1)
cat(sprintf("\nGap at the III/IV boundary: cavities at PC1<= %.2f; next active at %.2f; threshold mean=%.2f sits in a %.2f-wide gap\n",
            s[3], s[4], mu, s[4] - s[3]))

# ---------- figure ----------
png(file.path(OUT_DIR, "CJFR-cavity-fig.png"),
    width = 9, height = 5.2, units = "in", res = 170)
par(mar = c(4.2, 4.6, 3.2, 0.8), mgp = c(2.6, 0.7, 0))
plot(NA, xlim = c(-1.4, 5.3), ylim = c(-0.01, 0.48), xlab = "", ylab = "", axes = FALSE)
rect(par("usr")[1], 0, par("usr")[2], 0.15,
     col = adjustcolor("#8e6b3f", alpha.f = 0.10), border = NA)
abline(v = mu, col = "#c0392b", lty = 2, lwd = 1.3)
text(mu + 0.05, 0.02, "III/IV threshold\n(PC1 mean)", cex = 0.6, col = "#c0392b", adj = 0)
text(5.28, 0.075, "dry core\n(cavity-like)", cex = 0.6, col = "#6b4f2f", adj = 1)
for (i in seq_len(nrow(dm))) {
  r <- dm[i, ]
  sz <- 90 + r$percent_damaged * 10          # matplotlib area (pt^2)
  points(r$pc1, r$cma, pch = 21,
         bg = adjustcolor(SPECIES_COLS[[r$sp]], alpha.f = 0.9),
         col = if (r$crack) "black" else "white",
         lwd = if (r$crack) 2 else 0.6, cex = sqrt(sz) / 6.5)
  text(r$pc1 + 0.09, r$cma + 0.012, sprintf("%d", r$tree), cex = 0.55, adj = c(0, 0))
}
axis(1, cex.axis = 0.8)
axis(2, cex.axis = 0.8, las = 1)
box()
title(xlab = "ERT PC1  (low = drier/normal → cavity side; high = wet → active)", cex.lab = 0.85)
title(ylab = "CMA — central moisture accumulation\n(low = dry core)", cex.lab = 0.85)
title(main = "The 15 structurally-damaged trees: cavity vs active\n(marker size = % structural loss; black outline = crack detected)",
      cex.main = 0.85)
legend("topleft",
       legend = c(names(SPECIES_COLS), "crack detected"),
       pch = 21,
       pt.bg = c(unname(SPECIES_COLS), "grey"),
       col = c(rep("white", length(SPECIES_COLS)), "black"),
       pt.lwd = c(rep(1, length(SPECIES_COLS)), 2),
       pt.cex = 1.1, cex = 0.62, bty = "n")
invisible(dev.off())
cat("saved figure\n")
