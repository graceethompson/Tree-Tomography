# R port of capstone.py — site-direction robustness across axes and thresholds.
# Run from repo root: Rscript analysis/revision/scripts/capstone.R
source("analysis/revision/scripts/revision_common.R")

sd_pop <- function(x) sqrt(mean((x - mean(x))^2))  # numpy ndarray.std(), ddof=0

m <- load_merged()
bp <- build_pc1(m, "species")
pc <- bp$pc1
m$pc1 <- pc
m$struct <- m$percent_damaged > 1
mb <- m$site == "BGS" & !m$struct
me <- m$site == "EMS" & !m$struct
nb <- sum(mb); ne <- sum(me)
# species-standardized wetness (higher=wetter)
xdev <- xdev_species_median(m)
# pooled PC1
pcp <- build_pc1(m, "pooled")$pc1

sdv <- sd_pop(pc); sdx <- sd_pop(xdev); sdp <- sd_pop(pcp)
ks <- c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1)
# each axis: (values where higher=more anomalous/wet, thresholds spanning -1..+1 SD)
axes <- list(
  "PC1 (within-sp)" = list(pc, mean(pc) + ks * sdv),
  "PC1 (pooled)"    = list(pcp, mean(pcp) + ks * sdp),
  "species-median resistivity" = list(xdev, 0 + ks * sdx),
  "absolute resistivity" = list(-m$mean,
    sort(-quantile(m$mean, c(0.35, 0.42, 0.5, 0.58, 0.65), type = 7, names = FALSE))),
  "CMA" = list(m$cma, c(0.25, 0.30, 0.33, 0.36, 0.40))
)

cat("SITE DIRECTION ROBUSTNESS: incipient among structurally-sound trees, BGS vs EMS\n")
cat(sprintf("%-28s%5s%15s%16s%16s\n", "axis", "#thr", "BGS>EMS holds", "BGS% range", "EMS% range"))
total <- 0; held <- 0
for (name in names(axes)) {
  vals <- axes[[name]][[1]]; ths <- axes[[name]][[2]]
  bs <- numeric(0); es <- numeric(0); ok <- 0
  for (t0 in ths) {
    anom <- vals > t0
    b <- sum(anom & mb) / nb * 100
    e <- sum(anom & me) / ne * 100
    bs <- c(bs, b); es <- c(es, e); ok <- ok + (b >= e)
  }
  total <- total + length(ths); held <- held + ok
  cat(sprintf("%-28s%5d%8d/%-6d%6.0f-%-9.0f%6.0f-%-9.0f\n",
              name, length(ths), ok, length(ths), min(bs), max(bs), min(es), max(es)))
}
cat(sprintf("\nACROSS ALL %d axis×threshold combinations: BGS>=EMS incipient in %d/%d = %.0f%%\n",
            total, held, total, held / total * 100))

# structural prevalence stability by SoT threshold
cat("\nStructural-loss prevalence by site, across SoT thresholds:\n")
for (tau in c(1, 2, 3, 5)) {
  s <- m$percent_damaged > tau
  b <- mean(s[m$site == "BGS"]) * 100
  e <- mean(s[m$site == "EMS"]) * 100
  cat(sprintf("  SoT>%d%%: BGS %.0f%% vs EMS %.0f%%  (EMS>BGS: %s)\n",
              tau, b, e, if (e > b) "True" else "False"))
}

# species severity ordering stability (mean % damage among decayed, by species)
cat("\nSpecies severity (mean % damage among decayed trees) — ordering:\n")
dd <- m[m$percent_damaged > 0, ]
sv <- sort(round(tapply(dd$percent_damaged, dd$sp, mean), 1), decreasing = TRUE)
cat("   {", paste(sprintf("'%s': %.1f", names(sv), sv), collapse = ", "), "}\n", sep = "")
cnt <- table(dd$sp)
cat("  decayed n per species: {", paste(sprintf("'%s': %d", names(cnt), cnt), collapse = ", "), "}\n", sep = "")
