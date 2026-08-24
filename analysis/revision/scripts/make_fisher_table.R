# Writes output/CJFR-fisher-tests.csv (Table S-FISHER): Fisher exact tests of
# the wetland-upland incipient contrast at each candidate ERT threshold, the
# whole-sample version, and the generalists-only subset. Same computations as
# statistical_tests.R, formatted as a table.
# Run from repo root: Rscript analysis/revision/scripts/make_fisher_table.R
source("analysis/revision/scripts/revision_common.R")

m <- load_merged()
pc <- build_pc1(m, "species")$pc1
m$pc1 <- pc
m$struct <- m$percent_damaged > 1
sdv <- sd(pc)
sound <- !m$struct

row_for <- function(label, denom, a, n1, c, n2) {
  b <- n1 - a; d <- n2 - c
  or <- (a * d) / (b * c)
  p <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value
  data.frame(contrast = label, denominator = denom,
             bgs_incipient = a, bgs_n = n1, bgs_pct = round(a / n1 * 100),
             ems_incipient = c, ems_n = n2, ems_pct = round(c / n2 * 100),
             odds_ratio = ifelse(is.finite(or), round(or, 2), Inf),
             p_value = round(p, 3))
}

rows <- list()
for (cs in list(list(-0.5, "mean - 0.5 SD"), list(0, "mean (published)"),
                list(0.5, "mean + 0.5 SD"), list(1, "mean + 1 SD"))) {
  t0 <- mean(pc) + cs[[1]] * sdv
  inc <- sound & (pc > t0)
  rows[[length(rows) + 1]] <- row_for(
    paste0("PC1 threshold ", cs[[2]]), "structurally sound trees",
    sum(inc & m$site == "BGS"), sum(sound & m$site == "BGS"),
    sum(inc & m$site == "EMS"), sum(sound & m$site == "EMS"))
}

inc0 <- sound & (pc > mean(pc))
rows[[length(rows) + 1]] <- row_for(
  "PC1 threshold mean (published)", "all trees",
  sum(inc0 & m$site == "BGS"), sum(m$site == "BGS"),
  sum(inc0 & m$site == "EMS"), sum(m$site == "EMS"))

g <- m[m$sp %in% c("A.rubrum", "T.canadensis"), ]
incg <- (!g$struct) & (g$pc1 > mean(pc))
rows[[length(rows) + 1]] <- row_for(
  "PC1 threshold mean (published)", "generalist species only (all trees)",
  sum(incg & g$site == "BGS"), sum(g$site == "BGS"),
  sum(incg & g$site == "EMS"), sum(g$site == "EMS"))

out <- do.call(rbind, rows)
write.csv(out, "analysis/revision/output/CJFR-fisher-tests.csv", row.names = FALSE)
print(out)
cat("wrote analysis/revision/output/CJFR-fisher-tests.csv\n")
