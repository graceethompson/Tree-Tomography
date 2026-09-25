# Table S5: sampled DBH relative to stand-level size distributions (R1 #11).
# BGS reference = June 2025 variable-radius (prism) survey of the stand
#   (data/supplementary/BGS_VRP_2025.csv); live trees, each weighted by
#   1/basal area so the distribution is per tree, not per unit basal area.
# EMS reference = 2019 Harvard Forest ForestGEO census, live trees (main stem)
#   >= 10 cm DBH. Orwig, Foster & Ellison 2023, HF253 v6,
#   doi:10.6073/pasta/818789a882a318c1d7f3fc43a2289e12 — download
#   hf253-06-stems-2019.csv to data/external/ (not tracked in git).
suppressMessages(library(dplyr))
ours <- read.csv("data/Tree_ID_info.csv")
ours$sp <- c(rm = "acerru", bg = "nysssy", ro = "querru", hem = "tsugca")[ours$species]
ours$site <- ifelse(ours$plot == "BGS", "BGS", "EMS")
latin <- c(acerru = "A. rubrum", nysssy = "N. sylvatica", querru = "Q. rubra", tsugca = "T. canadensis")

vrp <- read.csv("data/supplementary/BGS_VRP_2025.csv") %>% filter(status == "alive", !is.na(dbh), dbh >= 10)
cen <- read.csv("data/external/hf253-06-stems-2019.csv") %>% filter(status == "A", !is.na(dbh)) %>%
  group_by(tree.id) %>% slice_max(dbh, n = 1, with_ties = FALSE) %>% ungroup() %>% filter(dbh >= 10)

wq <- function(x, w, p) { o <- order(x); x <- x[o]; w <- cumsum(w[o]) / sum(w); sapply(p, function(pp) x[which(w >= pp)[1]]) }
wpct <- function(x, w, q) sapply(q, function(qq) 100 * sum(w[x <= qq]) / sum(w))

rows <- list()
for (site in c("BGS", "EMS")) for (s in c("acerru", "nysssy", "querru", "tsugca")) {
  o <- ours$dbh[ours$sp == s & ours$site == site]
  if (!length(o)) next
  if (site == "BGS") { x <- vrp$dbh[vrp$species == s]; w <- 1 / x^2; ref <- "BGS prism survey 2025 (live trees >= 10 cm, weighted by trees per ha)" }
  else { x <- cen$dbh[cen$sp == s]; w <- rep(1, length(x)); ref <- "ForestGEO census 2019 (live trees >= 10 cm)" }
  qq <- wq(x, w, c(.25, .5, .75)); pc <- wpct(x, w, o)
  rows[[length(rows) + 1]] <- data.frame(site = site, species = latin[[s]], reference = ref, reference_n = length(x),
    ref_p25_cm = round(qq[1]), ref_median_cm = round(qq[2]), ref_p75_cm = round(qq[3]),
    sampled_n = length(o), sampled_dbh_min_cm = min(o), sampled_dbh_max_cm = max(o),
    sampled_percentile_min = round(min(pc)), sampled_percentile_max = round(max(pc)))
}
tab <- bind_rows(rows)
write.csv(tab, "analysis/revision/output/CJFR-dbh-reference.csv", row.names = FALSE)
print(tab[, -3])
