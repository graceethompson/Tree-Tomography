# SI table: stand context for the sampled trees (R1 #11).
# For each focal species at each site: share of live basal area, DBH quartiles
# of the stand, and the sampled DBH range.
# BGS = June 2025 variable-radius (prism) survey of the stand, 40 plots
#   (data/supplementary/BGS_VRP_2025.csv). In a prism tally every tree
#   represents the same basal area, so basal-area share = tally share; for the
#   stem-size quartiles each tree is weighted by its expansion factor
#   (trees per ha = BAF / tree basal area, proportional to 1/DBH^2).
# EMS = 2019 Harvard Forest ForestGEO census, live stems (Orwig, Foster &
#   Ellison 2023, HF253 v6, doi:10.6073/pasta/818789a882a318c1d7f3fc43a2289e12;
#   download hf253-06-stems-2019.csv to data/external/, not tracked in git).
# Quartiles use live trees >= 10 cm DBH at both sites.
suppressMessages(library(dplyr))
ours <- read.csv("data/Tree_ID_info.csv")
ours$sp <- c(rm = "acerru", bg = "nysssy", ro = "querru", hem = "tsugca")[ours$species]
ours$site <- ifelse(ours$plot == "BGS", "BGS", "EMS")
latin <- c(acerru = "A. rubrum", nysssy = "N. sylvatica", querru = "Q. rubra", tsugca = "T. canadensis")
vrp <- read.csv("data/supplementary/BGS_VRP_2025.csv") %>% filter(status == "alive", !is.na(dbh))
cen <- read.csv("data/external/hf253-06-stems-2019.csv") %>% filter(status == "A", !is.na(dbh))
ba_share <- list(BGS = with(vrp, tapply(dbh, species, length) / nrow(vrp)),
                 EMS = with(cen %>% mutate(ba = (dbh / 2)^2), tapply(ba, sp, sum) / sum(ba)))
q_bgs <- vrp %>% filter(dbh >= 10)
q_ems <- cen %>% group_by(tree.id) %>% slice_max(dbh, n = 1, with_ties = FALSE) %>% ungroup() %>% filter(dbh >= 10)
wq <- function(x, w, p) { o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w); sapply(p, function(pp) x[which(cw >= pp)[1]]) }
rows <- list()
for (site in c("BGS", "EMS")) for (s in c("acerru", "nysssy", "querru", "tsugca")) {
  o <- ours$dbh[ours$sp == s & ours$site == site]; if (!length(o)) next
  if (site == "BGS") { x <- q_bgs$dbh[q_bgs$species == s]; w <- 1 / x^2 } else { x <- q_ems$dbh[q_ems$sp == s]; w <- rep(1, length(x)) }
  qq <- wq(x, w, c(.25, .5, .75))
  rows[[length(rows) + 1]] <- data.frame(site = site, species = latin[[s]],
    pct_live_basal_area = round(100 * ba_share[[site]][[s]]),
    stand_dbh_p25_cm = round(qq[1]), stand_dbh_p50_cm = round(qq[2]), stand_dbh_p75_cm = round(qq[3]),
    n_sampled = length(o), sampled_dbh_min_cm = round(min(o) + 1e-9, 1), sampled_dbh_max_cm = round(max(o) + 1e-9, 1))
}
tab <- bind_rows(rows) %>% arrange(site, -pct_live_basal_area)
write.csv(tab, "analysis/revision/output/CJFR-stand-context.csv", row.names = FALSE)
print(tab, row.names = FALSE)
