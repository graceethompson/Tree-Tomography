# Per-tree category assignments under each classification scheme (R port of
# make_assignments_table.py). Writes analysis/revision/output/scheme_assignments.csv.
# Run from repo root: Rscript analysis/revision/scripts/make_assignments_table.R
source("analysis/revision/scripts/revision_common.R")
options(scipen = 12)  # keep small values (e.g. cma 0.0005) in decimal notation

m <- load_merged()
bp <- build_pc1(m, "species")
pc <- bp$pc1
m$pc1 <- round(pc, 3)
m$struct <- m$percent_damaged > 1
m$xdev <- round(xdev_species_median(m), 3)

m$cat_pc1 <- four_cat(pc > mean(pc), m$struct)
m$cat_species_median <- four_cat(m$xdev > 0, m$struct)
m$cat_absolute <- four_cat(m$mean < median(m$mean), m$struct)
m$cat_6cell <- ifelse(!m$struct & m$xdev > 0.5, "II",
              ifelse(!m$struct & m$xdev >= -0.5, "I-II",
              ifelse(!m$struct, "I",
              ifelse(m$struct & m$xdev < -0.5, "IV",
              ifelse(m$struct & m$xdev <= 0.5, "III-IV", "III")))))

out <- m[, c("tree", "sp", "site", "dbh", "percent_damaged", "mean", "cma", "pc1", "xdev",
             "cat_pc1", "cat_species_median", "cat_absolute", "cat_6cell")]
write.csv(out, file.path(OUT_DIR, "scheme_assignments.csv"), row.names = FALSE, quote = FALSE)
cat(sprintf("wrote %d rows; 3-way identical (pc1/median/absolute): %d/57\n", nrow(out),
            sum(out$cat_pc1 == out$cat_species_median & out$cat_species_median == out$cat_absolute)))
