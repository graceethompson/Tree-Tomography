# Per-tree SI table (merges the former sensor table and scheme table):
# species, site, DBH, sensor count/paths/spacing, recorded % damage, mean
# resistivity, species-normalized PC1, and the decay category under the
# manuscript's PC1 scheme, the species-median scheme and the absolute scheme.
suppressMessages(library(dplyr))
sens <- read.csv("analysis/revision/output/CJFR-sensor-counts.csv")
sch  <- read.csv("analysis/revision/output/scheme_assignments.csv")
tab <- sch %>% select(tree, species = sp, site, dbh, percent_damaged, mean_resistivity = mean, pc1,
                      category_pc1 = cat_pc1, category_species_median = cat_species_median, category_absolute = cat_absolute) %>%
  inner_join(sens %>% select(tree, sensor_count, paths, spacing), by = "tree") %>%
  mutate(mean_resistivity = round(mean_resistivity), pc1 = round(pc1, 2), spacing = round(spacing, 1)) %>%
  select(tree, species, site, dbh, sensor_count, paths, spacing, percent_damaged, mean_resistivity, pc1,
         category_pc1, category_species_median, category_absolute) %>% arrange(site, species, tree)
stopifnot(nrow(tab) == 57)
write.csv(tab, "analysis/revision/output/CJFR-tree-table.csv", row.names = FALSE)
cat("wrote CJFR-tree-table.csv:", nrow(tab), "rows x", ncol(tab), "cols\n"); print(head(tab, 3))
