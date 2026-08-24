# Sensor counts, stem size, and SoT effective resolution (Reviewer 1, comment #8).
#
# The reviewer asks that sensor counts and stem diameters be reported for all
# trees, and that we demonstrate whether stem size and sensor count covary with
# site, species, or the percent-damage estimate (SoT reconstruction error is
# resolution-dependent, and resolution scales with sensor count and stem size).
#
# Sensor counts were read from the numbered sensor positions on every PiCUS
# tomogram (analysis/revision/data/sot_image_annotations.csv; the same file's
# 'Solid wood %' header values verify against data/Tree_ID_info.csv in 56/57
# trees). Sensor spacing (circumference / sensor count) is used as the
# effective-resolution proxy.
#
# Run from the repo root: Rscript analysis/revision/scripts/sensor_counts.R

source("analysis/revision/scripts/revision_common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

ann <- read.csv("analysis/revision/data/sot_image_annotations.csv")
ann$tree <- suppressWarnings(as.integer(sub("_.*", "", ann$filename)))
main_ann <- ann[!is.na(ann$tree), ]
# one sensor count per tree (duplicate scans agree except where noted in the csv)
sc <- aggregate(sensor_count ~ tree, main_ann, function(x) max(x))

df <- merge(load_main(), sc, by = "tree")
df$struct <- df$percent_damaged > 1
df$spacing <- pi * df$dbh / df$sensor_count   # cm between sensors
df$paths <- df$sensor_count * (df$sensor_count - 1) / 2  # sensor-pair travel paths

cat("=== Sensor counts, all 57 main trees ===\n")
print(table(sensors = df$sensor_count))
cat("Sensor-pair paths per tomogram: 7 sensors -> 21 paths; 8 -> 28\n")
print(table(paths = df$paths))
cat("\nBy site:\n"); print(table(df$site, df$sensor_count))
cat("\nBy species:\n"); print(table(df$sp, df$sensor_count))

cat("\nFisher exact, sensor count (7 vs 8) x site: p =",
    signif(fisher.test(table(df$site, df$sensor_count))$p.value, 3), "\n")

cat("\n=== Stem size ===\n")
cat("DBH by site (cm):\n")
print(aggregate(dbh ~ site, df, function(x) round(c(mean = mean(x), sd = sd(x),
                                                    min = min(x), max = max(x)), 1)))
tt <- t.test(dbh ~ site, data = df)
cat(sprintf("Welch t-test DBH BGS vs EMS: t=%.2f, p=%.4f (EMS larger)\n",
            tt$statistic, tt$p.value))

cat("\n=== Effective resolution: sensor spacing (cm of circumference per sensor) ===\n")
print(aggregate(spacing ~ site, df, function(x) round(c(mean = mean(x), sd = sd(x),
                                                        min = min(x), max = max(x)), 1)))
print(aggregate(spacing ~ sp, df, function(x) round(c(mean = mean(x), sd = sd(x)), 1)))
st <- t.test(spacing ~ site, data = df)
cat(sprintf("Welch t-test spacing BGS vs EMS: t=%.2f, p=%.4f\n", st$statistic, st$p.value))

cat("\n=== Does resolution covary with the damage estimate? ===\n")
r1 <- cor.test(df$spacing, df$percent_damaged, method = "spearman", exact = FALSE)
cat(sprintf("spacing vs percent_damaged: Spearman rho=%.3f, p=%.3f\n",
            r1$estimate, r1$p.value))
r2 <- cor.test(df$dbh, df$percent_damaged, method = "spearman", exact = FALSE)
cat(sprintf("dbh vs percent_damaged:     Spearman rho=%.3f, p=%.3f\n",
            r2$estimate, r2$p.value))
g <- glm(struct ~ spacing, family = binomial, df)
cat(sprintf("structural loss (>1%%) ~ spacing logistic: beta=%.3f, p=%.3f\n",
            coef(g)[2], summary(g)$coefficients[2, 4]))

# does the site percent-damage contrast survive controlling for spacing?
a <- anova(lm(percent_damaged ~ spacing + site, df))
cat("\nANOVA percent_damaged ~ spacing + site (Type I, spacing first):\n")
print(round(as.data.frame(a), 4))

# smallest resolvable defect: a defect spanning less than about one sensor
# spacing is poorly constrained; express spacing as % of cross-section area
# for a circular defect of diameter = spacing/2 (illustrative)
df$defect_frac <- 100 * (df$spacing / 2)^2 / df$dbh^2   # (d/2)^2 / D^2 = area ratio
cat(sprintf("\nIllustrative smallest-resolvable-defect area (circle of diameter = spacing/2):\n"))
print(round(summary(df$defect_frac), 2))
cat("-> compare with the 1% structural threshold\n")

# hemlock validation set
hem <- ann[is.na(ann$tree), ]
cat("\n=== Validation hemlocks (per scan) ===\n")
print(table(sensors = hem$sensor_count))

# ---- independent cross-check of the visually-read sensor counts ----
# The sensor markers are small olive-green rings (dark green, B < 45) that are
# distinguishable from the legend green and the tomogram's decay-green class
# (both have B ~ 86). Cluster those pixels and count clusters; decay-zone
# boundary pixels can create false extra clusters on heavily damaged trees, so
# a pixel count ABOVE the visual read on a damaged tree is expected noise
# (verified by eye for 351 and 583), while a count BELOW would flag a missed
# sensor.
suppressPackageStartupMessages(library(jpeg))
count_dots <- function(path) {
  img <- readJPEG(path)
  R <- img[, , 1] * 255; G <- img[, , 2] * 255; B <- img[, , 3] * 255
  mask <- G >= 85 & G <= 175 & R < 115 & B < 45 & (G - R) > 25
  mask[1:32, ] <- FALSE                       # header bar
  idx <- which(mask, arr.ind = TRUE)
  if (nrow(idx) < 5) return(0)
  cl <- cutree(hclust(dist(idx), method = "single"), h = 6)
  sum(table(cl) >= 8)
}
main_files <- main_ann$filename
dots <- vapply(main_files,
               function(f) count_dots(file.path("images/main_SoT", f)), numeric(1))
cmp <- data.frame(filename = main_files,
                  visual = main_ann$sensor_count, pixel = dots)
cat("\n=== Pixel dot-count cross-check (main survey scans) ===\n")
cat("agree:", sum(cmp$visual == cmp$pixel), "of", nrow(cmp),
    "| pixel < visual (possible missed sensor):", sum(cmp$pixel < cmp$visual), "\n")
if (any(cmp$visual != cmp$pixel))
  print(cmp[cmp$visual != cmp$pixel, ], row.names = FALSE)

# ---- per-tree table for the supplement ----
out <- df[order(df$site, df$sp, df$tree),
          c("tree", "sp", "site", "dbh", "sensor_count", "paths", "spacing", "percent_damaged")]
out$spacing <- round(out$spacing, 1)
write.csv(out, file.path(OUT_DIR, "CJFR-sensor-counts.csv"), row.names = FALSE)

# ---- figure ----
p1 <- ggplot(df, aes(dbh, spacing, colour = site, shape = factor(sensor_count))) +
  geom_point(size = 2.6, alpha = 0.85) +
  scale_colour_manual(values = c(BGS = "#1f77b4", EMS = "#d9822b")) +
  labs(x = "DBH (cm)", y = "sensor spacing (cm)",
       shape = "sensors", colour = "site",
       title = "A") +
  theme_classic(base_size = 10)

set.seed(3)
p2 <- ggplot(df, aes(spacing, percent_damaged, colour = site)) +
  geom_jitter(height = 0.3, width = 0, size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = c(BGS = "#1f77b4", EMS = "#d9822b")) +
  labs(x = "sensor spacing (cm)", y = "SoT % damaged") +
  theme_classic(base_size = 10)

# Panel A (spacing vs DBH) dropped: it repeats Table S-SENSORS content.
# The figure is now the single covariance panel (damage vs spacing).
ggsave(file.path(OUT_DIR, "CJFR-sensor-counts.png"), p2,
       width = 5.8, height = 4.4, dpi = 170)
cat("\nsaved CJFR-sensor-counts.png and CJFR-sensor-counts.csv\n")
