# How were the green (intermediate-velocity) zones handled, and does it matter?
# (Reviewer 1, reproducibility comment #2.)
#
# PiCUS tomograms color each pixel by sonic velocity class: browns = solid wood,
# green = intermediate, violet = decayed, blue = cavity. The manuscript describes
# percent damage as "non-brown colors combined", but the recorded numbers say
# otherwise: in data/Tree_ID_info.csv, percent_solid_wood matches the "Solid
# wood: X %" printed by PiCUS on every tomogram (verified against the image
# headers, analysis/revision/data/sot_image_annotations.csv), and
# percent_solid_wood + percent_damaged < 100 for most damaged trees. The gap is
# the green zone: the recorded percent_damaged EXCLUDED green (violet+blue only).
#
# This script (1) quantifies the green share per tree from PiCUS's own
# tabulation, (2) independently re-measures all color-class areas from the
# tomogram pixels, and (3) re-runs the structural classification and the site
# contrast under the three green treatments (excluded / counted as damaged /
# counted as solid).
#
# Run from the repo root: Rscript analysis/revision/scripts/green_zone.R

source("analysis/revision/scripts/revision_common.R")
suppressPackageStartupMessages({
  library(jpeg)
  library(ggplot2)
  library(patchwork)
})

df <- load_main()
df$green_pct <- 100 - df$percent_solid_wood - df$percent_damaged

cat("=== PiCUS's own tabulation (Tree_ID_info.csv) ===\n")
cat("percent_solid_wood = PiCUS 'Solid wood %' (brown only; header-verified 56/57)\n")
cat("percent_damaged    = violet + blue (green EXCLUDED)\n")
cat("green share = 100 - solid - damaged\n\n")
cat("Green share distribution (57 trees):\n")
print(summary(df$green_pct))
cat("Trees with green > 5%:", sum(df$green_pct > 5),
    "; > 10%:", sum(df$green_pct > 10), "\n\n")

# ---- independent pixel re-measurement from the tomogram images ----
# images/main_SoT also holds three empty exports (366, 428, 677: ellipse only)
# from trees that are not part of the 57-tree study set; tree 412's first file
# is likewise empty, so its 26vii24 rescan is used.
scan_for <- function(tree) {
  if (tree == 412) return("images/main_SoT/412_26vii24.jpg")
  find_scan("images/main_SoT", tree)
}

classify_pixels <- function(path) {
  img <- readJPEG(path)
  h <- dim(img)[1]
  body <- img[35:h, , , drop = FALSE]
  rgbm <- rbind(as.vector(body[, , 1]), as.vector(body[, , 2]), as.vector(body[, , 3]))
  hsv <- rgb2hsv(rgbm, maxColorValue = 1)
  H <- hsv[1, ] * 360; S <- hsv[2, ]; V <- hsv[3, ]
  cls <- rep(NA_character_, length(H))
  keep <- V >= 0.15 & S >= 0.12
  ellipse <- keep & S > 0.85 & H >= 220 & H <= 255        # pure-blue outline
  rednum <- keep & (H < 10 | H > 345) & S > 0.5 & V > 0.4 # red sensor numbers
  ok <- keep & !ellipse & !rednum
  cls[ok & H >= 10 & H < 75] <- "brown"
  cls[ok & H >= 75 & H < 180] <- "green"
  cls[ok & H >= 180 & H <= 235 & S <= 0.85] <- "blue"
  cls[ok & H > 235 & H <= 345] <- "violet"
  tab <- table(factor(cls, levels = c("brown", "green", "violet", "blue")))
  100 * tab / sum(tab)
}

cat("=== Independent pixel re-measurement (hue classification) ===\n")
pix <- do.call(rbind, lapply(df$tree, function(t) {
  p <- classify_pixels(scan_for(t))
  data.frame(tree = t, pix_brown = p["brown"], pix_green = p["green"],
             pix_violet = p["violet"], pix_blue = p["blue"], row.names = NULL)
}))
m <- merge(df, pix, by = "tree", all.x = TRUE)
m$pix_damage_ex_green <- m$pix_violet + m$pix_blue
m$pix_damage_in_green <- m$pix_green + m$pix_violet + m$pix_blue

ok <- !is.na(m$pix_brown)
cat(sprintf("Trees measured: %d of 57\n", sum(ok)))
cat(sprintf("pixel violet+blue vs recorded percent_damaged: r = %.3f, MAE = %.1f pts\n",
            cor(m$pix_damage_ex_green[ok], m$percent_damaged[ok]),
            mean(abs(m$pix_damage_ex_green[ok] - m$percent_damaged[ok]))))
cat(sprintf("pixel non-brown vs 100 - percent_solid_wood:   r = %.3f, MAE = %.1f pts\n",
            cor(m$pix_damage_in_green[ok], 100 - m$percent_solid_wood[ok]),
            mean(abs(m$pix_damage_in_green[ok] - (100 - m$percent_solid_wood[ok])))))
cat(sprintf("pixel green vs PiCUS green share:              r = %.3f, MAE = %.1f pts\n\n",
            cor(m$pix_green[ok], m$green_pct[ok]),
            mean(abs(m$pix_green[ok] - m$green_pct[ok]))))

# ---- sensitivity: the three green treatments ----
treatments <- list(
  "green EXCLUDED (as recorded)" = df$percent_damaged,
  "green as DAMAGED (non-brown)" = 100 - df$percent_solid_wood,
  "green as SOLID"               = df$percent_damaged   # identical by construction
)
cat("=== Structural classification & site contrast under green treatments ===\n")
cat("(green-as-solid equals the recorded values: damaged classes never included green)\n\n")
for (nm in names(treatments)[1:2]) {
  dmg <- treatments[[nm]]
  for (tau in c(1, 5)) {
    s <- dmg > tau
    b <- mean(s[df$site == "BGS"]) * 100; e <- mean(s[df$site == "EMS"]) * 100
    cat(sprintf("%-30s tau>%d%%: n=%2d flagged | BGS %2.0f%% vs EMS %2.0f%%\n",
                nm, tau, sum(s), b, e))
  }
  a1 <- anova(lm(dmg ~ site, data = df))
  cat(sprintf("%-30s one-way site ANOVA: F=%.2f, p=%.3f\n\n",
              nm, a1$`F value`[1], a1$`Pr(>F)`[1]))
}

flip <- (treatments[[2]] > 1) != (df$percent_damaged > 1)
cat("Trees whose structural-loss call (>1%) flips if green counts as damage:",
    sum(flip), "of 57\n")
print(df[flip, c("tree", "sp", "site", "percent_solid_wood", "percent_damaged")],
      row.names = FALSE)

# ---- figure ----
pal <- c(brown = "#8a5a2b", green = "#2eaf5d", violet = "#a6188f", blue = "#3a7ebf")
long <- m[ok, c("tree", "sp", "percent_damaged", "pix_green", "pix_violet", "pix_blue")]
long <- long[order(-long$percent_damaged, -long$pix_green), ]
long$tree_f <- factor(long$tree, levels = long$tree)
stack <- reshape(long, direction = "long",
                 varying = c("pix_green", "pix_violet", "pix_blue"),
                 v.names = "pct", timevar = "class",
                 times = c("green", "violet", "blue"))
stack$class <- factor(stack$class, levels = c("blue", "violet", "green"))

p1 <- ggplot(m[ok, ], aes(percent_damaged, pix_damage_ex_green)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(aes(colour = "violet+blue (green excluded)"), size = 2, alpha = 0.8) +
  geom_point(aes(y = pix_damage_in_green, colour = "non-brown (green included)"),
             size = 2, alpha = 0.8) +
  scale_colour_manual(values = c("violet+blue (green excluded)" = "#a6188f",
                                 "non-brown (green included)" = "#2eaf5d")) +
  labs(x = "recorded percent_damaged", y = "pixel-measured damaged area (%)",
       colour = NULL,
       title = "A") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")

p2 <- ggplot(stack, aes(tree_f, pct, fill = class)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = pal[c("blue", "violet", "green")]) +
  labs(x = "tree (sorted by recorded damage)", y = "% of section",
       fill = "velocity class",
       title = "B") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5.5),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "CJFR-green-zone.png"), p1 | p2,
       width = 12.5, height = 5.2, dpi = 170)
cat("\nsaved", file.path(OUT_DIR, "CJFR-green-zone.png"), "\n")
