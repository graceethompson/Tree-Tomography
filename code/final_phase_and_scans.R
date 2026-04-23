library(tidyverse)
library(grid)
library(jpeg)
library(patchwork)
library(readxl)

# ============================================================================
# 1. DATA PREP — same as phase_image_panels.R
# ============================================================================

setwd("/Users/jongewirtzman/My Drive/Research/Tomography/Tree-Tomography")

tree_info <- read_csv("data/Tree_ID_info.csv", show_col_types = FALSE)
ert_data  <- read_csv("data/ERT_application_results.csv", show_col_types = FALSE)

tree_info <- tree_info %>%
  mutate(site = case_when(
    plot %in% c("C1", "C2", "D1", "D2", "E1", "E5", "F1", "G1", "H1") ~ "EMS",
    TRUE ~ "BGS"
  ))

tree_info <- tree_info %>% mutate(tree = as.character(tree))
ert_data  <- ert_data %>% mutate(tree = as.character(tree))

# TRAINING SET: original 57 trees
dat_train <- tree_info %>%
  inner_join(ert_data, by = "tree") %>%
  mutate(dataset = "training")

# VALIDATION SET: 12 hemlock trees
hem_ert <- read_csv("data/hemlock/validation_summary.csv", show_col_types = FALSE)
hem_sot <- read_csv("data/hemlock/SOT_results.csv", show_col_types = FALSE)

hem_sot_dbh <- hem_sot %>%
  filter(str_detect(Filename, "_DBH\\.jpg")) %>%
  mutate(tree_id = str_replace(Filename, "_DBH\\.jpg", "")) %>%
  group_by(tree_id) %>%
  summarise(percent_damaged = mean(pct_damaged, na.rm = TRUE), .groups = "drop")

hem_val <- hem_ert %>%
  rename(tree_id = tree_id,
         median = Median, sd = SD, cv = CV,
         gini = Gini, entropy = Entropy, cma = CMA,
         radialgradient = RadialGradient) %>%
  mutate(mean = if ("Mean" %in% names(.)) Mean else 1000 / Conductance) %>%
  left_join(hem_sot_dbh, by = "tree_id") %>%
  mutate(
    tree = tree_id,
    species = "hem",
    site = "HF",
    plot = "HF",
    percent_solid_wood = NA_real_,
    decay = if_else(percent_damaged > 0, "present", "absent"),
    filename = paste0(tree_id, "_DBH.jpg"),
    npixels = NA_real_,
    dataset = "validation"
  )

dat <- bind_rows(dat_train, hem_val) %>%
  mutate(structural_loss = percent_damaged,
         abs_cma = abs(cma),
         abs_radgrad = abs(radialgradient),
         neg_mean   = -mean,
         neg_median = -median)

sot_threshold <- 1

# ============================================================================
# 2. COMPUTE PC1 (species-normalized, training-fitted PCA)
# ============================================================================

pca_metrics <- c("mean", "median", "sd", "cv", "gini", "entropy",
                 "cma", "radialgradient")

# Species means/sds from TRAINING data
train_spp_stats <- dat %>%
  filter(dataset == "training") %>%
  group_by(species) %>%
  summarise(across(all_of(pca_metrics),
                   list(mu = ~ mean(., na.rm = TRUE),
                        sigma = ~ sd(., na.rm = TRUE))),
            .groups = "drop")

# Species-normalize ALL trees using training stats
pca_spp_normed <- dat %>%
  select(tree, species, dataset, all_of(pca_metrics)) %>%
  left_join(train_spp_stats, by = "species")

for (m in pca_metrics) {
  mu_col <- paste0(m, "_mu")
  sd_col <- paste0(m, "_sigma")
  pca_spp_normed[[m]] <- (pca_spp_normed[[m]] - pca_spp_normed[[mu_col]]) /
    pca_spp_normed[[sd_col]]
}

pca_input_all <- pca_spp_normed %>% select(all_of(pca_metrics)) %>% as.matrix()
pca_input_all[is.nan(pca_input_all)] <- 0

# Fit PCA on training rows only
train_rows <- which(dat$dataset == "training")
pca_input_train <- pca_input_all[train_rows, ]
pca_fit <- prcomp(pca_input_train, center = FALSE, scale. = FALSE)

cat("PCA variance explained:\n")
print(summary(pca_fit)$importance[, 1:4])
cat("\nPC1 loadings:\n")
print(round(pca_fit$rotation[, 1], 3))

# Project ALL trees onto training PCA
all_scores <- pca_input_all %*% pca_fit$rotation
dat$pc1 <- all_scores[, 1]

# Flip so high = wet/anomalous (mean loading should be negative after flip)
if (pca_fit$rotation["mean", 1] > 0) {
  dat$pc1 <- -dat$pc1
  cat("Flipped PC1 so high = wet/anomalous\n")
}

ert_threshold <- mean(dat$pc1[dat$dataset == "training"], na.rm = TRUE)

cat("\nPC1 threshold (training mean):", round(ert_threshold, 3), "\n")

# ============================================================================
# 2b. PCA BIPLOT — scores + loading arrows
# ============================================================================

library(scales)

spp_labels <- c(hem = "Hemlock", rm = "Red maple",
                bg = "Black gum", ro = "Red oak")

# Scores (all trees, PC1 vs PC2)
pc1_scores <- all_scores[, 1]
pc2_scores <- all_scores[, 2]

# Apply same flip to PC2 display if PC1 was flipped
pc1_flip <- if (pca_fit$rotation["mean", 1] > 0) -1 else 1
pc1_scores <- pc1_scores * pc1_flip

biplot_df <- dat %>%
  mutate(PC1 = pc1_scores,
         PC2 = pc2_scores,
         species_label = spp_labels[species])

# Loadings — scale arrows to fit within score range
loadings <- as.data.frame(pca_fit$rotation[, 1:2])
loadings$metric <- rownames(loadings)
names(loadings)[1:2] <- c("PC1", "PC2")
loadings$PC1 <- loadings$PC1 * pc1_flip  # match flip

# Scale arrows: find a multiplier so arrows are visible relative to scores
arrow_scale <- min(max(abs(pc1_scores)), max(abs(pc2_scores))) * 0.8 /
               max(sqrt(loadings$PC1^2 + loadings$PC2^2))
loadings$PC1 <- loadings$PC1 * arrow_scale
loadings$PC2 <- loadings$PC2 * arrow_scale

# Pretty metric labels
metric_labels <- c(mean = "Mean Resistivity", median = "Median Resistivity", sd = "SD", cv = "CV",
                   gini = "Gini", entropy = "Entropy",
                   cma = "CMA", radialgradient = "RadGrad")
loadings$label <- metric_labels[loadings$metric]

# Variance explained
ve <- summary(pca_fit)$importance[2, 1:2] * 100

spp_shapes <- c(Hemlock = 16, "Red maple" = 17, "Black gum" = 15, "Red oak" = 18)

p_biplot <- ggplot(biplot_df, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = species_label), size = 3, alpha = 0.7, color = "grey30") +
  # Loading arrows
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "#B22222", linewidth = 0.7, inherit.aes = FALSE) +
  geom_text(data = loadings,
            aes(x = PC1 * 1.12, y = PC2 * 1.12, label = label),
            color = "#B22222", size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  scale_shape_manual(name = "Species", values = spp_shapes) +
  labs(x = paste0("PC1 (", round(ve[1], 1), "% variance)\nhigh = wet / anomalous"),
       y = paste0("PC2 (", round(ve[2], 1), "% variance)")) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave("output/figures/pca_biplot.pdf", p_biplot, width = 9, height = 7)
ggsave("output/figures/pca_biplot.png", p_biplot, width = 9, height = 7, dpi = 300, bg = "white")
cat("Saved: output/figures/pca_biplot.pdf and .png\n")

# Assign quadrants
dat <- dat %>%
  mutate(quadrant = case_when(
    structural_loss <= sot_threshold & pc1 <= ert_threshold ~ "I: Sound",
    structural_loss <= sot_threshold & pc1 >  ert_threshold ~ "II: Incipient",
    structural_loss >  sot_threshold & pc1 >  ert_threshold ~ "III: Active",
    structural_loss >  sot_threshold & pc1 <= ert_threshold ~ "IV: Cavity"
  ))

cat("\nQuadrant counts:\n")
print(table(dat$quadrant, dat$dataset))

# ============================================================================
# 3. FINAL PHASE DIAGRAM — all species, symmetrical, quadrant labels
# ============================================================================

# Signed sqrt transform for x-axis (PC1 spans negative to positive)
signed_sqrt_trans <- trans_new(
  name = "signed_sqrt",
  transform = function(x) sign(x) * sqrt(abs(x)),
  inverse = function(x) sign(x) * x^2
)

# Make x-axis symmetrical around threshold
pc1_range <- range(dat$pc1, na.rm = TRUE)
pc1_max_dev <- max(abs(pc1_range - ert_threshold))
pc1_lim <- ert_threshold + c(-1, 1) * pc1_max_dev * 1.15

# Y-axis: use sqrt transform to spread out low values
# Limits in original scale (sqrt applied by ggplot)
sot_max <- max(dat$structural_loss, na.rm = TRUE) * 1.1
sot_lim <- c(-1.5, sot_max)

dat <- dat %>%
  mutate(species_label = spp_labels[species],
         is_validation = dataset == "validation")

# Muted, CB-safe severity ramp: steel blue → gold → sienna → brick
quad_fill <- c("I: Sound" = "#4E79A7", "II: Incipient" = "#E5C460",
               "III: Active" = "#D4873F", "IV: Cavity" = "#C4524E")
# Darker versions for text labels
quad_text <- c("I: Sound" = "#355570", "II: Incipient" = "#A89030",
               "III: Active" = "#9A5F28", "IV: Cavity" = "#8C3535")

p_final <- ggplot(dat, aes(x = pc1, y = structural_loss)) +
  # Quadrant shading
  annotate("rect",
           xmin = pc1_lim[1], xmax = ert_threshold,
           ymin = -2, ymax = sot_threshold,
           fill = quad_fill[1], alpha = 0.07) +
  annotate("rect",
           xmin = ert_threshold, xmax = pc1_lim[2],
           ymin = -2, ymax = sot_threshold,
           fill = quad_fill[2], alpha = 0.07) +
  annotate("rect",
           xmin = ert_threshold, xmax = pc1_lim[2],
           ymin = sot_threshold, ymax = sot_max,
           fill = quad_fill[3], alpha = 0.07) +
  annotate("rect",
           xmin = pc1_lim[1], xmax = ert_threshold,
           ymin = sot_threshold, ymax = sot_max,
           fill = quad_fill[4], alpha = 0.07) +
  # Threshold lines
  geom_hline(yintercept = sot_threshold, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = ert_threshold, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  # Quadrant labels at axis extremes
  annotate("text", x = pc1_lim[1], y = sot_lim[1],
           label = "I: Sound", hjust = 0, vjust = 0,
           color = quad_text[1], fontface = "bold", size = 4.5) +
  annotate("text", x = pc1_lim[2], y = sot_lim[1],
           label = "II: Incipient", hjust = 1, vjust = 0,
           color = quad_text[2], fontface = "bold", size = 4.5) +
  annotate("text", x = pc1_lim[2], y = sot_max,
           label = "III: Active", hjust = 1, vjust = 1,
           color = quad_text[3], fontface = "bold", size = 4.5) +
  annotate("text", x = pc1_lim[1], y = sot_max,
           label = "IV: Cavity", hjust = 0, vjust = 1,
           color = quad_text[4], fontface = "bold", size = 4.5) +
  # All points: shape by species, same size, dark grey, slight vertical jitter
  geom_point(aes(shape = species_label),
             size = 3.5, alpha = 0.75, color = "grey20",
             position = position_jitter(width = 0, height = 0.1, seed = 42)) +
  scale_shape_manual(name = "Species", values = spp_shapes) +
  scale_x_continuous(trans = signed_sqrt_trans,
                     breaks = c(-4, -2, -1, 0, 1, 2, 4)) +
  scale_y_continuous(trans = signed_sqrt_trans,
                     breaks = c(0, 1, 5, 10, 20, 30, 40, 50)) +
  coord_cartesian(xlim = pc1_lim, ylim = sot_lim) +
  labs(x = "ERT PC1 (species-normalized, sqrt scale)\nhigh = wet / anomalous",
       y = "Structural Loss (%)\n(SoT percent damaged, sqrt scale)") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("output/figures/final_phase_diagram.pdf", p_final, width = 9, height = 8)
ggsave("output/figures/final_phase_diagram.png", p_final, width = 9, height = 8, dpi = 300, bg = "white")
cat("\nSaved: output/figures/final_phase_diagram.pdf and .png\n")

# ============================================================================
# 4. COMPUTE PC1 FOR HEMLOCK PER-HEIGHT DATA (for scans_by_height figures)
# ============================================================================

# Load per-height hemlock ERT data
hem_ert_raw <- read_csv("data/hemlock/ERT_results_2026-03-05.csv", show_col_types = FALSE)

hem_ert_height <- hem_ert_raw %>%
  mutate(
    base = str_remove(Filename, "\\.jpg$"),
    height = str_extract(base, "[^_]+$"),
    tree_id = str_remove(base, "_[^_]+$")
  ) %>%
  select(-base) %>%
  group_by(tree_id, height) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  rename(mean = Mean, median = Median, sd = SD, cv = CV,
         gini = Gini, entropy = Entropy, cma = CMA,
         radialgradient = RadialGradient) %>%
  mutate(abs_cma = abs(cma),
         abs_radgrad = abs(radialgradient))

# Species-normalize using TRAINING hemlock stats (same as main analysis)
hem_stats <- train_spp_stats %>% filter(species == "hem")

for (m in pca_metrics) {
  mu_val <- hem_stats[[paste0(m, "_mu")]]
  sd_val <- hem_stats[[paste0(m, "_sigma")]]
  hem_ert_height[[paste0(m, "_norm")]] <- (hem_ert_height[[m]] - mu_val) / sd_val
}

# Project onto training PCA
hem_pca_input <- hem_ert_height %>%
  select(all_of(paste0(pca_metrics, "_norm"))) %>%
  as.matrix()
colnames(hem_pca_input) <- pca_metrics
hem_pca_input[is.nan(hem_pca_input)] <- 0

hem_scores <- hem_pca_input %*% pca_fit$rotation
hem_ert_height$pc1 <- hem_scores[, 1]

# Apply same flip
if (pca_fit$rotation["mean", 1] > 0) {
  hem_ert_height$pc1 <- -hem_ert_height$pc1
}

# Load core moisture
mc <- read_excel("data/hemlock/MC_Tomo_paper_Jon.xlsx", sheet = 1)
names(mc) <- c("tree_id", "moisture")
hem_ert_height <- hem_ert_height %>% left_join(mc, by = "tree_id")

# Load SoT per height
hem_sot_all <- hem_sot %>%
  mutate(
    base = str_remove(Filename, "\\.jpg$"),
    height = str_extract(base, "[^_]+$"),
    tree_id = str_remove(base, "_[^_]+$")
  ) %>%
  select(tree_id, height, pct_damaged) %>%
  group_by(tree_id, height) %>%
  summarise(pct_damaged = mean(pct_damaged, na.rm = TRUE), .groups = "drop")

hem_ert_height <- hem_ert_height %>%
  left_join(hem_sot_all, by = c("tree_id", "height"))

# Trees with any height data (include those missing a height — blank space for missing)
multi_height_trees <- hem_ert_height %>%
  group_by(tree_id) %>%
  filter(n_distinct(height) >= 2) %>%
  ungroup()

cat("\nTrees with 2+ heights:", n_distinct(multi_height_trees$tree_id), "\n")

# Sort by core moisture
tree_order <- multi_height_trees %>%
  select(tree_id, moisture) %>%
  distinct() %>%
  arrange(moisture)

cat("Tree order by moisture:", paste(tree_order$tree_id, collapse = ", "), "\n")

# ============================================================================
# 5. HELPER: Create scans-by-height figure
# ============================================================================

hem_ert_img_dir <- "/Users/jongewirtzman/My Drive/Research/Tomography/Hemlock_Tomography/ERT JPEGS_Absolute"

make_scans_by_height <- function(tree_data, tree_order, metric_rows,
                                  output_file, title_text = "") {
  # tree_data: data frame with tree_id, height, moisture, metrics, pc1
  # tree_order: ordered tree_id vector
  # metric_rows: list of list(label, values, fmt) — one per summary row
  # output_file: PDF path

  n_trees <- length(tree_order)
  heights <- c("Upper", "DBH", "Lower")
  n_heights <- length(heights)
  n_metric_rows <- length(metric_rows)

  # Layout dimensions
  img_row_h <- 0.22  # fraction of page per image row
  metric_row_h <- 0.04  # fraction per metric row
  top_margin <- 0.06
  bottom_margin <- 0.02
  left_margin <- 0.10  # for row labels
  right_margin <- 0.02

  total_img_h <- n_heights * img_row_h
  total_metric_h <- n_metric_rows * metric_row_h
  total_h <- total_img_h + total_metric_h

  # Scale if needed
  avail_h <- 1 - top_margin - bottom_margin
  scale <- min(1, avail_h / total_h)
  img_row_h <- img_row_h * scale
  metric_row_h <- metric_row_h * scale

  col_w <- (1 - left_margin - right_margin) / n_trees

  pdf(output_file, width = 3.5 * n_trees, height = 14)
  grid.newpage()

  # Title
  pushViewport(viewport(y = 1 - top_margin / 2, height = top_margin, just = "center"))
  grid.text(title_text, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()

  # Tree ID headers
  for (j in seq_along(tree_order)) {
    tid <- tree_order[j]
    x_center <- left_margin + (j - 0.5) * col_w
    pushViewport(viewport(x = x_center, y = 1 - top_margin,
                          width = col_w, height = 0.03, just = "top"))
    grid.text(tid, gp = gpar(fontsize = 11, fontface = "bold"))
    popViewport()
  }

  # Image rows
  for (i in seq_along(heights)) {
    h <- heights[i]
    y_top <- 1 - top_margin - 0.03 - (i - 1) * img_row_h

    # Row label
    pushViewport(viewport(x = left_margin / 2,
                          y = y_top - img_row_h / 2,
                          width = left_margin, height = img_row_h))
    grid.text(h, gp = gpar(fontsize = 12, fontface = "bold"))
    popViewport()

    # Horizontal separator line
    if (i > 1) {
      grid.lines(x = c(left_margin, 1 - right_margin),
                 y = c(y_top, y_top),
                 gp = gpar(col = "blue", lwd = 1.5))
    }

    for (j in seq_along(tree_order)) {
      tid <- tree_order[j]
      x_center <- left_margin + (j - 0.5) * col_w

      img_file <- file.path(hem_ert_img_dir, paste0(tid, "_", h, ".jpg"))

      pushViewport(viewport(x = x_center, y = y_top - img_row_h / 2,
                            width = col_w * 0.92, height = img_row_h * 0.92))

      if (file.exists(img_file)) {
        tryCatch({
          img <- readJPEG(img_file)
          grid.raster(img, interpolate = TRUE)
        }, error = function(e) NULL)
      }

      popViewport()

      # Vertical separator
      if (j > 1) {
        x_sep <- left_margin + (j - 1) * col_w
        grid.lines(x = c(x_sep, x_sep),
                   y = c(y_top, y_top - img_row_h),
                   gp = gpar(col = "grey80", lwd = 0.5))
      }
    }
  }

  # Metric summary rows below images
  metric_y_start <- 1 - top_margin - 0.03 - n_heights * img_row_h

  for (k in seq_along(metric_rows)) {
    mr <- metric_rows[[k]]
    y_center <- metric_y_start - (k - 0.5) * metric_row_h

    # Horizontal line above
    grid.lines(x = c(left_margin * 0.3, 1 - right_margin),
               y = c(y_center + metric_row_h / 2, y_center + metric_row_h / 2),
               gp = gpar(col = "grey60", lwd = 0.8))

    # Row label
    pushViewport(viewport(x = left_margin / 2, y = y_center,
                          width = left_margin, height = metric_row_h))
    grid.text(mr$label, gp = gpar(fontsize = 9, fontface = "bold"))
    popViewport()

    # Values with color coding
    vals <- mr$values
    fmt <- mr$fmt
    # Color scale: blue-white or custom
    val_range <- range(vals, na.rm = TRUE)

    for (j in seq_along(tree_order)) {
      x_center_j <- left_margin + (j - 0.5) * col_w
      v <- vals[j]

      # Color: interpolate from light to dark blue based on value
      if (diff(val_range) > 0) {
        frac <- (v - val_range[1]) / diff(val_range)
      } else {
        frac <- 0.5
      }
      # Light grey to dark blue
      bg_col <- rgb(1 - 0.6 * frac, 1 - 0.6 * frac, 1 - 0.2 * frac)

      pushViewport(viewport(x = x_center_j, y = y_center,
                            width = col_w, height = metric_row_h))
      grid.rect(gp = gpar(fill = bg_col, col = "grey80", lwd = 0.5))
      grid.text(sprintf(fmt, v), gp = gpar(fontsize = 10,
                                            col = if (frac > 0.7) "white" else "black"))
      popViewport()
    }
  }

  dev.off()
  cat("Saved:", output_file, "\n")
}

# ============================================================================
# 6. VERSION A: fig_scans_by_height with PC1 ADDED to existing metrics
# ============================================================================

# Get DBH-level values for each tree (for summary rows)
dbh_data <- multi_height_trees %>% filter(height == "DBH")

# Align to tree_order
dbh_ordered <- tibble(tree_id = tree_order$tree_id) %>%
  left_join(dbh_data, by = "tree_id")

metric_rows_add <- list(
  list(label = "Core %\nMoisture",
       values = tree_order$moisture,
       fmt = "%.0f%%"),
  list(label = "ERT Mean\nConductance\n(mS)",
       values = 1000 / dbh_ordered$mean,
       fmt = "%.1f"),
  list(label = "ERT CV",
       values = dbh_ordered$cv,
       fmt = "%.2f"),
  list(label = "ERT PC1",
       values = dbh_ordered$pc1,
       fmt = "%.2f")
)

make_scans_by_height(multi_height_trees, tree_order$tree_id,
                      metric_rows_add,
                      "output/hemlock_figures/fig_scans_by_height_with_pc1.pdf",
                      "ERT Scans by Height (sorted by core moisture)")

# ============================================================================
# 7. VERSION B: fig_scans_by_height with PC1 REPLACING other ERT metrics
# ============================================================================

metric_rows_replace <- list(
  list(label = "Core %\nMoisture",
       values = tree_order$moisture,
       fmt = "%.0f%%"),
  list(label = "ERT PC1",
       values = dbh_ordered$pc1,
       fmt = "%.2f")
)

make_scans_by_height(multi_height_trees, tree_order$tree_id,
                      metric_rows_replace,
                      "output/hemlock_figures/fig_scans_by_height_pc1_only.pdf",
                      "ERT Scans by Height (sorted by core moisture)")

# ============================================================================
# 8. FIG BEST PREDICTOR: PC1 vs Core Moisture (hemlock validation)
# ============================================================================

library(ggrepel)

# Use DBH-level PC1 for each hemlock tree (all 12)
hem_dbh <- hem_ert_height %>%
  filter(height == "DBH") %>%
  select(tree_id, moisture, pc1, mean) %>%
  mutate(conductance = 1000 / mean)  # mS (inverse of resistivity)

fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

val_theme <- theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.tag = element_text(face = "bold"))

# --- Panel A (left): Mean conductance vs moisture ---
ct_cond   <- cor.test(hem_dbh$moisture, hem_dbh$conductance, method = "pearson")
ct_cond_s <- cor.test(hem_dbh$moisture, hem_dbh$conductance, method = "spearman", exact = FALSE)
lm_cond   <- lm(conductance ~ moisture, data = hem_dbh)
lm_cond_sum <- summary(lm_cond)

ann_cond <- paste0(
  "r = ", round(ct_cond$estimate, 2), " (", fmt_p(ct_cond$p.value), ")\n",
  "\u03c1 = ", round(ct_cond_s$estimate, 2), " (", fmt_p(ct_cond_s$p.value), ")\n",
  "R\u00b2 = ", round(lm_cond_sum$r.squared, 2)
)

pA <- ggplot(hem_dbh, aes(x = moisture, y = conductance)) +
  geom_smooth(method = "lm", se = TRUE, color = "#B2182B",
              fill = "#B2182B", alpha = 0.15, linewidth = 0.8) +
  geom_point(size = 3, color = "#B2182B") +
  geom_text_repel(aes(label = tree_id), size = 2.8, max.overlaps = 20) +
  annotate("text",
           x = min(hem_dbh$moisture) + diff(range(hem_dbh$moisture)) * 0.02,
           y = max(hem_dbh$conductance) - diff(range(hem_dbh$conductance)) * 0.02,
           label = ann_cond,
           hjust = 0, vjust = 1, size = 3.2, color = "grey25", lineheight = 1.2) +
  labs(x = "Core Moisture (%)",
       y = "Mean Conductance (mS)",
       tag = "A") +
  val_theme

# --- Panel B (right): PC1 vs moisture ---
ct_p  <- cor.test(hem_dbh$moisture, hem_dbh$pc1, method = "pearson")
ct_s  <- cor.test(hem_dbh$moisture, hem_dbh$pc1, method = "spearman", exact = FALSE)
lm_pc1 <- lm(pc1 ~ moisture, data = hem_dbh)
lm_pc1_sum <- summary(lm_pc1)

ann_pc1 <- paste0(
  "r = ", round(ct_p$estimate, 2), " (", fmt_p(ct_p$p.value), ")\n",
  "\u03c1 = ", round(ct_s$estimate, 2), " (", fmt_p(ct_s$p.value), ")\n",
  "R\u00b2 = ", round(lm_pc1_sum$r.squared, 2)
)

pB <- ggplot(hem_dbh, aes(x = moisture, y = pc1)) +
  geom_smooth(method = "lm", se = TRUE, color = "#2166AC",
              fill = "#2166AC", alpha = 0.15, linewidth = 0.8) +
  geom_point(size = 3, color = "#2166AC") +
  geom_text_repel(aes(label = tree_id), size = 2.8, max.overlaps = 20) +
  annotate("text",
           x = min(hem_dbh$moisture) + diff(range(hem_dbh$moisture)) * 0.02,
           y = max(hem_dbh$pc1) - diff(range(hem_dbh$pc1)) * 0.02,
           label = ann_pc1,
           hjust = 0, vjust = 1, size = 3.2, color = "grey25", lineheight = 1.2) +
  labs(x = "Core Moisture (%)",
       y = "ERT PC1 (species-normalized)\nhigh = wet / anomalous",
       tag = "B") +
  val_theme

fig_bp <- pA + pB

ggsave("output/hemlock_figures/fig_best_predictor_pc1.pdf", fig_bp,
       width = 12, height = 5.5, device = cairo_pdf)
ggsave("output/hemlock_figures/fig_best_predictor_pc1.png", fig_bp,
       width = 12, height = 5.5, dpi = 300, bg = "white")
cat("Saved: output/hemlock_figures/fig_best_predictor_pc1.pdf and .png\n")

# ============================================================================
# 9. QUADRANT DISTRIBUTION BY SPECIES AND SITE (training set only)
# ============================================================================

dat_train_q <- dat %>%
  filter(dataset == "training") %>%
  mutate(species_label = spp_labels[species],
         quadrant = factor(quadrant,
                           levels = c("I: Sound", "II: Incipient",
                                      "III: Active", "IV: Cavity")))

quad_colors <- quad_fill  # same CB-safe palette as phase diagram

# --- By Species ---
spp_counts <- dat_train_q %>%
  count(species_label, quadrant, .drop = FALSE) %>%
  group_by(species_label) %>%
  mutate(pct = n / sum(n) * 100,
         total = sum(n)) %>%
  ungroup()

# --- By Site ---
site_counts <- dat_train_q %>%
  count(site, quadrant, .drop = FALSE) %>%
  group_by(site) %>%
  mutate(pct = n / sum(n) * 100,
         total = sum(n)) %>%
  ungroup()

# --- Percent stacked charts ---
p_spp_pct <- ggplot(spp_counts, aes(x = species_label, y = pct, fill = quadrant)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  geom_text(aes(label = ifelse(pct > 5, paste0(round(pct), "%"), ""),
                color = quadrant),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(name = "Decay Phase", values = quad_colors) +
  scale_color_manual(values = c("I: Sound" = "white", "II: Incipient" = "grey20",
                                "III: Active" = "white", "IV: Cavity" = "white")) +
  labs(x = NULL, y = "Percent of Trees") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
  )

p_site_pct <- ggplot(site_counts, aes(x = site, y = pct, fill = quadrant)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.5) +
  geom_text(aes(label = ifelse(pct > 5, paste0(round(pct), "%"), ""),
                color = quadrant),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(name = "Decay Phase", values = quad_colors) +
  scale_color_manual(values = c("I: Sound" = "white", "II: Incipient" = "grey20",
                                "III: Active" = "white", "IV: Cavity" = "white")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right"
  )

# Combined: phase diagram on top, distribution plots on bottom
fig_combined <- p_final / (p_spp_pct | p_site_pct) +
  plot_layout(heights = c(2, 1))

ggsave("output/figures/quadrant_distribution.pdf", fig_combined, width = 12, height = 14)
ggsave("output/figures/quadrant_distribution.png", fig_combined, width = 12, height = 14, dpi = 300, bg = "white")
cat("Saved: output/figures/quadrant_distribution.pdf\n")

cat("\n=== ALL DONE ===\n")
