library(tidyverse)
library(grid)
library(jpeg)
library(png)
library(patchwork)

# ============================================================================
# 1. DATA PREP (same as phase_diagram_analysis.R)
# ============================================================================

tree_info <- read_csv("data/Tree_ID_info.csv", show_col_types = FALSE)
ert_data  <- read_csv("data/ERT_application_results.csv", show_col_types = FALSE)

tree_info <- tree_info %>%
  mutate(site = case_when(
    plot %in% c("C1", "C2", "D1", "D2", "E1", "E5", "F1", "G1", "H1") ~ "EMS",
    TRUE ~ "BGS"
  ))

tree_info <- tree_info %>% mutate(tree = as.character(tree))
ert_data  <- ert_data %>% mutate(tree = as.character(tree))

# TRAINING SET: original 57 trees (define axes + thresholds from these)
dat_train <- tree_info %>%
  inner_join(ert_data, by = "tree") %>%
  mutate(dataset = "training")

# VALIDATION SET: 12 Hemlock Tomography trees (project onto training-defined axes)
hem_ert <- read_csv("data/hemlock/validation_summary.csv", show_col_types = FALSE)
hem_sot <- read_csv("data/hemlock/SOT_results.csv", show_col_types = FALSE)

hem_sot_dbh <- hem_sot %>%
  filter(str_detect(Filename, "_DBH\\.jpg")) %>%
  mutate(tree_id = str_replace(Filename, "_DBH\\.jpg", "")) %>%
  group_by(tree_id) %>%
  summarise(percent_damaged = mean(pct_damaged, na.rm = TRUE), .groups = "drop")

hem_val <- hem_ert %>%
  rename(tree_id = tree_id,
         mean = Mean, median = Median, sd = SD, cv = CV,
         gini = Gini, entropy = Entropy, cma = CMA,
         radialgradient = RadialGradient) %>%
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

# Combine for plotting, but normalizations computed from training set only
dat <- bind_rows(dat_train, hem_val) %>%
  mutate(structural_loss = percent_damaged,
         abs_cma = abs(cma),
         abs_radgrad = abs(radialgradient),
         neg_mean   = -mean,
         neg_median = -median)

cat("Training set:", sum(dat$dataset == "training"), "trees\n")
cat("Validation set (hemlock):", sum(dat$dataset == "validation"), "trees\n")
cat("Total:", nrow(dat), "trees\n")

sot_threshold <- 1

# Sound-wood baseline (from TRAINING set only — no detected damage)
sound_trees <- dat %>% filter(percent_damaged == 0, dataset == "training")

# ============================================================================
# 2. IMAGE PATHS
# ============================================================================

base_dir <- "images"
sot_dir  <- file.path(base_dir, "main_SoT")
ert_norm_dir <- file.path(base_dir, "main_ERT")
ert_abs_dir  <- file.path(base_dir, "main_ERT")

# Hemlock image directories (separate study)
hem_base <- "/Users/jongewirtzman/My Drive/Research/Tomography/Hemlock_Tomography"
hem_sot_dir <- file.path(hem_base, "SoTs_Percent_Calculations")
hem_ert_norm_dir <- file.path(hem_base, "ERT JPEGS_Absolute")

# Map tree IDs to filenames — use different paths for hemlock vs original
dat <- dat %>%
  mutate(
    sot_path = if_else(dataset == "validation",
                       file.path(hem_sot_dir, filename),
                       file.path(sot_dir, filename)),
    ert_norm_path = if_else(dataset == "validation",
                            file.path(hem_ert_norm_dir, filename),
                            file.path(ert_norm_dir, filename)),
    ert_abs_path = if_else(dataset == "validation",
                           file.path(hem_base, "ERT JPEGS_Absolute", filename),
                           file.path(ert_abs_dir, filename)),
    sot_exists = file.exists(sot_path),
    ert_exists = file.exists(ert_norm_path)
  )

cat("SoT images found:", sum(dat$sot_exists), "of", nrow(dat), "\n")
cat("ERT images found:", sum(dat$ert_exists), "of", nrow(dat), "\n")

# ============================================================================
# 3. DEFINE CANDIDATE AXES TO PLOT
# ============================================================================

# =========================================================================
# NOTE: All normalizations use TRAINING SET statistics only,
# then apply those same means/SDs to the validation hemlocks.
# Hemlock validation trees are normalized using the training hemlock stats
# (species = "hem" in training set) since site = "HF" has no training data.
# For species×site z-norms, hemlock validation uses species-only stats.
# =========================================================================

# Helper: compute group means/sds from training, apply to all
znorm_from_training <- function(dat, metric_col, group_vars) {
  train_stats <- dat %>%
    filter(dataset == "training") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(mu = mean(.data[[metric_col]], na.rm = TRUE),
              sigma = sd(.data[[metric_col]], na.rm = TRUE),
              .groups = "drop")

  dat %>%
    left_join(train_stats, by = group_vars) %>%
    mutate(z_value = if_else(!is.na(sigma) & sigma > 0,
                             (.data[[metric_col]] - mu) / sigma,
                             .data[[metric_col]] - mu)) %>%
    pull(z_value)
}

# --- Candidate A: neg_median, species×site z-norm (moisture LEVEL, high=wet) ---
# For validation hemlocks (site=HF), fall back to species-only normalization
dat$neg_median_spp_site_z <- znorm_from_training(dat, "neg_median", c("species", "site"))
# Fill NAs (hemlock validation with no training match for species×site) using species-only
na_mask <- is.na(dat$neg_median_spp_site_z)
if (any(na_mask)) {
  dat$neg_median_spp_site_z[na_mask] <-
    znorm_from_training(dat, "neg_median", "species")[na_mask]
}

# Threshold from training set only
ert_thresh_A <- mean(dat$neg_median_spp_site_z[dat$dataset == "training"], na.rm = TRUE)

# --- Candidate B: cv, species×site z-norm (moisture HETEROGENEITY) ---
dat$cv_spp_site_z <- znorm_from_training(dat, "cv", c("species", "site"))
na_mask <- is.na(dat$cv_spp_site_z)
if (any(na_mask)) {
  dat$cv_spp_site_z[na_mask] <- znorm_from_training(dat, "cv", "species")[na_mask]
}

# Find threshold using Otsu
find_otsu <- function(values) {
  values <- values[!is.na(values)]
  best <- -Inf; best_t <- median(values)
  for (t in quantile(values, seq(0.1, 0.9, 0.02))) {
    w0 <- mean(values <= t); w1 <- 1 - w0
    if (w0 > 0 & w1 > 0) {
      bv <- w0 * w1 * (mean(values[values <= t]) - mean(values[values > t]))^2
      if (bv > best) { best <- bv; best_t <- t }
    }
  }
  best_t
}

ert_thresh_B <- find_otsu(dat$cv_spp_site_z[dat$dataset == "training"])

# --- Candidate C: gini, species z-norm (heterogeneity) ---
dat$gini_spp_z <- znorm_from_training(dat, "gini", "species")
ert_thresh_C <- mean(dat$gini_spp_z[dat$dataset == "training"], na.rm = TRUE)

# --- Candidate D: cma, species×site z-norm (spatial structure) ---
dat$cma_spp_site_z <- znorm_from_training(dat, "cma", c("species", "site"))
na_mask <- is.na(dat$cma_spp_site_z)
if (any(na_mask)) {
  dat$cma_spp_site_z[na_mask] <- znorm_from_training(dat, "cma", "species")[na_mask]
}
ert_thresh_D <- mean(dat$cma_spp_site_z[dat$dataset == "training"], na.rm = TRUE)

# --- Candidate E: neg_mean, Lutz species baseline (high=wet) ---
# Baseline from training sound trees only
sound_baselines_spp <- sound_trees %>%
  group_by(species) %>%
  summarise(bl_neg_mean = mean(neg_mean, na.rm = TRUE), .groups = "drop")

dat <- dat %>%
  left_join(sound_baselines_spp, by = "species", suffix = c("", ".bl")) %>%
  mutate(neg_mean_lutz_spp = neg_mean - bl_neg_mean)

ert_thresh_E <- find_otsu(dat$neg_mean_lutz_spp[dat$dataset == "training"])

# --- Candidate F: PCA PC1 on species-normalized ERT metrics ---
# FIT PCA on training data only, then PROJECT validation hemlocks
pca_metrics <- c("mean", "median", "sd", "cv", "gini", "entropy",
                 "abs_cma", "abs_radgrad")

# Compute species means/sds from TRAINING data
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

cat("\nPCA on species-normalized ERT metrics (training set):\n")
print(summary(pca_fit)$importance[, 1:4])
cat("\nPC loadings:\n")
print(round(pca_fit$rotation[, 1:3], 3))

# Project ALL trees (training + validation) onto training PCA
all_scores <- pca_input_all %*% pca_fit$rotation

dat$pc1_spp <- all_scores[, 1]
dat$pc2_spp <- all_scores[, 2]
dat$pc3_spp <- all_scores[, 3]

# Flip so high = wet/anomalous
if (pca_fit$rotation["mean", 1] > 0) {
  dat$pc1_spp <- -dat$pc1_spp
  cat("Flipped PC1 so high = wet/anomalous\n")
}

ert_thresh_F <- mean(dat$pc1_spp[dat$dataset == "training"], na.rm = TRUE)
ert_thresh_G <- mean(dat$pc2_spp[dat$dataset == "training"], na.rm = TRUE)

# --- Candidate H: PCA PC1 on species×site normalized (training-defined) ---
# Use species-only normalization (same as above) since validation has no site match
# This makes H identical to F in this setup, so skip duplicating
ert_thresh_H <- ert_thresh_F
dat$pc1_sppsite <- dat$pc1_spp
dat$pc2_sppsite <- dat$pc2_spp

# Assign quadrants for each candidate
assign_quad <- function(sot, ert, sot_t, ert_t) {
  case_when(
    sot <= sot_t & ert <= ert_t ~ "I: Sound",
    sot <= sot_t & ert >  ert_t ~ "II: Incipient",
    sot >  sot_t & ert >  ert_t ~ "III: Active",
    sot >  sot_t & ert <= ert_t ~ "IV: Cavity"
  )
}

dat <- dat %>%
  mutate(
    quad_A = assign_quad(structural_loss, neg_median_spp_site_z, sot_threshold, ert_thresh_A),
    quad_B = assign_quad(structural_loss, cv_spp_site_z, sot_threshold, ert_thresh_B),
    quad_C = assign_quad(structural_loss, gini_spp_z, sot_threshold, ert_thresh_C),
    quad_D = assign_quad(structural_loss, cma_spp_site_z, sot_threshold, ert_thresh_D),
    quad_E = assign_quad(structural_loss, neg_mean_lutz_spp, sot_threshold, ert_thresh_E),
    quad_F = assign_quad(structural_loss, pc1_spp, sot_threshold, ert_thresh_F),
    quad_G = assign_quad(structural_loss, pc2_spp, sot_threshold, ert_thresh_G),
    quad_H = assign_quad(structural_loss, pc1_sppsite, sot_threshold, ert_thresh_H)
  )

# ============================================================================
# 4. HELPER: Read JPEG and make a grob
# ============================================================================

read_image_grob <- function(path) {
  if (!file.exists(path)) return(nullGrob())
  tryCatch({
    img <- readJPEG(path)
    rasterGrob(img, interpolate = TRUE)
  }, error = function(e) nullGrob())
}

# ============================================================================
# 5. CREATE IMAGE PANEL PLOTS
# ============================================================================
# For each candidate axis, make a PDF with pages organized as:
# - One page per quadrant
# - Within each page: grid of trees, sorted by species
# - Each tree cell shows: SoT image (left), ERT image (right), tree label

make_image_panel <- function(dat, quad_col, ert_col, ert_thresh, candidate_name,
                             output_file) {

  dat$quadrant <- dat[[quad_col]]
  dat$ert_value <- dat[[ert_col]]

  quad_order <- c("I: Sound", "II: Incipient", "III: Active", "IV: Cavity")
  spp_order <- c("hem", "rm", "bg", "ro")
  spp_labels <- c(hem = "T. canadensis", rm = "A. rubrum",
                  bg = "N. sylvatica", ro = "Q. rubra")

  pdf(output_file, width = 20, height = 24)

  for (q in quad_order) {
    trees_in_q <- dat %>%
      filter(quadrant == q) %>%
      mutate(species = factor(species, levels = spp_order)) %>%
      arrange(species, tree)

    n_trees <- nrow(trees_in_q)
    if (n_trees == 0) {
      plot.new()
      title(paste0(candidate_name, " — ", q, " (no trees)"), cex.main = 1.5)
      next
    }

    # Layout: up to 4 columns (each column = SoT + ERT side by side)
    ncols <- min(4, n_trees)
    nrows <- ceiling(n_trees / ncols)

    # Set up the page
    par(mfrow = c(1, 1), mar = c(0, 0, 2, 0))
    plot.new()
    title(paste0(candidate_name, "\n", q, " (n=", n_trees, ")"),
          cex.main = 1.8, font.main = 2)

    # Now draw images in a grid using grid graphics
    grid.newpage()

    # Title
    pushViewport(viewport(y = 0.97, height = 0.06, just = "top"))
    grid.text(paste0(candidate_name, "  |  ", q, "  (n=", n_trees, ")"),
              gp = gpar(fontsize = 18, fontface = "bold"))
    popViewport()

    # Image grid area
    pushViewport(viewport(y = 0.47, height = 0.88, just = "center"))

    # Each tree gets a cell with SoT on left, ERT on right, label on top
    cell_w <- 1 / ncols
    cell_h <- 1 / nrows

    for (idx in seq_len(n_trees)) {
      row_idx <- ceiling(idx / ncols)
      col_idx <- ((idx - 1) %% ncols) + 1

      tree_row <- trees_in_q[idx, ]

      # Cell viewport
      x_center <- (col_idx - 0.5) * cell_w
      y_center <- 1 - (row_idx - 0.5) * cell_h

      pushViewport(viewport(x = x_center, y = y_center,
                            width = cell_w * 0.95, height = cell_h * 0.95))

      # Label at top of cell
      pushViewport(viewport(y = 0.95, height = 0.10, just = "top"))
      label_txt <- paste0("Tree ", tree_row$tree, " (",
                          spp_labels[as.character(tree_row$species)], ")  |  ",
                          "SoT loss: ", tree_row$structural_loss, "%  |  ",
                          ert_col, ": ", round(tree_row$ert_value, 2))
      grid.text(label_txt, gp = gpar(fontsize = 8, fontface = "bold"))
      popViewport()

      # SoT image (left half)
      pushViewport(viewport(x = 0.25, y = 0.45, width = 0.48, height = 0.80))
      if (file.exists(tree_row$sot_path)) {
        tryCatch({
          img <- readJPEG(tree_row$sot_path)
          grid.raster(img, interpolate = TRUE)
        }, error = function(e) {
          grid.text("SoT\nNot found", gp = gpar(fontsize = 7, col = "grey50"))
        })
      } else {
        grid.text("SoT\nNot found", gp = gpar(fontsize = 7, col = "grey50"))
      }
      # Label below image
      pushViewport(viewport(y = -0.02, height = 0.05, just = "top"))
      grid.text("SoT", gp = gpar(fontsize = 7, col = "grey40"))
      popViewport()
      popViewport()

      # ERT image (right half)
      pushViewport(viewport(x = 0.75, y = 0.45, width = 0.48, height = 0.80))
      if (file.exists(tree_row$ert_norm_path)) {
        tryCatch({
          img <- readJPEG(tree_row$ert_norm_path)
          grid.raster(img, interpolate = TRUE)
        }, error = function(e) {
          grid.text("ERT\nNot found", gp = gpar(fontsize = 7, col = "grey50"))
        })
      } else {
        grid.text("ERT\nNot found", gp = gpar(fontsize = 7, col = "grey50"))
      }
      pushViewport(viewport(y = -0.02, height = 0.05, just = "top"))
      grid.text("ERT (norm)", gp = gpar(fontsize = 7, col = "grey40"))
      popViewport()
      popViewport()

      # Border around cell
      grid.rect(gp = gpar(col = "grey70", fill = NA, lwd = 0.5))

      popViewport()  # cell viewport
    }

    popViewport()  # image grid area
  }

  dev.off()
  cat("Saved:", output_file, "\n")
}

# ============================================================================
# 6. GENERATE PANELS FOR EACH CANDIDATE
# ============================================================================

cat("\n--- Generating image panels ---\n\n")

# A: neg_median, species×site z-norm (moisture LEVEL, high=wet)
make_image_panel(dat, "quad_A", "neg_median_spp_site_z", ert_thresh_A,
                 "A: -Median Resistivity (spp×site z-norm) [high=wet]",
                 "output/panels/panel_A_neg_median_spp_site_znorm.pdf")

# B: cv, species×site z-norm (moisture HETEROGENEITY)
make_image_panel(dat, "quad_B", "cv_spp_site_z", ert_thresh_B,
                 "B: CV of Resistivity (spp×site z-norm) [high=heterogeneous]",
                 "output/panels/panel_B_cv_spp_site_znorm.pdf")

# C: gini, species z-norm (heterogeneity, simpler normalization)
make_image_panel(dat, "quad_C", "gini_spp_z", ert_thresh_C,
                 "C: Gini of Resistivity (spp z-norm) [high=heterogeneous]",
                 "output/panels/panel_C_gini_spp_znorm.pdf")

# D: cma, species×site z-norm (spatial structure)
make_image_panel(dat, "quad_D", "cma_spp_site_z", ert_thresh_D,
                 "D: Center-of-Mass Anomaly (spp×site z-norm)",
                 "output/panels/panel_D_cma_spp_site_znorm.pdf")

# E: neg_mean, Lutz species baseline (moisture level, Lutz-style, high=wet)
make_image_panel(dat, "quad_E", "neg_mean_lutz_spp", ert_thresh_E,
                 "E: -Mean Resistivity (Lutz spp baseline) [high=wet]",
                 "output/panels/panel_E_neg_mean_lutz_spp.pdf")

# F: PC1 from species-normalized PCA (composite axis)
make_image_panel(dat, "quad_F", "pc1_spp", ert_thresh_F,
                 "F: PC1 (species-normalized PCA) [high=anomalous]",
                 "output/panels/panel_F_pc1_spp_norm.pdf")

# G: PC2 from species-normalized PCA
make_image_panel(dat, "quad_G", "pc2_spp", ert_thresh_G,
                 "G: PC2 (species-normalized PCA)",
                 "output/panels/panel_G_pc2_spp_norm.pdf")

# H: PC1 from species×site normalized PCA
make_image_panel(dat, "quad_H", "pc1_sppsite", ert_thresh_H,
                 "H: PC1 (spp×site normalized PCA) [high=anomalous]",
                 "output/panels/panel_H_pc1_sppsite_norm.pdf")

# ============================================================================
# 7. SPECIES-FACETED PHASE DIAGRAMS
# ============================================================================
# For each top candidate, make a scatter plot faceted by species so you can
# see if the quadrant assignments make sense within each species

cat("\n--- Generating species-faceted phase diagrams ---\n\n")

spp_labels_expr <- c(
  hem = "T. canadensis (hemlock)",
  rm  = "A. rubrum (red maple)",
  bg  = "N. sylvatica (black gum)",
  ro  = "Q. rubra (red oak)"
)

quad_colors <- c(
  "I: Sound"      = "#2ca02c",
  "II: Incipient" = "#ff7f0e",
  "III: Active"   = "#d62728",
  "IV: Cavity"    = "#7f7f7f"
)

make_species_phase_plot <- function(dat, quad_col, ert_col, ert_thresh,
                                    sot_thresh, title_str) {
  dat$quadrant <- dat[[quad_col]]
  dat$ert_value <- dat[[ert_col]]
  dat$species_label <- spp_labels_expr[as.character(dat$species)]

  ggplot(dat, aes(x = ert_value, y = structural_loss)) +
    annotate("rect",
             xmin = -Inf, xmax = ert_thresh, ymin = -Inf, ymax = sot_thresh,
             fill = "#2ca02c", alpha = 0.06) +
    annotate("rect",
             xmin = ert_thresh, xmax = Inf, ymin = -Inf, ymax = sot_thresh,
             fill = "#ff7f0e", alpha = 0.06) +
    annotate("rect",
             xmin = ert_thresh, xmax = Inf, ymin = sot_thresh, ymax = Inf,
             fill = "#d62728", alpha = 0.06) +
    annotate("rect",
             xmin = -Inf, xmax = ert_thresh, ymin = sot_thresh, ymax = Inf,
             fill = "#7f7f7f", alpha = 0.06) +
    geom_hline(yintercept = sot_thresh, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = ert_thresh, linetype = "dashed", color = "grey50") +
    # Training points: filled circles; Validation: open triangles with black border
    geom_point(data = ~ filter(., dataset == "training"),
               aes(color = quadrant), size = 3, alpha = 0.8, shape = 16) +
    geom_point(data = ~ filter(., dataset == "validation"),
               aes(color = quadrant), size = 4, alpha = 0.9, shape = 17, stroke = 1) +
    geom_text(aes(label = tree), size = 2.5, vjust = -1, color = "grey20") +
    scale_color_manual(values = quad_colors) +
    facet_wrap(~ species_label, scales = "free_x") +
    labs(title = title_str,
         subtitle = "Circles = training (57 trees), Triangles = validation hemlock (12 trees)",
         x = ert_col,
         y = "Damaged Wood (%)") +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      strip.text = element_text(size = 10, face = "italic"),
      legend.position = "bottom"
    )
}

pdf("output/figures/phase_diagrams_by_species.pdf", width = 14, height = 10)

print(make_species_phase_plot(dat, "quad_A", "neg_median_spp_site_z", ert_thresh_A,
                              sot_threshold, "A: -Median Resistivity (spp×site z-norm) [high=wet]"))
print(make_species_phase_plot(dat, "quad_B", "cv_spp_site_z", ert_thresh_B,
                              sot_threshold, "B: CV (spp×site z-norm) [high=heterogeneous]"))
print(make_species_phase_plot(dat, "quad_D", "cma_spp_site_z", ert_thresh_D,
                              sot_threshold, "D: CMA (spp×site z-norm)"))
print(make_species_phase_plot(dat, "quad_E", "neg_mean_lutz_spp", ert_thresh_E,
                              sot_threshold, "E: -Mean Resistivity (Lutz spp baseline) [high=wet]"))
print(make_species_phase_plot(dat, "quad_F", "pc1_spp", ert_thresh_F,
                              sot_threshold, "F: PC1 (species-normalized PCA) [high=anomalous]"))
print(make_species_phase_plot(dat, "quad_G", "pc2_spp", ert_thresh_G,
                              sot_threshold, "G: PC2 (species-normalized PCA)"))
print(make_species_phase_plot(dat, "quad_H", "pc1_sppsite", ert_thresh_H,
                              sot_threshold, "H: PC1 (spp×site normalized PCA) [high=anomalous]"))

dev.off()
cat("Saved: phase_diagrams_by_species.pdf\n")

cat("\n=== ALL PANELS DONE ===\n")
