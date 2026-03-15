library(tidyverse)
library(cluster)      # silhouette
library(mclust)       # GMM for threshold detection
library(factoextra)   # fviz_silhouette if needed
library(patchwork)    # combining plots

# ============================================================================
# 1. DATA PREP
# ============================================================================

tree_info <- read_csv("Tree_ID_info.csv", show_col_types = FALSE)
ert_data  <- read_csv("ERT_application_results.csv", show_col_types = FALSE)

# Add site designation
tree_info <- tree_info %>%
  mutate(site = case_when(
    plot %in% c("C1", "C2", "D1", "D2", "E1", "E5", "F1", "G1", "H1") ~ "EMS",
    TRUE ~ "BGS"
  ))

# Merge — this is the TRAINING set (axes + thresholds defined from these 57 trees)
dat <- tree_info %>%
  inner_join(ert_data, by = "tree")

# SoT axis: use percent_damaged directly (not 100 - percent_solid_wood)
# percent_damaged excludes edge artifacts that inflate non-solid-wood area
# Also create negated level metrics: higher = wetter (more moisture)
# Raw ERT pixel values: high = high resistivity = DRY
# So we negate mean/median so that high values = high moisture anomaly
dat <- dat %>%
  mutate(structural_loss = percent_damaged,
         abs_cma = abs(cma),
         abs_radgrad = abs(radialgradient),
         neg_mean   = -mean,
         neg_median = -median,
         neg_radgrad = -radialgradient)

cat("Training dataset:", nrow(dat), "trees\n")
cat("Species:", paste(unique(dat$species), collapse = ", "), "\n")
cat("Sites:", paste(unique(dat$site), collapse = ", "), "\n\n")

# ============================================================================
# 2. DEFINE SOUND-WOOD BASELINE POPULATION
# ============================================================================
# Trees with no detected damage for Lutz-style normalization
# Using percent_damaged == 0 (no internal damage detected by SoT)
sound_trees <- dat %>% filter(percent_damaged == 0)
cat("Sound-wood baseline trees (0% damaged):", nrow(sound_trees), "of", nrow(dat), "\n\n")

# ============================================================================
# 3. GENERATE ALL CANDIDATE ERT AXES
# ============================================================================

# Raw ERT metrics to consider
# Level metrics are NEGATED so high = wet (high moisture anomaly)
# Heterogeneity metrics (sd, cv, gini, entropy) are NOT negated (high = more variable)
# Spatial metrics: cma and abs variants kept as-is; neg_radgrad added
raw_metrics <- c("neg_mean", "neg_median", "mean", "median",
                 "sd", "cv", "gini", "entropy",
                 "cma", "abs_cma", "radialgradient", "neg_radgrad",
                 "abs_radgrad")

# --- 3a. Raw metrics ---
raw_long <- dat %>%
  select(tree, species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "value") %>%
  mutate(transform = "raw")

# --- 3b. Species-normalized (z-score within species) ---
species_norm <- dat %>%
  select(tree, species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  group_by(species, metric) %>%
  mutate(value = (raw_value - mean(raw_value, na.rm = TRUE)) /
           sd(raw_value, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-raw_value) %>%
  mutate(transform = "species_znorm")

# --- 3c. Site-normalized (z-score within site) ---
site_norm <- dat %>%
  select(tree, species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  group_by(site, metric) %>%
  mutate(value = (raw_value - mean(raw_value, na.rm = TRUE)) /
           sd(raw_value, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-raw_value) %>%
  mutate(transform = "site_znorm")

# --- 3d. Species+site normalized ---
spp_site_norm <- dat %>%
  select(tree, species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  group_by(species, site, metric) %>%
  mutate(value = if (n() > 2 & sd(raw_value, na.rm = TRUE) > 0) {
    (raw_value - mean(raw_value, na.rm = TRUE)) / sd(raw_value, na.rm = TRUE)
  } else {
    raw_value - mean(raw_value, na.rm = TRUE)  # center only if can't z-score
  }) %>%
  ungroup() %>%
  select(-raw_value) %>%
  mutate(transform = "species_site_znorm")

# --- 3e. Lutz-style: subtract sound-wood baseline (global) ---
# For each metric, compute mean of sound-wood trees, subtract it
sound_baselines <- sound_trees %>%
  select(all_of(raw_metrics)) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "baseline")

lutz_global <- dat %>%
  select(tree, species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  left_join(sound_baselines, by = "metric") %>%
  mutate(value = raw_value - baseline) %>%
  select(-raw_value, -baseline) %>%
  mutate(transform = "lutz_global")

# --- 3f. Lutz-style: subtract sound-wood baseline (per species) ---
sound_baselines_spp <- sound_trees %>%
  select(species, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  group_by(species, metric) %>%
  summarise(baseline = mean(raw_value, na.rm = TRUE), .groups = "drop")

lutz_species <- dat %>%
  select(tree, species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  left_join(sound_baselines_spp, by = c("species", "metric")) %>%
  mutate(value = raw_value - baseline) %>%
  select(-raw_value, -baseline) %>%
  mutate(transform = "lutz_species")

# --- 3g. Lutz-style: subtract sound-wood baseline (per species+site) ---
sound_baselines_spp_site <- sound_trees %>%
  select(species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  group_by(species, site, metric) %>%
  summarise(baseline = mean(raw_value, na.rm = TRUE),
            n_baseline = n(), .groups = "drop")

lutz_spp_site <- dat %>%
  select(tree, species, site, all_of(raw_metrics)) %>%
  pivot_longer(cols = all_of(raw_metrics),
               names_to = "metric", values_to = "raw_value") %>%
  left_join(sound_baselines_spp_site, by = c("species", "site", "metric")) %>%
  mutate(value = if_else(!is.na(baseline) & n_baseline >= 2,
                         raw_value - baseline,
                         NA_real_)) %>%
  select(-raw_value, -baseline, -n_baseline) %>%
  mutate(transform = "lutz_species_site")

# --- Combine all single-metric candidates ---
all_candidates <- bind_rows(
  raw_long, species_norm, site_norm, spp_site_norm,
  lutz_global, lutz_species, lutz_spp_site
) %>%
  mutate(axis_id = paste0(metric, "__", transform))

cat("Total single-metric candidates:", n_distinct(all_candidates$axis_id), "\n")

# ============================================================================
# 4. COMPOSITE AXES: PRODUCTS, PCA
# ============================================================================

# --- 4a. Products of heterogeneity × spatial structure ---
# Use neg_mean for products so high = wet × heterogeneous
products <- dat %>%
  mutate(
    cv_x_abscma       = cv * abs_cma,
    gini_x_abscma     = gini * abs_cma,
    cv_x_absradgrad    = cv * abs_radgrad,
    gini_x_absradgrad  = gini * abs_radgrad,
    sd_x_abscma        = sd * abs_cma,
    cv_x_negmean       = cv * neg_mean,
    gini_x_negmean     = gini * neg_mean,
    sd_x_cv            = sd * cv,
    negmean_x_abscma   = neg_mean * abs_cma
  ) %>%
  select(tree, species, site,
         cv_x_abscma, gini_x_abscma, cv_x_absradgrad,
         gini_x_absradgrad, sd_x_abscma, cv_x_negmean,
         gini_x_negmean, sd_x_cv, negmean_x_abscma) %>%
  pivot_longer(cols = -c(tree, species, site),
               names_to = "metric", values_to = "value") %>%
  mutate(transform = "product",
         axis_id = paste0(metric, "__product"))

# --- 4b. Species-normalized products ---
products_spp_norm <- products %>%
  group_by(species, metric) %>%
  mutate(value = (value - mean(value, na.rm = TRUE)) /
           sd(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(transform = "product_species_znorm",
         axis_id = paste0(metric, "__product_species_znorm"))

# --- 4c. PCA on raw ERT metrics ---
pca_raw_mat <- dat %>%
  select(all_of(c("mean", "median", "sd", "cv", "gini", "entropy",
                   "abs_cma", "abs_radgrad"))) %>%
  scale()

pca_raw <- prcomp(pca_raw_mat, center = FALSE, scale. = FALSE)
cat("\nPCA on raw ERT metrics - variance explained:\n")
print(summary(pca_raw)$importance[, 1:4])

pca_raw_scores <- as_tibble(pca_raw$x[, 1:4]) %>%
  mutate(tree = dat$tree, species = dat$species, site = dat$site) %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "metric", values_to = "value") %>%
  mutate(transform = "pca_raw",
         axis_id = paste0(metric, "__pca_raw"))

# --- 4d. PCA on species-normalized metrics ---
pca_spp_mat <- dat %>%
  select(tree, species, all_of(c("mean", "median", "sd", "cv", "gini",
                                  "entropy", "abs_cma", "abs_radgrad"))) %>%
  group_by(species) %>%
  mutate(across(-tree, ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE))) %>%
  ungroup() %>%
  select(-tree, -species) %>%
  as.matrix()

# Handle any NaN from constant columns within species
pca_spp_mat[is.nan(pca_spp_mat)] <- 0

pca_spp <- prcomp(pca_spp_mat, center = FALSE, scale. = FALSE)
cat("\nPCA on species-normalized ERT metrics - variance explained:\n")
print(summary(pca_spp)$importance[, 1:4])

pca_spp_scores <- as_tibble(pca_spp$x[, 1:4]) %>%
  mutate(tree = dat$tree, species = dat$species, site = dat$site) %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "metric", values_to = "value") %>%
  mutate(transform = "pca_species_norm",
         axis_id = paste0(metric, "__pca_species_norm"))

# --- 4e. Lutz-style products ---
lutz_product_dat <- dat %>%
  left_join(sound_baselines_spp %>%
              pivot_wider(names_from = metric, values_from = baseline,
                          names_prefix = "bl_"),
            by = "species") %>%
  mutate(
    lutz_cv       = cv - bl_cv,
    lutz_gini     = gini - bl_gini,
    lutz_abscma   = abs_cma - bl_abs_cma,
    lutz_absrad   = abs_radgrad - bl_abs_radgrad,
    lutz_cv_x_abscma   = lutz_cv * lutz_abscma,
    lutz_gini_x_abscma = lutz_gini * lutz_abscma
  ) %>%
  select(tree, species, site,
         lutz_cv_x_abscma, lutz_gini_x_abscma) %>%
  pivot_longer(cols = -c(tree, species, site),
               names_to = "metric", values_to = "value") %>%
  mutate(transform = "lutz_product",
         axis_id = paste0(metric, "__lutz_product"))

# --- Combine all candidates ---
all_candidates <- bind_rows(
  all_candidates,
  products, products_spp_norm,
  pca_raw_scores, pca_spp_scores,
  lutz_product_dat
)

cat("\nTotal candidate ERT axes:", n_distinct(all_candidates$axis_id), "\n\n")

# ============================================================================
# 5. DEFINE SoT AXIS AND THRESHOLD
# ============================================================================

sot_vals <- sort(unique(dat$structural_loss))
cat("Structural loss (percent_damaged) values:", head(sot_vals, 20), "\n")
cat("Trees with 0% damage:", sum(dat$structural_loss == 0), "\n")
cat("Trees with >0% damage:", sum(dat$structural_loss > 0), "\n\n")

# Threshold: any detected damage (percent_damaged > 0) is meaningful
# since percent_damaged specifically identifies internal structural decay
# Use a small threshold just above 0 to separate sound from decaying
sot_threshold <- 1  # trees with >=1% detected damage are "decaying"

# ============================================================================
# 6. THRESHOLD METHODS FOR ERT AXIS (independent of SoT)
# ============================================================================

find_ert_thresholds <- function(values, method = "all") {
  values <- values[!is.na(values)]
  thresholds <- list()

  # Method 1: Median
  thresholds$median <- median(values)

  # Method 2: Mean
  thresholds$mean <- mean(values)

  # Method 3: Otsu's method (maximize between-class variance)
  # Try candidate thresholds at each observed value
  otsu_best <- -Inf
  otsu_thresh <- median(values)
  candidates <- quantile(values, probs = seq(0.1, 0.9, by = 0.02))
  for (t in candidates) {
    w0 <- mean(values <= t)
    w1 <- 1 - w0
    if (w0 > 0 & w1 > 0) {
      mu0 <- mean(values[values <= t])
      mu1 <- mean(values[values > t])
      between_var <- w0 * w1 * (mu0 - mu1)^2
      if (between_var > otsu_best) {
        otsu_best <- between_var
        otsu_thresh <- t
      }
    }
  }
  thresholds$otsu <- otsu_thresh

  # Method 4: GMM with 2 components
  tryCatch({
    gmm <- Mclust(values, G = 2, verbose = FALSE)
    if (!is.null(gmm)) {
      # Threshold at the crossover point between the two Gaussians
      means_sorted <- sort(gmm$parameters$mean)
      thresholds$gmm <- mean(means_sorted)  # midpoint of two means
    }
  }, error = function(e) {})

  # Method 5: 75th percentile
  thresholds$q75 <- quantile(values, 0.75)

  # Method 6: Mean + 1 SD
  thresholds$mean_plus_sd <- mean(values) + sd(values)

  return(thresholds)
}

# ============================================================================
# 7. EVALUATE EACH CANDIDATE AXIS × THRESHOLD COMBINATION
# ============================================================================

# Prepare structural loss for each tree
sot_by_tree <- dat %>% select(tree, structural_loss, species, site,
                               percent_solid_wood, percent_damaged, decay)

# Function to assign quadrants given SoT and ERT values + thresholds
assign_quadrants <- function(sot, ert, sot_thresh, ert_thresh) {
  case_when(
    sot <= sot_thresh & ert <= ert_thresh ~ "I_sound",
    sot <= sot_thresh & ert >  ert_thresh ~ "II_incipient",
    sot >  sot_thresh & ert >  ert_thresh ~ "III_active",
    sot >  sot_thresh & ert <= ert_thresh ~ "IV_cavity",
    TRUE ~ NA_character_
  )
}

# Function to compute silhouette score for a quadrant assignment
compute_silhouette <- function(sot_scaled, ert_scaled, quadrants) {
  valid <- !is.na(quadrants) & !is.na(sot_scaled) & !is.na(ert_scaled)
  if (sum(valid) < 10) return(NA_real_)

  quads <- as.integer(factor(quadrants[valid]))
  n_groups <- length(unique(quads))
  if (n_groups < 2) return(NA_real_)

  coords <- cbind(sot_scaled[valid], ert_scaled[valid])
  d <- dist(coords)
  sil <- silhouette(quads, d)
  mean(sil[, "sil_width"])
}

# Function to compute balance metric (how evenly distributed across quadrants)
compute_balance <- function(quadrants) {
  quadrants <- quadrants[!is.na(quadrants)]
  n <- length(quadrants)
  if (n == 0) return(NA_real_)
  counts <- table(factor(quadrants, levels = c("I_sound", "II_incipient",
                                                "III_active", "IV_cavity")))
  props <- as.numeric(counts) / n
  # Entropy of distribution (max = log(4) for uniform)
  props <- props[props > 0]
  -sum(props * log(props)) / log(4)
}

# Function to check ecological coherence
# ERT axis: high values = high moisture anomaly (level metrics are negated)
# Expected trajectory through phase space:
#   I (sound):     low SoT, low ERT moisture anomaly
#   II (incipient): low SoT, elevated ERT moisture (decay starting, accumulating moisture)
#   III (active):   high SoT, high ERT moisture (structural loss + wet)
#   IV (cavity):    high SoT, low ERT moisture (wood gone, air-filled = dry)
compute_coherence <- function(sot, ert, quadrants) {
  df <- tibble(sot = sot, ert = ert, q = quadrants) %>%
    filter(!is.na(q))

  if (nrow(df) < 4) return(list(score = NA_real_, details = NULL))

  means <- df %>%
    group_by(q) %>%
    summarise(mean_sot = mean(sot), mean_ert = mean(ert), n = n(),
              .groups = "drop")

  quad_means <- setNames(
    means$mean_sot[match(c("I_sound", "II_incipient", "III_active", "IV_cavity"),
                         means$q)],
    c("I", "II", "III", "IV")
  )

  ert_means <- setNames(
    means$mean_ert[match(c("I_sound", "II_incipient", "III_active", "IV_cavity"),
                         means$q)],
    c("I", "II", "III", "IV")
  )

  score <- 0
  checks <- 0

  # SoT expectations: I and II should have less damage than III and IV
  if (!is.na(quad_means["I"]) & !is.na(quad_means["III"])) {
    if (quad_means["I"] < quad_means["III"]) score <- score + 1
    checks <- checks + 1
  }
  if (!is.na(quad_means["II"]) & !is.na(quad_means["III"])) {
    if (quad_means["II"] < quad_means["III"]) score <- score + 1
    checks <- checks + 1
  }
  if (!is.na(quad_means["I"]) & !is.na(quad_means["IV"])) {
    if (quad_means["I"] < quad_means["IV"]) score <- score + 1
    checks <- checks + 1
  }

  # ERT moisture expectations (high = wet):
  # II > I (incipient decay has more moisture than sound)
  # III > I (active decay has more moisture than sound)
  # III > IV (active decay is wetter than cavity)
  # IV ~ I or IV < III (cavity is dry, similar to sound on moisture axis)
  if (!is.na(ert_means["II"]) & !is.na(ert_means["I"])) {
    if (ert_means["II"] > ert_means["I"]) score <- score + 1
    checks <- checks + 1
  }
  if (!is.na(ert_means["III"]) & !is.na(ert_means["I"])) {
    if (ert_means["III"] > ert_means["I"]) score <- score + 1
    checks <- checks + 1
  }
  if (!is.na(ert_means["III"]) & !is.na(ert_means["IV"])) {
    if (ert_means["III"] > ert_means["IV"]) score <- score + 1
    checks <- checks + 1
  }

  coherence_score <- if (checks > 0) score / checks else NA_real_

  list(score = coherence_score, details = means)
}

# --- Run the evaluation ---
results <- tibble()

axis_ids <- unique(all_candidates$axis_id)
cat("Evaluating", length(axis_ids), "candidate ERT axes...\n")

for (aid in axis_ids) {

  cand <- all_candidates %>%
    filter(axis_id == aid) %>%
    select(tree, value)

  eval_df <- sot_by_tree %>%
    inner_join(cand, by = "tree")

  if (nrow(eval_df) < 10 || all(is.na(eval_df$value))) next

  # Get thresholds
  thresholds <- find_ert_thresholds(eval_df$value)

  for (thresh_name in names(thresholds)) {
    ert_thresh <- thresholds[[thresh_name]]
    if (is.na(ert_thresh)) next

    # Assign quadrants
    eval_df$quadrant <- assign_quadrants(
      eval_df$structural_loss, eval_df$value,
      sot_threshold, ert_thresh
    )

    # Skip if fewer than 2 quadrants populated
    n_quads <- n_distinct(eval_df$quadrant[!is.na(eval_df$quadrant)])
    if (n_quads < 2) next

    # Scale for silhouette
    sot_s <- scale(eval_df$structural_loss)[, 1]
    ert_s <- scale(eval_df$value)[, 1]

    sil <- compute_silhouette(sot_s, ert_s, eval_df$quadrant)
    bal <- compute_balance(eval_df$quadrant)
    coh <- compute_coherence(eval_df$structural_loss, eval_df$value,
                             eval_df$quadrant)

    quad_counts <- table(factor(eval_df$quadrant,
                                levels = c("I_sound", "II_incipient",
                                           "III_active", "IV_cavity")))

    results <- bind_rows(results, tibble(
      axis_id = aid,
      threshold_method = thresh_name,
      ert_threshold = ert_thresh,
      silhouette = sil,
      balance = bal,
      coherence = coh$score,
      n_I   = as.integer(quad_counts["I_sound"]),
      n_II  = as.integer(quad_counts["II_incipient"]),
      n_III = as.integer(quad_counts["III_active"]),
      n_IV  = as.integer(quad_counts["IV_cavity"]),
      n_total = nrow(eval_df)
    ))
  }
}

cat("Total combinations evaluated:", nrow(results), "\n\n")

# ============================================================================
# 8. RANK RESULTS
# ============================================================================

# Composite score: weight silhouette, balance, and coherence
# Penalize if any quadrant has 0 trees (missing a phase)
results <- results %>%
  mutate(
    n_empty_quads = (n_I == 0) + (n_II == 0) + (n_III == 0) + (n_IV == 0),
    has_all_quads = n_empty_quads == 0,
    # Composite score (higher = better)
    composite = case_when(
      is.na(silhouette) | is.na(coherence) ~ NA_real_,
      TRUE ~ 0.35 * silhouette + 0.30 * coherence + 0.25 * balance +
        0.10 * as.numeric(has_all_quads)
    )
  ) %>%
  arrange(desc(composite))

cat("=== TOP 30 CANDIDATE ERT AXES ===\n\n")
top30 <- results %>%
  filter(!is.na(composite)) %>%
  head(30) %>%
  select(axis_id, threshold_method, silhouette, balance, coherence,
         composite, n_I, n_II, n_III, n_IV)

print(top30, n = 30, width = 200)

# Also show top results that have all 4 quadrants populated
cat("\n=== TOP 20 WITH ALL 4 QUADRANTS POPULATED ===\n\n")
top20_full <- results %>%
  filter(has_all_quads, !is.na(composite)) %>%
  head(20) %>%
  select(axis_id, threshold_method, silhouette, balance, coherence,
         composite, n_I, n_II, n_III, n_IV)

print(top20_full, n = 20, width = 200)

# ============================================================================
# 9. PHASE DIAGRAM PLOTS FOR TOP CANDIDATES
# ============================================================================

# Helper function to make a phase diagram plot
make_phase_plot <- function(dat_merged, ert_col_name, ert_thresh, sot_thresh,
                            title_str, color_by_species = TRUE) {

  dat_merged$quadrant <- assign_quadrants(
    dat_merged$structural_loss, dat_merged$ert_value,
    sot_thresh, ert_thresh
  )

  quad_colors <- c(
    "I_sound"      = "#2ca02c",
    "II_incipient" = "#ff7f0e",
    "III_active"   = "#d62728",
    "IV_cavity"    = "#7f7f7f"
  )

  spp_shapes <- c(hem = 16, rm = 17, bg = 15, ro = 18)

  p <- ggplot(dat_merged, aes(x = ert_value, y = structural_loss)) +
    # Quadrant shading
    annotate("rect",
             xmin = -Inf, xmax = ert_thresh, ymin = -Inf, ymax = sot_thresh,
             fill = "#2ca02c", alpha = 0.08) +
    annotate("rect",
             xmin = ert_thresh, xmax = Inf, ymin = -Inf, ymax = sot_thresh,
             fill = "#ff7f0e", alpha = 0.08) +
    annotate("rect",
             xmin = ert_thresh, xmax = Inf, ymin = sot_thresh, ymax = Inf,
             fill = "#d62728", alpha = 0.08) +
    annotate("rect",
             xmin = -Inf, xmax = ert_thresh, ymin = sot_thresh, ymax = Inf,
             fill = "#7f7f7f", alpha = 0.08) +
    # Threshold lines
    geom_hline(yintercept = sot_thresh, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = ert_thresh, linetype = "dashed", color = "grey40") +
    # Quadrant labels
    annotate("text", x = -Inf, y = -Inf, label = "I: Sound",
             hjust = -0.1, vjust = -0.5, color = "#2ca02c", fontface = "bold", size = 3) +
    annotate("text", x = Inf, y = -Inf, label = "II: Incipient",
             hjust = 1.1, vjust = -0.5, color = "#ff7f0e", fontface = "bold", size = 3) +
    annotate("text", x = Inf, y = Inf, label = "III: Active",
             hjust = 1.1, vjust = 1.5, color = "#d62728", fontface = "bold", size = 3) +
    annotate("text", x = -Inf, y = Inf, label = "IV: Cavity",
             hjust = -0.1, vjust = 1.5, color = "#7f7f7f", fontface = "bold", size = 3)

  if (color_by_species) {
    spp_colors <- c(hem = "palegreen4", rm = "darkorange2",
                    bg = "steelblue3", ro = "orchid3")
    p <- p +
      geom_point(aes(color = species, shape = species), size = 3, alpha = 0.8) +
      scale_color_manual(values = spp_colors) +
      scale_shape_manual(values = spp_shapes)
  } else {
    p <- p +
      geom_point(aes(color = quadrant), size = 3, alpha = 0.8) +
      scale_color_manual(values = quad_colors)
  }

  p <- p +
    geom_text(aes(label = tree), size = 2, vjust = -1, color = "grey30") +
    labs(title = title_str,
         x = ert_col_name,
         y = "Damaged Wood (%)") +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(size = 9, face = "bold")
    )

  return(p)
}

# Plot top N candidates
n_to_plot <- min(20, nrow(top30))

pdf("phase_diagram_candidates.pdf", width = 14, height = 8)

for (i in seq_len(n_to_plot)) {
  row <- if (nrow(top20_full) > 0) top20_full[min(i, nrow(top20_full)), ]
         else top30[i, ]

  aid   <- row$axis_id
  tmeth <- row$threshold_method

  cand <- all_candidates %>%
    filter(axis_id == aid) %>%
    select(tree, value)

  plot_df <- sot_by_tree %>%
    inner_join(cand, by = "tree") %>%
    rename(ert_value = value)

  ert_thresh <- results %>%
    filter(axis_id == aid, threshold_method == tmeth) %>%
    pull(ert_threshold)

  title_txt <- sprintf("#%d: %s | thresh: %s (%.2f) | sil=%.3f bal=%.2f coh=%.2f",
                       i, aid, tmeth, ert_thresh,
                       row$silhouette, row$balance, row$coherence)

  p1 <- make_phase_plot(plot_df, aid, ert_thresh, sot_threshold,
                        paste0(title_txt, "\n(colored by species)"),
                        color_by_species = TRUE)
  p2 <- make_phase_plot(plot_df, aid, ert_thresh, sot_threshold,
                        paste0(title_txt, "\n(colored by quadrant)"),
                        color_by_species = FALSE)

  print(p1 + p2)
}

dev.off()
cat("\nSaved phase_diagram_candidates.pdf\n")

# ============================================================================
# 10. DETAILED SUMMARY OF TOP 5
# ============================================================================

cat("\n============================================================\n")
cat("DETAILED SUMMARY OF TOP 5 CANDIDATES\n")
cat("============================================================\n\n")

top5 <- if (nrow(top20_full) >= 5) top20_full[1:5, ] else top30[1:min(5, nrow(top30)), ]

for (i in seq_len(nrow(top5))) {
  row <- top5[i, ]
  aid   <- row$axis_id
  tmeth <- row$threshold_method

  cat(sprintf("--- #%d: %s (threshold: %s) ---\n", i, aid, tmeth))
  cat(sprintf("  Silhouette: %.4f\n", row$silhouette))
  cat(sprintf("  Balance: %.4f\n", row$balance))
  cat(sprintf("  Coherence: %.4f\n", row$coherence))
  cat(sprintf("  Composite: %.4f\n", row$composite))
  cat(sprintf("  Quadrant counts: I=%d, II=%d, III=%d, IV=%d\n",
              row$n_I, row$n_II, row$n_III, row$n_IV))

  # Show which trees fall in each quadrant
  cand <- all_candidates %>%
    filter(axis_id == aid) %>%
    select(tree, value)

  eval_df <- sot_by_tree %>%
    inner_join(cand, by = "tree") %>%
    rename(ert_value = value)

  ert_thresh <- results %>%
    filter(axis_id == aid, threshold_method == tmeth) %>%
    pull(ert_threshold)

  eval_df$quadrant <- assign_quadrants(
    eval_df$structural_loss, eval_df$ert_value,
    sot_threshold, ert_thresh
  )

  for (q in c("I_sound", "II_incipient", "III_active", "IV_cavity")) {
    trees_in_q <- eval_df %>% filter(quadrant == q) %>%
      arrange(tree) %>%
      mutate(label = paste0(tree, "(", species, ")"))
    cat(sprintf("  %s: %s\n", q, paste(trees_in_q$label, collapse = ", ")))
  }
  cat("\n")
}

# ============================================================================
# 11. WRITE FULL RESULTS TABLE
# ============================================================================

write_csv(results %>% arrange(desc(composite)),
          "phase_diagram_results.csv")
cat("Full results table saved to phase_diagram_results.csv\n")

# ============================================================================
# 12. BONUS: k-means clustering (data-driven, no predefined thresholds)
# ============================================================================

cat("\n============================================================\n")
cat("BONUS: k-MEANS CLUSTERING (k=4) ON TOP ERT AXES\n")
cat("============================================================\n\n")

pdf("phase_diagram_kmeans.pdf", width = 14, height = 8)

for (i in seq_len(min(10, nrow(top5)))) {
  row <- top5[i, ]
  aid <- row$axis_id

  cand <- all_candidates %>%
    filter(axis_id == aid) %>%
    select(tree, value)

  km_df <- sot_by_tree %>%
    inner_join(cand, by = "tree") %>%
    rename(ert_value = value) %>%
    filter(!is.na(ert_value))

  # Scale both axes
  km_mat <- km_df %>%
    mutate(sot_s = scale(structural_loss)[, 1],
           ert_s = scale(ert_value)[, 1]) %>%
    select(sot_s, ert_s) %>%
    as.matrix()

  set.seed(42)
  km <- kmeans(km_mat, centers = 4, nstart = 25)
  km_df$cluster <- factor(km$cluster)

  sil_km <- silhouette(km$cluster, dist(km_mat))
  avg_sil <- mean(sil_km[, "sil_width"])

  p <- ggplot(km_df, aes(x = ert_value, y = structural_loss)) +
    geom_point(aes(color = cluster, shape = species), size = 3, alpha = 0.8) +
    geom_text(aes(label = tree), size = 2, vjust = -1, color = "grey30") +
    scale_shape_manual(values = c(hem = 16, rm = 17, bg = 15, ro = 18)) +
    labs(title = sprintf("k-means (k=4) on %s | silhouette = %.3f", aid, avg_sil),
         x = aid,
         y = "Damaged Wood (%)") +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(size = 10, face = "bold"))

  print(p)
}

dev.off()
cat("Saved phase_diagram_kmeans.pdf\n")

cat("\n=== DONE ===\n")
