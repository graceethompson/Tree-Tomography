##############################################################################
# ERT Validation Analysis
# Compares ERT image-derived metrics to core moisture measurements
##############################################################################

library(tidyverse)
library(readxl)
library(ggrepel)
library(pheatmap)
library(factoextra)
library(patchwork)
library(cowplot)
library(grid)

# Set working directory
setwd("/Users/jongewirtzman/My Drive/Research/Tomography/Tree-Tomography/Hemlock ERT")

# ── Color palette ──────────────────────────────────────────────────────────
pal_height <- c("DBH" = "#2166AC", "Lower" = "#4DAC26", "Upper" = "#D6604D")
pal_accent <- "#2166AC"

# ── Base theme ─────────────────────────────────────────────────────────────
theme_pub <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 10),
    legend.position = "bottom"
  )
theme_set(theme_pub)

# ── Metric names (consistent labels) ──────────────────────────────────────
metric_names <- c("Conductance", "Median", "SD", "CV", "Gini", "Entropy", "CMA", "RadialGradient")
# Note: Conductance is in milliSiemens (mS) = 1000/Resistance

##############################################################################
# 1. DATA LOADING & CLEANING
##############################################################################

# Load ERT results
ert_raw <- read_csv("ERT_results_2026-03-05.csv", show_col_types = FALSE)

# Parse tree_id and height from Filename
ert_raw <- ert_raw %>%
  mutate(
    base = str_remove(Filename, "\\.jpg$"),
    height = str_extract(base, "[^_]+$"),
    tree_id = str_remove(base, "_[^_]+$")
  ) %>%
  select(-base)

# Deduplicate: keep last row per tree_id × height
ert <- ert_raw %>%
  group_by(tree_id, height) %>%
  slice_tail(n = 1) %>%
  ungroup()

cat("Unique tree×height combos after dedup:", nrow(ert), "\n")
cat("Heights per tree:\n")
print(ert %>% count(tree_id, name = "n_heights"))

# Load core moisture data
mc <- read_excel("MC_Tomo_paper_Jon.xlsx", sheet = 1)
names(mc) <- c("tree_id", "moisture")

# Merge ERT with moisture
ert <- ert %>% left_join(mc, by = "tree_id")

# Add conductance (1/Resistance)
ert <- ert %>% mutate(Conductance = 1000 / Mean)  # milliSiemens (mS)

# Subsets
ert_dbh   <- ert %>% filter(height == "DBH")
ert_lower <- ert %>% filter(height == "Lower")
ert_upper <- ert %>% filter(height == "Upper")

# Average across heights per tree
ert_avg <- ert %>%
  group_by(tree_id) %>%
  summarise(across(all_of(metric_names), mean), moisture = first(moisture),
            Conductance = mean(Conductance), .groups = "drop")

cat("\nDBH trees:", nrow(ert_dbh), "\n")
cat("Lower trees:", nrow(ert_lower), "\n")
cat("Upper trees:", nrow(ert_upper), "\n")

##############################################################################
# HELPER: Scatter plot with correlation stats
##############################################################################
make_scatter <- function(data, metric, x_var = "moisture",
                         show_label = TRUE, color = pal_accent,
                         x_lab = "Core Moisture (%)", point_size = 2.5) {

  # Correlation tests
  ct_p <- cor.test(data[[x_var]], data[[metric]], method = "pearson")
  ct_s <- cor.test(data[[x_var]], data[[metric]], method = "spearman",
                   exact = FALSE)

  r_val  <- round(ct_p$estimate, 3)
  p_pear <- ct_p$p.value
  rho    <- round(ct_s$estimate, 3)
  p_spear <- ct_s$p.value

  # Format p-values
  fmt_p <- function(p) {
    if (p < 0.001) return("p < 0.001")
    paste0("p = ", formatC(p, format = "f", digits = 3))
  }

  label_text <- paste0(
    "r = ", r_val, " (", fmt_p(p_pear), ")\n",
    "\u03C1 = ", rho, " (", fmt_p(p_spear), ")"
  )

  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[metric]])) +
    geom_point(size = point_size, color = color, alpha = 0.85)

  # Add regression line if significant
  if (p_pear < 0.05) {
    p <- p + geom_smooth(method = "lm", se = TRUE, color = color,
                         fill = color, alpha = 0.15, linewidth = 0.7)
  }

  if (show_label) {
    p <- p + geom_text_repel(aes(label = tree_id), size = 2.5,
                              max.overlaps = 20, segment.color = "grey60",
                              segment.size = 0.3)
  }

  # Annotation position: top-left
  y_range <- range(data[[metric]], na.rm = TRUE)
  x_range <- range(data[[x_var]], na.rm = TRUE)

  p <- p +
    annotate("text",
             x = x_range[1] + diff(x_range) * 0.02,
             y = y_range[2] - diff(y_range) * 0.02,
             label = label_text, hjust = 0, vjust = 1, size = 2.8,
             color = "grey30") +
    labs(x = x_lab, y = metric)

  return(p)
}

##############################################################################
# FIGURE 1: Scatter plots — ERT metrics vs Core Moisture at DBH
##############################################################################
cat("\n── Figure 1: DBH scatter plots ──\n")

plots_dbh <- map(metric_names, ~ make_scatter(ert_dbh, .x))
fig1 <- wrap_plots(plots_dbh, ncol = 4) +
  plot_annotation(title = "ERT Metrics vs. Core Moisture (DBH)")

ggsave("fig_scatter_dbh.pdf", fig1, width = 14, height = 8)
cat("Saved fig_scatter_dbh.pdf\n")

##############################################################################
# FIGURE 2: Scatter plots by height (DBH, Lower, Upper)
##############################################################################
cat("\n── Figure 2: Scatter plots by height ──\n")

make_scatter_height <- function(data, metric, height_label, color) {
  n <- nrow(data)
  if (n < 3) return(ggplot() + theme_void() + ggtitle(paste(height_label, "-", metric)))
  make_scatter(data, metric, color = color, point_size = 2) +
    ggtitle(height_label)
}

height_sets <- list(
  list(data = ert_dbh,   label = "DBH",   color = pal_height["DBH"]),
  list(data = ert_lower, label = "Lower", color = pal_height["Lower"]),
  list(data = ert_upper, label = "Upper", color = pal_height["Upper"])
)

all_panels <- list()
for (hs in height_sets) {
  row_plots <- map(metric_names, ~ make_scatter_height(hs$data, .x, hs$label, hs$color))
  all_panels <- c(all_panels, row_plots)
}

fig2 <- wrap_plots(all_panels, ncol = 8, nrow = 3) +
  plot_annotation(title = "ERT Metrics vs. Core Moisture by Measurement Height")

ggsave("fig_scatter_by_height.pdf", fig2, width = 22, height = 12)
cat("Saved fig_scatter_by_height.pdf\n")

##############################################################################
# FIGURE 3: Average-across-heights ERT metrics vs Core Moisture
##############################################################################
cat("\n── Figure 3: Average scatter plots ──\n")

plots_avg <- map(metric_names, ~ make_scatter(ert_avg, .x))
fig3 <- wrap_plots(plots_avg, ncol = 4) +
  plot_annotation(title = "Average ERT Metrics (Across Heights) vs. Core Moisture")

ggsave("fig_scatter_avg.pdf", fig3, width = 14, height = 8)
cat("Saved fig_scatter_avg.pdf\n")

##############################################################################
# FIGURE 4: Correlation heatmap
##############################################################################
cat("\n── Figure 4: Correlation heatmap ──\n")

compute_cors <- function(data, label) {
  results <- tibble()
  for (m in metric_names) {
    if (nrow(data) < 3) {
      results <- bind_rows(results, tibble(
        metric = m, group = label,
        pearson_r = NA, pearson_p = NA,
        spearman_rho = NA, spearman_p = NA
      ))
      next
    }
    cp <- cor.test(data$moisture, data[[m]], method = "pearson")
    cs <- cor.test(data$moisture, data[[m]], method = "spearman", exact = FALSE)
    results <- bind_rows(results, tibble(
      metric = m, group = label,
      pearson_r = cp$estimate, pearson_p = cp$p.value,
      spearman_rho = cs$estimate, spearman_p = cs$p.value
    ))
  }
  results
}

cor_results <- bind_rows(
  compute_cors(ert_dbh, "DBH"),
  compute_cors(ert_lower, "Lower"),
  compute_cors(ert_upper, "Upper"),
  compute_cors(ert_avg, "Average")
)

# Build matrices for heatmap
# Pearson matrix
pear_wide <- cor_results %>%
  select(metric, group, pearson_r) %>%
  pivot_wider(names_from = group, values_from = pearson_r, names_prefix = "Pearson\n") %>%
  column_to_rownames("metric")

spear_wide <- cor_results %>%
  select(metric, group, spearman_rho) %>%
  pivot_wider(names_from = group, values_from = spearman_rho, names_prefix = "Spearman\n") %>%
  column_to_rownames("metric")

heat_mat <- as.matrix(cbind(pear_wide, spear_wide))

# Significance stars
pear_p <- cor_results %>%
  select(metric, group, pearson_p) %>%
  pivot_wider(names_from = group, values_from = pearson_p, names_prefix = "Pearson\n") %>%
  column_to_rownames("metric")

spear_p <- cor_results %>%
  select(metric, group, spearman_p) %>%
  pivot_wider(names_from = group, values_from = spearman_p, names_prefix = "Spearman\n") %>%
  column_to_rownames("metric")

p_mat <- as.matrix(cbind(pear_p, spear_p))

# Display text: r/rho value + star if significant
disp_mat <- matrix("", nrow = nrow(heat_mat), ncol = ncol(heat_mat))
for (i in seq_len(nrow(heat_mat))) {
  for (j in seq_len(ncol(heat_mat))) {
    val <- heat_mat[i, j]
    pval <- p_mat[i, j]
    if (is.na(val)) {
      disp_mat[i, j] <- "—"
    } else {
      star <- ifelse(!is.na(pval) & pval < 0.05, "*", "")
      disp_mat[i, j] <- paste0(formatC(val, format = "f", digits = 2), star)
    }
  }
}

pdf("fig_correlation_heatmap.pdf", width = 10, height = 6)
pheatmap(
  heat_mat,
  display_numbers = disp_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "grey90",
  fontsize = 10,
  fontsize_number = 9,
  main = "Correlations: ERT Metrics vs. Core Moisture\n(* p < 0.05)",
  angle_col = 45
)
dev.off()
cat("Saved fig_correlation_heatmap.pdf\n")

##############################################################################
# FIGURE 5a: PCA biplot — DBH only
##############################################################################
cat("\n── Figure 5: PCA biplots ──\n")

make_pca_biplot <- function(data, title_text, filename) {
  pca_data <- data %>% select(all_of(metric_names))
  pca_res <- prcomp(pca_data, scale. = TRUE)

  # Scores
  scores <- as.data.frame(pca_res$x[, 1:2])
  scores$tree_id <- data$tree_id
  scores$moisture <- data$moisture
  if ("height" %in% names(data)) scores$height <- data$height

  # Loadings (scaled for visibility)
  loadings <- as.data.frame(pca_res$rotation[, 1:2])
  loadings$metric <- rownames(loadings)

  # Variance explained
  ve <- summary(pca_res)$importance[2, 1:2] * 100

  # Scale loadings for plotting
  scale_factor <- max(abs(scores$PC1), abs(scores$PC2)) * 0.8 /
    max(abs(loadings$PC1), abs(loadings$PC2))
  loadings$PC1 <- loadings$PC1 * scale_factor
  loadings$PC2 <- loadings$PC2 * scale_factor

  p <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = moisture), size = 3.5) +
    scale_color_viridis_c(name = "Core\nMoisture (%)") +
    geom_text_repel(aes(label = tree_id), size = 2.8, max.overlaps = 20) +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "grey40", linewidth = 0.5) +
    geom_text(data = loadings,
              aes(x = PC1 * 1.12, y = PC2 * 1.12, label = metric),
              color = "grey30", size = 3, fontface = "italic") +
    labs(
      x = paste0("PC1 (", round(ve[1], 1), "%)"),
      y = paste0("PC2 (", round(ve[2], 1), "%)"),
      title = title_text
    ) +
    coord_fixed()

  ggsave(filename, p, width = 8, height = 7)
  cat("Saved", filename, "\n")
}

make_pca_biplot(ert_dbh, "PCA of ERT Metrics — DBH (n = 12)", "fig_pca_biplot.pdf")

# 5b: PCA biplot — all heights pooled
make_pca_biplot(ert, "PCA of ERT Metrics — All Heights Pooled", "fig_pca_biplot_all.pdf")

##############################################################################
# FIGURE 6: Best-predictor regression
##############################################################################
cat("\n── Figure 6: Best predictor ──\n")

best_cors <- cor_results %>%
  filter(group == "DBH") %>%
  arrange(pearson_p) %>%
  slice(1)

best_metric <- best_cors$metric
cat("Best predictor (Pearson):", best_metric,
    "r =", round(best_cors$pearson_r, 3),
    "p =", round(best_cors$pearson_p, 4), "\n")

# Regression
lm_fit <- lm(as.formula(paste(best_metric, "~ moisture")), data = ert_dbh)
lm_sum <- summary(lm_fit)
r2 <- round(lm_sum$r.squared, 3)
pval <- lm_sum$coefficients[2, 4]
b0 <- round(coef(lm_fit)[1], 2)
b1 <- round(coef(lm_fit)[2], 2)

fmt_p_best <- if (pval < 0.001) "p < 0.001" else paste0("p = ", formatC(pval, format = "f", digits = 4))

# Spearman for annotation
ct_s_best <- cor.test(ert_dbh$moisture, ert_dbh[[best_metric]],
                      method = "spearman", exact = FALSE)
rho_best <- round(ct_s_best$estimate, 3)
p_spear_best <- ct_s_best$p.value
fmt_ps_best <- if (p_spear_best < 0.001) "p < 0.001" else paste0("p = ", formatC(p_spear_best, format = "f", digits = 4))

# Annotation (plain text with Unicode, rendered via cairo_pdf)
ann_label <- paste0(
  best_metric, " = ", b0, " + ", b1, " \u00d7 Moisture\n",
  "R\u00b2 = ", r2, ", ", fmt_p_best, "\n",
  "r = ", round(best_cors$pearson_r, 3),
  ",  \u03c1 = ", rho_best, " (", fmt_ps_best, ")"
)

fig6 <- ggplot(ert_dbh, aes(x = moisture, y = .data[[best_metric]])) +
  geom_smooth(method = "lm", se = TRUE, color = pal_accent,
              fill = pal_accent, alpha = 0.15, linewidth = 0.8) +
  geom_point(size = 3, color = pal_accent) +
  geom_text_repel(aes(label = tree_id), size = 3, max.overlaps = 20) +
  annotate("text",
           x = max(ert_dbh$moisture) - diff(range(ert_dbh$moisture)) * 0.02,
           y = max(ert_dbh[[best_metric]]) - diff(range(ert_dbh[[best_metric]])) * 0.02,
           label = ann_label,
           hjust = 1, vjust = 1, size = 3.5, color = "grey25",
           lineheight = 1.1) +
  labs(x = "Core Moisture (%)", y = best_metric,
       title = paste0("Best Single Predictor: ", best_metric))

ggsave("fig_best_predictor.pdf", fig6, width = 6, height = 5, device = cairo_pdf)
cat("Saved fig_best_predictor.pdf\n")

##############################################################################
# FIGURE 7: Rank-order comparison (bump/slope chart)
##############################################################################
cat("\n── Figure 7: Rank comparison ──\n")

# Find top 3 best-correlated metrics at DBH
top3 <- cor_results %>%
  filter(group == "DBH") %>%
  arrange(pearson_p) %>%
  slice(1:3)

rank_data <- ert_dbh %>%
  select(tree_id, moisture, all_of(top3$metric)) %>%
  mutate(rank_moisture = rank(-moisture))

for (m in top3$metric) {
  rank_data[[paste0("rank_", m)]] <- rank(-rank_data[[m]])
}

rank_long <- rank_data %>%
  select(tree_id, starts_with("rank_")) %>%
  pivot_longer(-tree_id, names_to = "variable", values_to = "rank") %>%
  mutate(variable = str_remove(variable, "rank_"),
         variable = factor(variable,
                           levels = c("moisture", top3$metric)))

fig7 <- ggplot(rank_long, aes(x = variable, y = rank, group = tree_id)) +
  geom_line(alpha = 0.5, color = "grey50") +
  geom_point(aes(color = variable), size = 3) +
  geom_text_repel(aes(label = tree_id), size = 2.5, max.overlaps = 30,
                  direction = "y") +
  scale_y_reverse(breaks = 1:12) +
  scale_color_manual(values = c("moisture" = "#E41A1C",
                                setNames(c("#377EB8", "#4DAF4A", "#984EA3"),
                                         top3$metric)),
                     guide = "none") +
  labs(x = "", y = "Rank (1 = highest value)",
       title = "Rank Comparison: Core Moisture vs. Top ERT Metrics (DBH)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("fig_rank_comparison.pdf", fig7, width = 8, height = 6)
cat("Saved fig_rank_comparison.pdf\n")

##############################################################################
# FIGURE 8: Within-tree vertical variation
##############################################################################
cat("\n── Figure 8: Vertical variation ──\n")

# Trees with multi-height data
multi_trees <- ert %>%
  group_by(tree_id) %>%
  filter(n_distinct(height) > 1) %>%
  ungroup()

# Order heights
multi_trees <- multi_trees %>%
  mutate(height = factor(height, levels = c("Lower", "DBH", "Upper")))

# Long format for faceting
multi_long <- multi_trees %>%
  select(tree_id, height, moisture, all_of(metric_names)) %>%
  pivot_longer(all_of(metric_names), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_names))

fig8 <- ggplot(multi_long, aes(x = height, y = value, group = tree_id, color = tree_id)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_color_brewer(palette = "Set2", name = "Tree") +
  labs(x = "Measurement Height", y = "Metric Value",
       title = "Within-Tree Vertical Variation of ERT Metrics") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

# Add moisture reference: horizontal line at each tree's core moisture
# Since metrics have different scales, add a small marker instead—
# annotate the DBH point with a subtle indicator
# Actually, add moisture as a secondary annotation: a dashed line per tree
# in a "Moisture" panel
moisture_panel <- multi_trees %>%
  select(tree_id, height, moisture) %>%
  distinct() %>%
  mutate(metric = "Core Moisture\n(reference)", value = moisture,
         metric = factor(metric))

# Combine
multi_long2 <- multi_long %>%
  mutate(metric = as.character(metric)) %>%
  bind_rows(moisture_panel %>% mutate(metric = as.character(metric))) %>%
  mutate(metric = factor(metric, levels = c(metric_names, "Core Moisture\n(reference)")))

fig8 <- ggplot(multi_long2, aes(x = height, y = value, group = tree_id, color = tree_id)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set2", name = "Tree") +
  labs(x = "Measurement Height", y = "Value",
       title = "Within-Tree Vertical Variation of ERT Metrics") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("fig_vertical_variation.pdf", fig8, width = 12, height = 10)
cat("Saved fig_vertical_variation.pdf\n")

##############################################################################
# FIGURE 9: ERT Metric Redundancy
##############################################################################
cat("\n── Figure 9: Metric redundancy ──\n")

cor_metrics <- cor(ert_dbh[, metric_names])

pdf("fig_metric_redundancy.pdf", width = 8, height = 7)
pheatmap(
  cor_metrics,
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "grey90",
  fontsize = 10,
  fontsize_number = 9,
  main = "Pairwise Correlations Among ERT Metrics (DBH Data)"
)
dev.off()
cat("Saved fig_metric_redundancy.pdf\n")

##############################################################################
# TABLE 10: Summary statistics CSV
##############################################################################
cat("\n── Table 10: Summary CSV ──\n")

avg_cols <- ert_avg %>%
  select(tree_id, all_of(metric_names)) %>%
  rename_with(~ paste0(.x, "_avg"), all_of(metric_names))

summary_table <- ert_dbh %>%
  select(tree_id, moisture, all_of(metric_names)) %>%
  left_join(avg_cols, by = "tree_id") %>%
  arrange(tree_id)

write_csv(summary_table, "validation_summary.csv")
cat("Saved validation_summary.csv\n")

cat("\n══ Analysis complete ══\n")
