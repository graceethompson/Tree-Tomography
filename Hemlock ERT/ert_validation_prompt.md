# ERT Validation Analysis — Claude Code Prompt

You are helping me build a validation analysis for an Electrical Resistivity Tomography (ERT) paper about tree moisture. I have two datasets to merge and analyze:

## Data

1. **ERT image-derived metrics** — `ERT_results_2026-03-05.csv` in the current directory. Columns: `Filename, Mean, Median, SD, CV, Gini, Entropy, CMA, RadialGradient, NPixels`. Each row is one processed ERT image. Filenames follow the pattern `HF_XXXXXX_Height.jpg` where Height is DBH, Lower, or Upper.

2. **Core moisture data** — `MC_Tomo_paper_Jon.xlsx` in the current directory, Sheet1. Columns: `ID` (e.g. HF_201290) and `Moisture content (%)`. These are directly-measured moisture values from increment cores taken at approximately DBH height on 12 hemlock trees.

## Data Cleaning

- Parse `tree_id` and `height` from the Filename column (split on last underscore, strip `.jpg`).
- There are duplicate rows for some tree×height combinations (from re-processing). For each unique tree_id × height combination, **keep only the last row** (most recent processing run).
- All 12 trees have DBH images. Only 6 trees (HF_021110, HF_031205, HF_061273, HF_070836, HF_101924, HF_978) have Lower and/or Upper images.

## Analyses to Produce

Generate all figures as publication-quality PDFs saved to the current directory. Use a clean, minimal style (no heavy gridlines, reasonable font sizes). Use a consistent color palette throughout. All figures should be clearly labeled with axis titles and, where relevant, tree IDs as point labels.

The ERT metrics to analyze are: Mean, Median, SD, CV, Gini, Entropy, CMA, RadialGradient (8 metrics total). The validation target is `Moisture content (%)` from the core data.

### 1. Scatter Plots: Each ERT Metric vs. Core Moisture (at DBH)
- One multi-panel figure (e.g., 2×4 grid) with 8 scatter plots.
- x-axis = core moisture (%), y-axis = ERT metric value.
- Label each point with the tree ID.
- Include Pearson r, Spearman ρ, and p-values as text annotations on each panel.
- Overlay a linear regression line with 95% CI shading where the correlation is significant (p < 0.05).
- Save as `fig_scatter_dbh.pdf`.

### 2. Scatter Plots: Each ERT Metric vs. Core Moisture (by Height)
- Same as above but with **three sets**: DBH-only, Lower-only, Upper-only.
- Since Lower and Upper only have 6 trees, use a different layout or combined figure.
- Could be a 3-row × 8-column panel figure, or three separate files — use your judgment for readability.
- Color/shape-code points by height.
- Save as `fig_scatter_by_height.pdf`.

### 3. Scatter Plots: Average-Across-Heights ERT Metrics vs. Core Moisture
- For each tree, average its ERT metrics across all available heights (some trees have 1 height, some have 3).
- Same 2×4 scatter panel as #1, using the averaged metrics.
- Save as `fig_scatter_avg.pdf`.

### 4. Correlation Heatmap
- A heatmap showing Pearson and Spearman correlations between each ERT metric and core moisture.
- Rows = ERT metrics, Columns = correlation type × height grouping (DBH, Lower, Upper, Average).
- Annotate cells with r or ρ values. Star significant correlations (p < 0.05).
- Save as `fig_correlation_heatmap.pdf`.

### 5. PCA of ERT Metrics, Colored by Core Moisture
- Standardize the 8 ERT metrics (z-score).
- Run PCA on the DBH-only data (n=12 trees).
- Biplot: PC1 vs PC2, points colored by a continuous color scale representing core moisture (%).
- Show loading vectors for each metric.
- Label points with tree IDs.
- Include variance explained on axes.
- Save as `fig_pca_biplot.pdf`.
- Also run PCA on all heights pooled (n=28 unique tree×height combos after dedup) and color by core moisture. Save as `fig_pca_biplot_all.pdf`.

### 6. Best-Predictor Regression
- Identify which single ERT metric has the strongest correlation with core moisture at DBH.
- Produce a clean single-panel scatter plot with regression line, equation, R², and p-value.
- Save as `fig_best_predictor.pdf`.

### 7. Rank-Order Comparison
- For DBH data: rank the 12 trees by core moisture and by each ERT metric.
- Produce a bump chart or slope graph showing how ranks compare between core moisture and the top 2-3 best-correlated ERT metrics.
- Save as `fig_rank_comparison.pdf`.

### 8. Within-Tree Vertical Variation
- For the 6 trees with multi-height data, plot each ERT metric across Lower → DBH → Upper.
- One line per tree, faceted by metric.
- Add a horizontal reference line or marker showing the core moisture value for context.
- Save as `fig_vertical_variation.pdf`.

### 9. ERT Metric Redundancy (Pairwise Correlations)
- Pairwise scatter matrix or correlation heatmap among the 8 ERT metrics (using DBH data).
- This shows which metrics are independently informative vs. redundant.
- Save as `fig_metric_redundancy.pdf`.

### 10. Summary Statistics Table
- Save a CSV (`validation_summary.csv`) with one row per tree containing: tree ID, core moisture, and all ERT metrics at DBH, plus the averaged-across-heights metrics.

## Technical Notes

- **Use R for everything.** Write one or more `.R` scripts to produce all outputs.
- Key packages: `tidyverse` (ggplot2, dplyr, tidyr, stringr, readr), `readxl`, `ggrepel` (point labels), `corrplot` or `pheatmap` (heatmaps), `ggfortify` or `factoextra` (PCA biplots), `patchwork` or `cowplot` (multi-panel layouts).
- Install any missing packages with `install.packages()` as needed.
- For correlations, use `cor.test()` with `method = "pearson"` and `method = "spearman"`.
- For PCA, use `prcomp()` with `scale. = TRUE`.
- Save all PDFs using `ggsave()` or `pdf()` / `dev.off()`. Use reasonable dimensions (e.g., 12×8 inches for multi-panel figures, 6×5 for single panels).
- Use `theme_minimal()` or `theme_classic()` as a base ggplot theme for a clean, publication-ready look.
- Use `ggrepel::geom_text_repel()` to label points with tree IDs without overlap.
