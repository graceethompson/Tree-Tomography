# Revision analysis (cjfr-2026-0202)

Sensitivity and robustness analyses supporting the major revision of
"Internal decay in living trees: a quantitative tomography framework and its
application in a temperate forest." All scripts read only the existing repo data
(`data/`, `images/`, plus `analysis/revision/data/`) and write to
`analysis/revision/output/`; nothing else in the repository is touched.

Each analysis script independently reconstructs the manuscript's ERT PC1
(species-normalized z-scores, PCA fit on the 57 main trees, oriented high = wet)
via `scripts/revision_common.R` and reproduces the published values (PC1 47.1% /
PC2 30.0% variance; classification 30/12/12/3), consistent with
`code/final_phase_and_scans.R`.

## How to run

R (>= 4.x) with dplyr, ggplot2, patchwork, betareg, emmeans, DHARMa, mclust,
magick, jpeg. Run each script from the repository root, e.g.
`Rscript analysis/revision/scripts/reanalysis.R`.

(These analyses were originally prototyped in Python; the R scripts replace the
prototypes 1:1 and were verified to reproduce the same statistics.)

## Scripts

Shared helper: `revision_common.R` (data loading, PC1 reconstruction,
classification, Otsu/GMM helpers, tomogram lookup).

### Threshold / classification sensitivity

| Script | Purpose |
|---|---|
| `reanalysis.R` | Reproduce reviewer-cited statistics: ANOVAs (one-way vs two-way), hurdle-model LRTs, classification counts, SoT damage distribution, DBH/site confound check, validation correlations, carbon denominators |
| `statistical_tests.R` | Fisher exact tests on the incipient site contrast (per threshold; generalists-only subset); calibration-anchor uncertainty (CI on the calibration line vs single-tree prediction interval) |
| `threshold.R`, `threshold_fig.R` | Unimodality of the ERT axis (GMM/BIC); threshold-rule menu (mean, median, Otsu, GMM, ±SD); continuous threshold sweep of the site contrast |
| `absolute.R` | Same analyses in absolute units (Ω·m; predicted moisture via the hemlock calibration); species-baseline confound of the absolute axis; literature/validation anchor comparison |
| `boundaries.R` | Boundary-by-boundary stability: flips at I↔II, II↔III, III↔IV, I↔IV under threshold perturbation; SoT structural threshold 1–10% |
| `compare.R` | Classification under alternative moisture axes (PC1, CMA, absolute resistivity); PC1–CMA orthogonality |
| `capstone.R` | Robustness grid: site-contrast direction across ~31 axis × threshold combinations; SoT-threshold stability of structural prevalence; species severity counts |
| `incipient.R` | Incipient-onset sensitivity on the species-median-referenced axis; comparison with structurally confirmed decay trees; 6-cell dead-band variant |
| `rules.R`, `flow.R` | Consolidated classification rules; decision-tree diagram |
| `cavity.R` | Active-vs-cavity (III/IV) inspection among structurally damaged trees; core-wetness criterion |
| `montage.R`, `contact.R`, `grouped.R`, `paired_montages.R` | Tomogram visualisations: per-category exemplars, contact sheets, and SoT+ERT paired montages grouped by scheme |
| `make_assignments_table.R` | Writes `output/scheme_assignments.csv`: every tree's category under each scheme |

### Reviewer-requested additions

| Script | Reviewer point | Purpose |
|---|---|---|
| `hurdle_figs.R` | #6, #7 | Occurrence-stage fitted probabilities over the jittered binary observations by species × site; severity-stage beta-regression fit over the individual decayed trees (n per species shown, incl. the single *N. sylvatica* reference); severity residual diagnostics (the single-reference tree has leverage 1 → undefined residual); DHARMa panel for the binomial stage |
| `sensor_counts.R` | #8 | Sensor counts + DBH for all 57 trees (from `data/sot_image_annotations.csv`, read off every tomogram); sensor spacing as effective-resolution proxy; covariance of resolution with site, species, and the damage estimate |
| `vertical_variation.R` | #10 | Multi-height validation scans (Lower/DBH/Upper): within-tree vertical variation in SoT damage and ERT resistivity; how often a BH-only plane misses decay flagged at another height |
| `green_zone.R` | repro #2 | How green (intermediate-velocity) zones were tabulated. Resolves the ambiguity: `percent_solid_wood` = PiCUS "Solid wood %" (brown only; header-verified 56/57), and the recorded `percent_damaged` EXCLUDED green (violet+blue only) — the manuscript's "non-brown colors combined" wording needs correcting. Includes an independent pixel-color re-measurement (r ≈ 0.98–0.99) and the sensitivity of classifications and the site contrast to counting green as damaged |

`data/sot_image_annotations.csv`: sensor count and PiCUS "Solid wood %" header
value read from every SoT tomogram (main survey + validation hemlocks).

## Outputs

`output/` — all generated figures (threshold sweep, absolute-unit panels,
boundary stability, binning comparison, 6-cell scheme, decision tree, scan
montages, hurdle-model panels, sensor-count panels, vertical variation,
green-zone panels) and tables (`CJFR-absolute-thresholds.csv`,
`scheme_assignments.csv`, `CJFR-sensor-counts.csv`).

## Headline numbers

- ERT anomaly axis is unimodal (2-component GMM fits worse by BIC); any threshold is
  an operational convention. SoT damage is effectively bimodal (0% or ≥6%), so the
  1% structural cut is stable (identical classifications 1–5% except 2 trees).
- The wetland>upland incipient direction holds across every resistivity-based
  axis × threshold combination and both normalizations; the contrast is marginally
  significant (Fisher exact p ≈ 0.10; generalists-only 32% vs 11%, p = 0.23) and is
  presented as exploratory (one stand per habitat type).
- Scheme agreement (identical 4-category call): PC1 vs species-median 79%,
  species-median vs absolute 82%, PC1 vs absolute 68%; all agree on all
  structural (sound vs damaged) calls.
- Absolute anchor from the hemlock calibration: moisture>100% ⇔ 402 Ω·m
  (95% CI on the stand-level anchor 275–529; single-tree PI ≈ −50 to 850).
- Sensor counts: 7 sensors for 50 trees; 8 for the seven largest *Q. rubra* at
  EMS — sensor count and spacing DO covary with site (spacing coarser at EMS,
  p = 0.013) and weakly with the damage estimate (ρ = 0.28, p = 0.04); with
  spacing entered first, the one-way site damage contrast attenuates to p = 0.08.
  Supports reframing the site contrast as exploratory.
- Vertical variation (6 multi-height hemlocks): within-tree SoT range averages
  3.2 points (max 8); the BH plane understates the section maximum by 1.7 points
  on average and misses 1 of 5 trees flagged >1% at some height. ERT resistivity
  rises with height in 5/5 trees (wetter near the base), consistent with basal
  decay columns tapering upward.
- Green zones: recorded damage excluded green (median green share 3%, up to 26%).
  Counting green as damage would flip the >1% structural call for 21 of 57 trees;
  the site percent-damage contrast is not sensitive to the choice (one-way ANOVA
  p = 0.027 recorded vs p = 0.006 with green).
- The published Fig. 4 severity panel (predict type="response") is correct, but
  the exploratory `code/decay_hurdle_analysis.R` double-transforms betareg
  emmeans with plogis(); its Panel C values are wrong and should not be reused.
