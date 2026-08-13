# Revision analysis (cjfr-2026-0202)

Sensitivity and robustness analyses supporting the major revision of
"Internal decay in living trees: a quantitative tomography framework and its
application in a temperate forest." All scripts read only the existing repo data
(`data/`, `images/`) and write to `analysis/revision/output/` and
`analysis/revision/tables/`; nothing else in the repository is touched.

Each script independently reconstructs the manuscript's ERT PC1 (species-normalized
z-scores, PCA fit on the 57 main trees, oriented high = wet) and reproduces the
published values (PC1 47.1% / PC2 30.0% variance; classification 30/12/12/3),
consistent with `code/final_phase_and_scans.R`.

## How to run

Python 3 with pandas, numpy, scipy, scikit-learn, statsmodels, matplotlib, pillow.
Run from the repository root, e.g. `python3 analysis/revision/scripts/reanalysis.py`.

## Scripts

| Script | Purpose |
|---|---|
| `reanalysis.py` | Reproduce reviewer-cited statistics: ANOVAs (one-way vs two-way), hurdle-model LRTs, classification counts, SoT damage distribution, DBH/site confound check, validation correlations, carbon denominators |
| `statistical_tests.py` | Fisher exact tests on the incipient site contrast (per threshold; generalists-only subset); calibration-anchor uncertainty (CI on the calibration line vs single-tree prediction interval) |
| `threshold.py`, `threshold_fig.py` | Unimodality of the ERT axis (GMM/BIC); threshold-rule menu (mean, median, Otsu, GMM, ±SD); continuous threshold sweep of the site contrast |
| `absolute.py` | Same analyses in absolute units (Ω·m; predicted moisture via the hemlock calibration); species-baseline confound of the absolute axis; literature/validation anchor comparison |
| `boundaries.py` | Boundary-by-boundary stability: flips at I↔II, II↔III, III↔IV, I↔IV under threshold perturbation; SoT structural threshold 1–10% |
| `compare.py` | Classification under alternative moisture axes (PC1, CMA, absolute resistivity); PC1–CMA orthogonality |
| `capstone.py` | Robustness grid: site-contrast direction across ~31 axis × threshold combinations; SoT-threshold stability of structural prevalence; species severity counts |
| `incipient.py` | Incipient-onset sensitivity on the species-median-referenced axis; comparison with structurally confirmed decay trees; 6-cell dead-band variant |
| `rules.py`, `flow.py` | Consolidated classification rules; decision-tree diagram |
| `cavity.py` | Active-vs-cavity (III/IV) inspection among structurally damaged trees; core-wetness criterion |
| `montage.py`, `contact.py`, `grouped.py`, `paired_montages.py` | Tomogram visualisations: per-category exemplars, contact sheets, and SoT+ERT paired montages grouped by scheme |
| `make_assignments_table.py` | Writes `tables/scheme_assignments.csv`: every tree's category under each scheme |

## Outputs

`output/` — generated figures (threshold sweep, absolute-unit panels, boundary
stability, binning comparison, 6-cell scheme, decision tree, scan montages).
`tables/` — `CJFR-absolute-thresholds.csv` (absolute threshold menu),
`scheme_assignments.csv` (per-tree assignments under each scheme).

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
