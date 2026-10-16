# Tree-Tomography

Code, working data and figures for Thompson et al., "Internal decay in living trees: a quantitative tomography framework and its application in a temperate forest" (preprint: https://doi.org/10.64898/2026.06.10.730433). The repository also holds the ERT image analysis app and the build script for the tomography data package deposited in the Harvard Forest Data Archive.

## Study design

### SoT and ERT measurements

Corresponding to an existing biometric plot network established around the footprint of an eddy covariance tower at Harvard Forest, we defined our upland site as EMS. Our wetland site was the Black Gum Swamp (BGS), a forested peatland within the same Prospect Hill tract. We selected four tree species common in either the EMS or BGS site, or both: the generalists red maple (_Acer rubrum_) and eastern hemlock (_Tsuga canadensis_), the upland specialist northern red oak (_Quercus rubra_) and the wetland specialist black gum (_Nyssa sylvatica_). Each tree was representative of the average size class (dbh) for its species, estimated from Harvard Forest's ForestGEO data. Sonic (SoT) and electrical resistivity (ERT) tomography was performed once at each tree. The analysed set is 57 trees: 10 wetland red maples, 10 upland red maples, 9 red oaks, 10 black gums, 9 wetland eastern hemlocks and 9 upland eastern hemlocks.

### Eastern hemlock water content

Twelve eastern hemlocks were scanned with ERT and SoT at breast height, six of them also at a lower and an upper height, to validate ERT-derived moisture metrics against core-extracted gravimetric moisture content.

## Repository structure

```
Tree-Tomography/
├── code/             main analyses and the main-text figures
├── analysis/         sensitivity and robustness analyses (supplement)
│   ├── scripts/      R scripts; common.R holds the shared helpers
│   ├── data/         sot_image_annotations.csv (sensor counts, PiCUS headers)
│   └── output/       supplementary figures and tables
├── app/ERT_App/      ERT image analysis app (Shiny) and batch processor
├── data/             tabular data read by the scripts
├── images/           SoT and ERT tomogram exports
├── output/           figures and tables written by code/
└── data_package/     build script and metadata for the archived data package
```

### `data/`

- `Tree_ID_info.csv`: the 57 study trees (plot, species, tree tag, dbh, SoT percent damaged and percent solid wood)
- `ERT_application_results.csv`: 8 ERT metrics per tree (mean, median, sd, cv, gini, entropy, cma, radial gradient), from the app
- `hemlock/`: the cored hemlocks
  - `validation_summary.csv`: ERT metrics and core moisture per tree (written by `code/ert_validation_analysis.R`)
  - `validation_phases.csv`: phase-diagram class of each cored hemlock (written by `code/final_phase_and_scans.R`)
  - `SOT_results.csv`: SoT percent damaged per scan (tree × height)
  - `ERT_results_2026-03-05.csv`: ERT metrics per scan (tree × height)
  - `hemlock_core_moisture.xlsx`: gravimetric moisture content of the cores
- `supplementary/BGS_VRP_2025.csv`: 2025 variable-radius (prism) survey of the Black Gum Swamp stand
- `external/` (not tracked): the ForestGEO census needed by `analysis/scripts/stand_context.R`. Download `hf253-06-stems-2019.csv` from Harvard Forest dataset HF253 (Orwig, Foster & Ellison 2023, https://doi.org/10.6073/pasta/818789a882a318c1d7f3fc43a2289e12) into this folder.

### `images/`

- `main_SoT/`, `main_ERT/`: tomograms of the study trees
- `main_ERT_normalized/`: ERT tomograms on the normalized colour scale
- `hemlock_SoT/`, `hemlock_ERT/`: tomograms of the cored hemlocks (`<tree>_<Lower|DBH|Upper>.jpg`)
- `app_screenshot.png`: screenshot used in the supplementary app figure

### `app/ERT_App/`

An R/Shiny application that extracts quantitative metrics from PiCUS ERT tomogram images. It maps pixel colours to calibrated resistivity values from the colour bar and computes 8 summary statistics (mean, median, SD, CV, Gini, entropy, CMA, radial gradient) and an 8-ring radial profile. It detects the PiCUS measurement polygon automatically or takes a manual outline, and exports results as CSV.

- `app.R`: the Shiny app. Live version: https://jgewirtzman-tree-tomography.share.connect.posit.cloud. Run locally with `shiny::runApp("app/ERT_App")`.
- `processing.R`: the image-processing functions shared by the app and the batch processor
- `batch.R`: headless batch processor, e.g. `Rscript app/ERT_App/batch.R --input images/main_ERT --output results` (`--help` lists the options)
- `manifest.json`: deployment manifest for Posit Connect Cloud

## Requirements

R 4.4 (tested with R 4.4.3). Install the packages once:

```r
install.packages(c(
  "tidyverse", "patchwork", "cowplot", "gridExtra", "scales", "ggrepel",
  "RColorBrewer", "factoextra", "pheatmap", "betareg", "emmeans", "broom",
  "DHARMa", "lme4", "lmPerm", "mclust", "cluster", "readxl", "magick",
  "jpeg", "png", "sp",
  # ERT app
  "shiny", "shinydashboard", "DT", "imager", "DescTools"))
# data package build only
install.packages("remotes"); remotes::install_github("EDIorg/EMLassemblyline")
```

## Running the analyses

Run every script from the repository root (open `Tree-Tomography.Rproj`, or `cd` into the repository), for example `Rscript code/combined_figures.R`. Scripts create the output folders they need. Each step takes seconds.

1. Main-text figures and the phase diagram

   ```sh
   Rscript code/fig_ert_workflow.R
   Rscript code/combined_figures.R
   Rscript code/final_phase_and_scans.R
   Rscript code/ert_validation_analysis.R
   ```

2. Descriptive statistics (ANOVAs, boxplots). `treedecay_app.R` uses the table that `treedecay_og.R` builds, so run them in one R session. `treedecay_og.R` writes all its figures and then stops in a final exploratory beta-regression block, hence the `try()`:

   ```sh
   Rscript -e 'try(source("code/treedecay_og.R")); source("code/treedecay_app.R")'
   ```

3. Other analyses in `code/`

   ```sh
   Rscript code/phase_diagram_analysis.R
   Rscript code/decay_presence_severity.R
   Rscript code/crop_ert_image.R
   ```

4. Supplementary analyses. Run `make_assignments_table.R` before `absolute.R`, and `sensor_counts.R` before `make_tree_table.R`. `stand_context.R` needs the HF253 census in `data/external/`.

   ```sh
   for s in reanalysis statistical_tests threshold threshold_fig make_assignments_table \
            absolute boundaries compare capstone incipient rules cavity make_fisher_table \
            hurdle_figs sensor_counts make_tree_table vertical_variation sot_scale_check \
            green_zone stand_context montage contact grouped paired_montages plates; do
     Rscript analysis/scripts/$s.R
   done
   ```

Some scripts also print the statistics reported in the text (for example `reanalysis.R`, `statistical_tests.R`, `threshold.R`).

### Exploratory scripts

These scripts are kept for completeness. No figure or table depends on them, and they do not run end to end with the current data:

- `code/treedecay_dbh.R`: decay against stem diameter (stops at an undefined object, `BGS_treedecayinfo`)
- `code/check_agreement.R`: agreement between SoT and ERT classifications
- `code/phase_image_panels.R`: tomogram panels arranged by phase-diagram position, written to `output/panels/` (not tracked), and `output/figures/phase_diagrams_by_species.pdf`. The hemlock panels also need the hemlock SoT percent-calculation exports, which are not in the repository (set `HEMLOCK_TOMO_DIR` to their folder).

`check_agreement.R` and `phase_image_panels.R` expect a `Mean` column in `data/hemlock/validation_summary.csv`, which now names it `Conductance`.

## Script → output map

Paths are relative to the repository root.

| Script | Outputs | Used for |
|---|---|---|
| `code/fig_ert_workflow.R` | `output/figures/fig_ert_workflow.{pdf,png}`, `fig_si_app_screenshot.{pdf,png}` | Fig. 1 (ERT workflow); supplementary app figure |
| `code/combined_figures.R` | `output/figures/fig2_sot_damaged_combined.*`, `fig3_ert_mean_cv_combined.*`, `fig4_hurdle_3panel.*` (and `fig4_hurdle_4panel.*`) | Figs. 2–4 (SoT damage, ERT mean and CV, hurdle model) |
| `code/final_phase_and_scans.R` | `output/figures/final_phase_diagram.*`, `pca_biplot.*`, `quadrant_distribution.*`; `output/hemlock_figures/fig_best_predictor_pc1.*`, `fig_scans_by_height_with_pc1.pdf`, `fig_scans_by_height_pc1_only.pdf` (large, not tracked); `data/hemlock/validation_phases.csv` | Phase diagram (Fig. 7), PCA biplot, quadrant distributions, hemlock validation scans |
| `code/ert_validation_analysis.R` | `output/hemlock_figures/fig_scatter_*.pdf`, `fig_correlation_heatmap.pdf`, `fig_pca_biplot*.pdf`, `fig_best_predictor.pdf`, `fig_rank_comparison.pdf`, `fig_vertical_variation.pdf`, `fig_metric_redundancy.pdf`; `data/hemlock/validation_summary.csv` | Hemlock ERT validation against core moisture |
| `code/treedecay_og.R` | `output/figures/sot_*.pdf` | SoT descriptive statistics and ANOVAs |
| `code/treedecay_app.R` | `output/figures/ert_*_by_species_site.pdf` | ERT descriptive statistics and ANOVAs |
| `code/phase_diagram_analysis.R` | `output/figures/phase_diagram_candidates.pdf`, `phase_diagram_kmeans.pdf`; `output/tables/phase_diagram_results.csv` | Search over 714 ERT axis × threshold combinations |
| `code/decay_presence_severity.R` | `output/figures/decay_presence_severity.pdf` | Earlier presence/severity figure |
| `code/crop_ert_image.R` | `output/cropped/309_18vii24_cropped.png` (default; takes input and output paths as arguments) | Example of cropping a tomogram to its measurement polygon |
| `analysis/scripts/threshold_fig.R` | `analysis/output/threshold-fig.png` | Threshold sensitivity of the ERT axis (SI) |
| `analysis/scripts/absolute.R` | `absolute-fig.png`, `absolute-thresholds.csv`, `absolute-anchor-classification.csv` | Classification in absolute units (SI) |
| `analysis/scripts/make_assignments_table.R` | `scheme_assignments.csv` | Each tree's category under each scheme (SI table) |
| `analysis/scripts/make_fisher_table.R` | `fisher-tests.csv` | Fisher exact tests of the site contrast (SI table) |
| `analysis/scripts/compare.R` | `binning-compare.png` | Classification under alternative ERT axes (SI) |
| `analysis/scripts/hurdle_figs.R` | `hurdle-fig.png` | Severity-stage residual diagnostics (SI) |
| `analysis/scripts/sensor_counts.R` | `sensor-counts.csv`, `sensor-counts.png` | Sensor counts and spacing (SI table and figure) |
| `analysis/scripts/vertical_variation.R` | `vertical-variation.png` | Vertical variation in the multi-height hemlocks (SI) |
| `analysis/scripts/boundaries.R` | `boundary-fig.png` | Stability of the structural threshold (SI) |
| `analysis/scripts/incipient.R` | `6cell.png` | Graded variant of the classification (SI) |
| `analysis/scripts/paired_montages.R` | `paired-PC1.png` | Paired SoT and ERT tomograms of all trees by category (SI) |
| `analysis/scripts/plates.R` | `plate-a.png` … `plate-d.png` | Tomogram plates, one per category (SI) |
| `analysis/scripts/stand_context.R` | `stand-context.csv` | Stand context of the sampled trees (SI table) |
| `analysis/scripts/make_tree_table.R` | `tree-table.csv` | Per-tree table: sensors, dbh, categories (SI table) |
| `analysis/scripts/reanalysis.R`, `statistical_tests.R`, `threshold.R`, `capstone.R`, `rules.R`, `sot_scale_check.R` | printed statistics (`threshold.R` also writes `additional/threshold-rules.png`) | Statistics reported in the text and supplement |
| `analysis/scripts/green_zone.R`, `cavity.R`, `flow.R`, `montage.R`, `contact.R`, `grouped.R` | `analysis/output/additional/*.png` (not tracked) | Supporting checks and visualisations, not in the paper |

`analysis/output/` paths in the table are relative to that folder. `analysis/scripts/common.R` reconstructs the ERT PC1 used throughout (species-normalized z-scores, PCA on the 57 study trees, oriented so that higher = wetter) and reproduces the values of `code/final_phase_and_scans.R` (PC1 47.1% and PC2 30.0% of variance; classification 30/12/12/3).

Files in `output/` that no current script writes: `output/hemlock_figures/fig_scans_by_height.pdf`, `fig_scans_by_moisture.pdf`, `quadrant_distribution.{pdf,png}`, and the 600-dpi exports of Figs. 1–4 in `output/figures/hires/`.

## Building the data package

```sh
Rscript data_package/build_data_package.R
```

The script copies the tables from `data/` and `analysis/data/`, converts `data/hemlock/hemlock_core_moisture.xlsx` to CSV, zips the tomogram images, and writes EML metadata from the text in `data_package/metadata/`. The package goes to `data_package/output/package/` (not tracked).

## Archived data

The tomography data and images are archived in the Harvard Forest Data Archive (dataset ID to be added). The copies in `data/` and `images/` are a working copy for running the code; cite the archived dataset for the data.

## Citation and license

Code: MIT (`LICENSE`). Please cite the archived release of this repository (Thompson and Gewirtzman 2026, Zenodo, https://doi.org/10.5281/zenodo.23421500 for the latest version; v1.0.0: https://doi.org/10.5281/zenodo.23421501) and Thompson et al. (2026), https://doi.org/10.64898/2026.06.10.730433.
