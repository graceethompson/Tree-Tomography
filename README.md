# Tree-Tomography

Natural history study of internal tree decay within and amongst dominant tree species present along an environmental gradient at the Harvard Forest.

## Study design

### SoT and ERT measurements

Corresponding to an existing biometric plot network established around the footprint of an eddy covariance tower at Harvard Forest, we defined our upland site as EMS. Our wetland site was the Black Gum Swamp, which is a forested peatland within the same Prospect Hill tract as the EMS site. We selected four tree species common in either the EMS or BGS site, or both. The generalist species were red maple (_Acer rubrum_) and eastern hemlock (_Tsuga canadensis_). The upland specialist species was northern red oak (_Quercus rubra_), and the wetland specialist species was black gum (_Nyssa sylvatica_). We selected 10 trees per species for observation. Each tree had to be representative of the average size class (dbh) for its respective species estimated using Harvard Forest's ForestGEO data. Sonic and electric resistivity tomography was performed once at each tree. Some electric resistivity tomograms were not fully uploaded at the end of the study, and the tree count at each site was as follows: 10 wetland red maples, 10 upland red maples, 9 red oaks, 10 black gums, 9 wetland eastern hemlocks, and 9 upland eastern hemlocks.

### Eastern hemlock water content

Twelve eastern hemlocks were scanned with ERT and SoT at breast height, six of them also at a lower and an upper height, to validate ERT-derived moisture metrics against core-extracted gravimetric moisture content.

## ERT Image Analysis App

An open-source R/Shiny application for extracting quantitative metrics from PiCUS ERT tomogram images. The app maps pixel colors to calibrated resistivity values via colorbar extraction and computes 8 summary statistics (mean, median, SD, CV, Gini, entropy, CMA, radial gradient) along with an 8-ring radial profile.

**Live app:** https://jgewirtzman-tree-tomography.share.connect.posit.cloud

**Features:**
- Automatic boundary detection of PiCUS blue measurement polygon
- Manual polygon delineation as alternative
- Log-scale colorbar calibration with user-defined control points
- Batch processing with CSV export

Source code: `app/ERT_App/app.R`

## Repository structure

### `app/`

- `ERT_App/app.R` — Shiny app for ERT image analysis (deployed to Posit Connect Cloud)
- `ERT_App/manifest.json` — Deployment manifest for Posit Connect Cloud

### `data/`

- `Tree_ID_info.csv` — 57 main study trees: plot, species, tree tag, dbh, SoT percent damaged/solid wood
- `ERT_application_results.csv` — 8 ERT metrics (mean, median, sd, cv, gini, entropy, cma, radial gradient) per tree
- `hemlock/` — Hemlock validation data
  - `validation_summary.csv` — 12 hemlock trees with ERT conductance metrics and core moisture
  - `SOT_results.csv` — SoT percent damaged per hemlock section
  - `ERT_results_2026-03-05.csv` — Per-height ERT metrics for hemlocks
  - `MC_Tomo_paper_Jon.xlsx` — Core moisture content data
- `supplementary/` — Inputs of the original data compilation (`code/Tree_Tomography_Data_Analysis.Rmd`) and the 2025 Black Gum Swamp prism survey (`BGS_VRP_2025.csv`)

### `code/`

- `final_phase_and_scans.R` — Main analysis: 2x2 decay phase diagram (species-normalized PCA on ERT metrics vs SoT structural loss), PCA biplot, quadrant distributions, hemlock scan panels
- `phase_diagram_analysis.R` — Systematic search across 714 ERT axis / threshold combinations
- `phase_image_panels.R` — Tomogram image panels arranged by phase diagram position
- `treedecay_og.R` — SoT descriptive statistics and species/site comparisons (boxplots, ANOVAs)
- `treedecay_app.R` — ERT descriptive statistics and species/site comparisons (boxplots, ANOVAs)
- `treedecay_relfreq.R` — Relative frequency analysis
- `ert_validation_analysis.R` — Hemlock ERT validation against core moisture
- `check_agreement.R` — Agreement analysis between SoT and ERT
- `decay_hurdle_analysis.R` — Hurdle model analysis of decay presence/severity
- `decay_presence_severity.R` — Decay presence/severity analysis with figures
- `fig_ert_workflow.R` — 3-panel ERT workflow figure (raw tomogram, calibration, quantitative output) and the app screenshot figure
- `combined_figures.R` — combined SoT, ERT and hurdle-model figures (`fig2_*`, `fig3_*`, `fig4_*` in `output/figures/`)
- `treedecay_dbh.R` — decay versus stem diameter
- `crop_ert_image.R` — crops an ERT tomogram to its measurement polygon (example output in `output/cropped/`)
- `Tree_Tomography_Data_Analysis.Rmd` — Data compilation (Max Lutz)

### `images/`

- `main_ERT/` — Main study ERT tomogram images
- `main_SoT/` — Main study SoT tomogram images
- `main_ERT_normalized/` — Normalized ERT images
- `hemlock_ERT/` — Hemlock ERT tomogram images
- `hemlock_SoT/` — Hemlock SoT tomogram images

### `output/`

- `figures/` — Main study figures (phase diagram, PCA biplot, quadrant distributions, SoT/ERT boxplots, ERT workflow figure)
- `hemlock_figures/` — Hemlock validation and scan figures
- `panels/` — Tomogram image panels by phase diagram position (large; not tracked, regenerate with `code/phase_image_panels.R`)
- `tables/` — Analysis result tables

## Running the analysis

All scripts use relative paths from the project root. Set your working directory to the repository root before running. The primary analysis pipeline is:

1. `treedecay_og.R` — SoT descriptive analysis
2. `treedecay_app.R` — ERT descriptive analysis
3. `final_phase_and_scans.R` — Phase diagram and PCA analysis
4. `fig_ert_workflow.R` and `combined_figures.R` — manuscript figures
5. `ert_validation_analysis.R` — hemlock validation figures

The sensitivity and robustness analyses of the supplement are in `analysis/revision/` (see its README). The hemlock tomogram panels of `phase_image_panels.R` also need the hemlock SoT percent-calculation exports, which are not in the repository (`HEMLOCK_TOMO_DIR`). The ForestGEO census used by `analysis/revision/scripts/stand_context.R` and `dbh_reference.R` is Harvard Forest dataset HF253 (`hf253-06-stems-2019.csv`), placed in `data/external/`.

R packages: tidyverse, patchwork, cowplot, gridExtra, scales, ggrepel, factoextra, pheatmap, betareg, emmeans, broom, DHARMa, lme4, lmPerm, mclust, cluster, DescTools, ROCR, TOC, readxl, magick, jpeg, png, imager, sp; for the app also shiny, shinydashboard and DT.

## Data

The tomography data and images in `data/` and `images/` are a working copy for running the code. The archived dataset is in the Harvard Forest Data Archive (dataset ID to be added); cite that dataset for the data.

## Citation and license

Code: MIT (`LICENSE`). Please cite the archived release (Zenodo DOI; see `CITATION.cff`) and Thompson et al. (2026), https://doi.org/10.64898/2026.06.10.730433.
