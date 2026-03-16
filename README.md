# Tree-Tomography

Natural history study of internal tree decay within and amongst dominant tree species present along an environmental gradient at the Harvard Forest.

## Study design

### SoT and ERT measurements

Corresponding to an existing biometric plot network established around the footprint of an eddy covariance tower at Harvard Forest, we defined our upland site as EMS. Our wetland site was the Black Gum Swamp, which is a forested peatland within the same Prospect Hill tract as the EMS site. We selected four tree species common in either the EMS or BGS site, or both. The generalist species were red maple (_Acer rubrum_) and eastern hemlock (_Tsuga canadensis_). The upland specialist species was northern red oak (_Quercus rubra_), and the wetland specialist species was black gum (_Nyssa sylvatica_). We selected 10 trees per species for observation. Each tree had to be representative of the average size class (dbh) for its respective species estimated using Harvard Forest's ForestGEO data. Sonic and electric resistivity tomography was performed once at each tree. Some electric resistivity tomograms were not fully uploaded at the end of the study, and the tree count at each site was as follows: 10 wetland red maples, 10 upland red maples, 9 red oaks, 10 black gums, 9 wetland eastern hemlocks, and 9 upland eastern hemlocks.

### Eastern hemlock water content

Twelve eastern hemlocks were measured at multiple heights with ERT to validate ERT-derived moisture metrics against core-extracted gravimetric moisture content.

## Repository structure

### `data/`

- `Tree_ID_info.csv` — 57 main study trees: plot, species, tree tag, dbh, SoT percent damaged/solid wood
- `ERT_application_results.csv` — 8 ERT metrics (mean, median, sd, cv, gini, entropy, cma, radial gradient) per tree
- `hemlock/` — Hemlock validation data
  - `validation_summary.csv` — 12 hemlock trees with ERT conductance metrics and core moisture
  - `SOT_results.csv` — SoT percent damaged per hemlock section
  - `ERT_results_2026-03-05.csv` — Per-height ERT metrics for hemlocks
  - `MC_Tomo_paper_Jon.xlsx` — Core moisture content data
- `supplementary/` — Supporting data files from earlier analyses

### `code/`

- `final_phase_and_scans.R` — Main analysis: 2x2 decay phase diagram (species-normalized PCA on ERT metrics vs SoT structural loss), PCA biplot, quadrant distributions, hemlock scan panels
- `phase_diagram_analysis.R` — Systematic search across 714 ERT axis / threshold combinations
- `phase_image_panels.R` — Tomogram image panels arranged by phase diagram position
- `treedecay_og.R` — SoT descriptive statistics and species/site comparisons (boxplots, ANOVAs)
- `treedecay_app.R` — ERT descriptive statistics and species/site comparisons (boxplots, ANOVAs)
- `treedecay_relfreq.R` — Relative frequency analysis
- `ert_validation_analysis.R` — Hemlock ERT validation against core moisture
- `check_agreement.R` — Agreement analysis between SoT and ERT
- `Tree_Tomography_Data_Analysis.Rmd` — Data compilation (Max Lutz)

### `images/`

- `main_ERT/` — Main study ERT tomogram images
- `main_SoT/` — Main study SoT tomogram images
- `main_ERT_normalized/` — Normalized ERT images
- `hemlock_ERT/` — Hemlock ERT tomogram images
- `hemlock_SoT/` — Hemlock SoT tomogram images

### `output/`

- `figures/` — Main study figures (phase diagram, PCA biplot, quadrant distributions, SoT/ERT boxplots)
- `hemlock_figures/` — Hemlock validation and scan figures
- `panels/` — Tomogram image panels by phase diagram position
- `tables/` — Analysis result tables

## Running the analysis

All scripts use relative paths from the project root. Set your working directory to the repository root before running. The primary analysis pipeline is:

1. `treedecay_og.R` — SoT descriptive analysis
2. `treedecay_app.R` — ERT descriptive analysis
3. `final_phase_and_scans.R` — Phase diagram and PCA analysis
