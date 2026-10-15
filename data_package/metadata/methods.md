# Study trees

The 60 trees of the stem methane study were selected on size alone (second and third quartiles of each species' stand DBH distribution), without regard to external signs of decay: ten trees of each species at each site where the species occurs (Black Gum Swamp: Nyssa sylvatica, Acer rubrum, Tsuga canadensis; EMS upland: Quercus rubra, A. rubrum, T. canadensis). Trees are identified by the tag numbers of the stem flux study. Tree locations, species and DBH are also given in the stem flux dataset (trees.csv).

# Sonic and electrical resistance tomography

Each tree was measured at breast height (1.37 m) between 2 July and 2 August 2024 (three trees were scanned twice) with the PiCUS system (Argus Electronic GmbH, Rostock, Germany): the PiCUS Sonic Tomograph for SoT and the TreeTronic 3 for ERT. For SoT, at least eight sensors were placed around the stem (sensor counts read from each tomogram in sot_sensor_counts.csv) and sound travel times between all sensor pairs were used to compute the velocity tomogram; structural loss is the percentage of the cross-section in non-brown (low-velocity) classes, and solid wood the PiCUS "solid wood %". For ERT, electrodes were inserted through the bark at the same height and resistance between electrode pairs was used to compute the resistivity tomogram. Three trees have incomplete SoT scans; 57 trees have both SoT and ERT.

# ERT image analysis

Eight resistivity metrics (mean, median, standard deviation, coefficient of variation, Gini coefficient, Shannon entropy, central moisture accumulation and radial gradient) were extracted from each absolute-scale ERT tomogram with the ERT image-analysis application of the Tree-Tomography repository (Thompson and Gewirtzman 2026, https://github.com/graceethompson/Tree-Tomography): the application detects the measured cross-section, maps pixel colours to resistivity by calibrating the colour bar on a log scale, and computes the metrics over the pixels inside the cross-section (ERT_application_results.csv).

# Hemlock validation set

Twelve eastern hemlocks in the Harvard Forest CTFS-ForestGEO plot (HF253; hemlock IDs are HF_ followed by the ForestGEO tag) were scanned with SoT and ERT at breast height on 30 August and 11 September 2024, six of them also at a lower and an upper height. An increment core was taken from each tree just below breast height (about 1 m) for gravimetric moisture content, expressed on an oven-dry mass basis, so values can exceed 100 % (hemlock_core_moisture.csv). ERT metrics of the hemlock scans (hemlock_ERT_results.csv; some scans were analysed more than once and are listed as exported) and SoT damage (hemlock_SoT_results.csv) are given per image.

# File naming

Tomogram files of the 60 study trees are named tree_date (e.g. 412_22vii24 = tree 412, 22 July 2024; _2 marks a second scan of the same tree). The file 420_30ii24 is a second scan of tree 420 made on 30 July 2024 (the date recorded by the instrument), so "ii" in its name stands for "vii". Hemlock files are HF_tag_height (Lower, DBH, Upper).

# Code

The analyses of the companion papers are reproduced by the Tree-Tomography repository (archived on Zenodo) and the tree-flux-2025 repository (https://github.com/jgewirtzman/tree-flux-2025).
