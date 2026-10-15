# ============================================================
# data_package/build_data_package.R
# Assemble the tomography data package for the Harvard Forest Data Archive (on EDI):
# sonic (SoT) and electrical resistance (ERT) tomography of the 60 stem-flux trees and of
# the 12 hemlocks of the moisture validation set. A separate dataset from the stem flux
# package; both papers (Thompson et al., CJFR; Gewirtzman et al., JGR-B) cite it.
#
# Run from the repository root. Sources: data/, images/ and analysis/revision/data/ of this
# repository; metadata text and column definitions in data_package/metadata/.
# Output : data_package/output/package/ (EML templates in data_package/output/).
# Requires: EMLassemblyline, readxl
# ============================================================
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
suppressPackageStartupMessages(library(EMLassemblyline))

PACKAGE_ID <- "knb-lter-hfr.0.2"   # placeholder until the HF archive issues the ID
REPO <- "."
META <- "data_package/metadata"
TEMPLATE_DIR <- "data_package/output"; PACKAGE_DIR <- file.path(TEMPLATE_DIR, "package")
unlink(TEMPLATE_DIR, recursive = TRUE); dir.create(PACKAGE_DIR, recursive = TRUE)

# ---- tables -------------------------------------------------------------------
tables <- c(
  Tree_ID_info.csv                = file.path(REPO, "data", "Tree_ID_info.csv"),
  ERT_application_results.csv     = file.path(REPO, "data", "ERT_application_results.csv"),
  tomography_results_compiled.csv = file.path(REPO, "data", "tomography_results_compiled.csv"),
  sot_sensor_counts.csv           = file.path(REPO, "analysis", "revision", "data", "sot_image_annotations.csv"),
  hemlock_ERT_results.csv         = file.path(REPO, "data", "hemlock", "ERT_results_2026-03-05.csv"),
  hemlock_SoT_results.csv         = file.path(REPO, "data", "hemlock", "SOT_results.csv"))
stopifnot(all(file.exists(tables)))
for (k in names(tables)) {   # copied as recorded; byte-order marks removed so the column names are clean
  x <- readLines(tables[[k]], encoding = "UTF-8", warn = FALSE); x[1] <- sub("^﻿", "", x[1])
  writeLines(x, file.path(PACKAGE_DIR, k), useBytes = TRUE)
}
mc <- as.data.frame(readxl::read_excel(file.path(REPO, "data", "hemlock", "MC_Tomo_paper_Jon.xlsx")))
stopifnot(identical(names(mc), c("ID", "Moisture content (%)")))
write.csv(data.frame(tree_id = mc$ID, moisture_content_pct = mc$`Moisture content (%)`),
          file.path(PACKAGE_DIR, "hemlock_core_moisture.csv"), row.names = FALSE)
tables <- c(names(tables), "hemlock_core_moisture.csv")

# ---- attribute templates from the dictionaries ---------------------------------------------
units <- c("ohm m" = "ohmMeter", cm = "centimeter", "%" = "percent")
for (f in tables) {
  d <- read.csv(file.path(PACKAGE_DIR, f), stringsAsFactors = FALSE, check.names = FALSE)
  dict <- read.csv(file.path(META, "dictionaries", f), stringsAsFactors = FALSE, check.names = FALSE)
  if (!setequal(names(d), dict$column)) stop(f, ": dictionary and table differ: ",
    paste(c(setdiff(names(d), dict$column), setdiff(dict$column, names(d))), collapse = ", "))
  dict <- dict[match(names(d), dict$column), ]
  ids <- c("tree", "Tree", "tree_id")
  cls <- vapply(names(d), function(v) if (is.numeric(d[[v]]) && !v %in% ids) "numeric" else "character", "")
  un <- ifelse(is.na(dict$unit), "", dict$unit)
  u <- ifelse(cls == "numeric", ifelse(nzchar(un), unname(units[un]), "dimensionless"), "")
  if (anyNA(u)) stop(f, ": unit not mapped")
  write.table(data.frame(attributeName = names(d), attributeDefinition = dict$description, class = cls, unit = u,
                         dateTimeFormatString = "", missingValueCode = "NA", missingValueCodeExplanation = "not available"),
              file.path(TEMPLATE_DIR, paste0("attributes_", sub("\\.csv$", "", f), ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
file.copy(file.path(META, c("abstract.md", "methods.md", "keywords.txt", "personnel.txt", "intellectual_rights.txt",
                            "geographic_coverage.txt", "custom_units.txt")), TEMPLATE_DIR, overwrite = TRUE)

# ---- tomogram images, zipped (folder per image type; tests and duplicate copies left out) ----
zip_images <- function(zipname, folders) {
  stage <- file.path(tempdir(), sub("\\.zip$", "", zipname)); unlink(stage, recursive = TRUE)
  n <- 0
  for (dest in names(folders)) {
    f <- list.files(file.path(REPO, "images", folders[[dest]]), pattern = "\\.jpg$", full.names = TRUE)
    dir.create(file.path(stage, dest), recursive = TRUE); file.copy(f, file.path(stage, dest), copy.date = TRUE); n <- n + length(f)
  }
  old <- setwd(stage); on.exit(setwd(old))
  zip(normalizePath(file.path(old, PACKAGE_DIR, zipname), mustWork = FALSE), list.files(".", recursive = TRUE), flags = "-q")
  message("  ", zipname, ": ", n, " images")
}
zip_images("tomograms_study_trees.zip", c(ERT_absolute = "main_ERT", ERT_normalized = "main_ERT_normalized", SoT = "main_SoT"))
zip_images("tomograms_hemlock.zip", c(ERT = "hemlock_ERT", SoT = "hemlock_SoT"))

table_info <- rbind(
  c("Tomography tree table", "SoT results of each of the 60 study trees: plot, species, tag, DBH, crack detected, percent solid and damaged wood."),
  c("ERT image analysis results", "Eight ERT resistivity metrics for each ERT tomogram of the study trees, as exported by the ERT image-analysis application (Tree-Tomography repository)."),
  c("Tomography metrics (compiled)", "The ERT metrics and SoT solid and damaged percentages of each study tree in one table (as used by the stem flux study)."),
  c("SoT sensor counts", "Number of SoT sensors and PiCUS solid-wood percentage read from each SoT tomogram (study trees and hemlocks)."),
  c("Hemlock ERT results", "Eight ERT resistivity metrics for each ERT tomogram of the hemlock validation set (up to three heights per tree)."),
  c("Hemlock SoT results", "SoT damaged percentage for each SoT tomogram of the hemlock validation set."),
  c("Hemlock core moisture", "Gravimetric moisture content of the increment core of each hemlock in the validation set."))
other_info <- rbind(
  c("tomograms_study_trees.zip", "Tomograms of the 60 study trees",
    "JPEG exports of the PiCUS tomograms at breast height: ERT_absolute (ERT, absolute resistivity scale), ERT_normalized (ERT, scale normalized to each tree) and SoT (sonic). File names are tree_date (e.g. 412_22vii24 = tree 412, 22 July 2024)."),
  c("tomograms_hemlock.zip", "Tomograms of the hemlock validation set",
    "JPEG exports of the ERT and SoT tomograms of the 12 validation hemlocks. File names are HF_tag_height (Lower, DBH, Upper)."))

eml <- make_eml(
  path = TEMPLATE_DIR, data.path = PACKAGE_DIR, eml.path = PACKAGE_DIR,
  dataset.title = "Sonic and electrical resistance tomography of upland and wetland trees at Harvard Forest 2024, with a hemlock moisture validation set",
  temporal.coverage = c("2024-07-02", "2024-09-11"),   # scan dates (Zeit field of the PiCUS .pit files)
  maintenance.description = "Completed.",
  data.table = tables, data.table.name = table_info[, 1], data.table.description = table_info[, 2],
  data.table.quote.character = rep('"', length(tables)),
  other.entity = other_info[, 1], other.entity.name = other_info[, 2], other.entity.description = other_info[, 3],
  user.id = "jgewirtzman", user.domain = "EDI", package.id = PACKAGE_ID, write.file = TRUE, return.obj = TRUE)
message("EML: ", file.path(PACKAGE_DIR, paste0(PACKAGE_ID, ".xml")))
