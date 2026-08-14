# Ported from flow.py -- decision-tree flowchart of the decay classification
# rules and thresholds. Run from the repo root:
#   Rscript analysis/revision/scripts/flow.R
if (!l10n_info()[["UTF-8"]]) try(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)
library(grid)
source("analysis/revision/scripts/revision_common.R")

out <- file.path(OUT_DIR, "CJFR-decision-tree.png")
png(out, width = 13.5, height = 9, units = "in", res = 150, type = "cairo")
grid.newpage()
# data coordinates 0-100 x 0-100 as in the Python script, plus headroom for the title
pushViewport(viewport(xscale = c(0, 100), yscale = c(-1, 104)))

box <- function(x, y, w, h, text, fc, ec = "#333333", fs = 9.5, tc = "#111111", bold = FALSE) {
  pad <- 1.4   # approximates FancyBboxPatch round,pad=0.6
  grid.roundrect(x = unit(x, "native"), y = unit(y, "native"),
                 width = unit(w + pad, "native"), height = unit(h + pad, "native"),
                 r = unit(2, "mm"), gp = gpar(fill = fc, col = ec, lwd = 1.7))
  grid.text(text, x = unit(x, "native"), y = unit(y, "native"),
            gp = gpar(fontsize = fs, col = tc, lineheight = 1.25,
                      fontface = if (bold) "bold" else "plain"))
}
arrowf <- function(x1, y1, x2, y2, label = NULL, lx = 0, ly = 0, col = "#444444") {
  grid.lines(x = unit(c(x1, x2), "native"), y = unit(c(y1, y2), "native"),
             arrow = arrow(type = "closed", angle = 22, length = unit(3, "mm")),
             gp = gpar(col = col, fill = col, lwd = 1.8))
  if (!is.null(label))
    grid.text(label, x = unit((x1 + x2) / 2 + lx, "native"),
              y = unit((y1 + y2) / 2 + ly, "native"),
              gp = gpar(fontsize = 8.5, col = col, fontface = "bold"))
}
note <- function(x, y, text, fs = 7, col = "#555555") {
  grid.text(text, x = unit(x, "native"), y = unit(y, "native"),
            gp = gpar(fontsize = fs, col = col, fontface = "italic", lineheight = 1.25))
}

grid.text("Decay classification — decision rules & thresholds",
          x = unit(50, "native"), y = unit(102.3, "native"),
          gp = gpar(fontsize = 13, fontface = "bold"))

# input
box(50, 95, 54, 7,
    "Each tree:  SoT % structural damage  +  ERT resistivity metrics (8, incl. mean ρ and CMA)",
    "#eef2f7", fs = 9)
arrowf(50, 91.5, 50, 88)
# split 1: SoT
box(50, 84, 40, 8,
    "AXIS 1 — STRUCTURAL LOSS  (SoT)\nIs SoT % damage  >  1% ?",
    "#d6e4f0", fs = 9.5, bold = TRUE)
note(50, 77.5,
     "crisp / bimodal: 0% or ≥6%, empty gap between  ·  stable for any cut in 1–5% (2 trees change)",
     fs = 7.3)

# NO branch -> moisture anomaly (I/II)
arrowf(31, 84, 16, 70, "NO  (no structural loss)", lx = -3, ly = 6, col = "#26456e")
box(16, 66, 30, 9,
    "AXIS 2a — MOISTURE ANOMALY\nspecies-normalized ERT PC1\n(continuum — graded)",
    "#e8f0e2", fs = 8.8, bold = TRUE)
# graded outcomes
arrowf(16, 61.5, 10, 50)
box(10, 45, 17, 8, "I — No decay\nPC1 ≤ mean (0)",
    "#dbe4ef", fs = 8.5, tc = "#26456e", bold = TRUE)
arrowf(16, 61.5, 27, 50)
box(30, 45, 20, 10,
    "II — Incipient\n\"possible\": 0 < PC1 ≤ 1.0\n\"confident\": PC1 > 1.0",
    "#e5efd6", fs = 8.2, tc = "#4d7f00", bold = TRUE)
note(20, 37,
     "≈ 30 normal / 7 possible / 5 confident (of 39 sound trees)\nincipient ≠ decay: could be early decay OR bacterial wetwood")

# YES branch -> core wetness (III/IV)
arrowf(69, 84, 84, 70, "YES  (structural loss present)", lx = 6, ly = 6, col = "#7f4f00")
box(84, 66, 30, 9,
    "AXIS 2b — CORE WETNESS\nis the damaged core WET or DRY ?",
    "#f3e6d6", fs = 8.8, bold = TRUE)
arrowf(84, 61.5, 72, 50)
box(69, 44, 26, 11,
    "III — Active decay\nWET core:\nmean ρ < 300 Ω·m  OR  CMA ≥ 0.15",
    "#f0e2cf", fs = 8, tc = "#7f4f00", bold = TRUE)
arrowf(84, 61.5, 95, 50)
box(92, 44, 15, 11,
    "IV — Cavity\nDRY core:\nmean ρ ≥ 300\nAND CMA < 0.15",
    "#f0d6d6", fs = 7.8, tc = "#7f0000", bold = TRUE)
note(84, 35.5,
     "physically-grounded refinement of the composite cut\n(reclassifies 893: wet core → Active, not Cavity)")

# sensitivity footer band
box(50, 22, 92, 12, paste0(
  "SENSITIVITY / STATUS\n",
  "• SoT 1% threshold: ROBUST — reclassifies only 2/57 trees across 1–5%; report equipment detection floor.\n",
  "• PC1 anomaly cut (I↔II): OPERATIONAL CONVENTION on a continuum — ±0.5 SD (±0.95) reclassifies 11–14 trees; ",
  "report graded bands + full sweep. Direction (BGS>EMS) invariant to threshold & normalization.\n",
  "• Core-wetness cut (III↔IV): PHYSICALLY grounded but thresholds (300 Ω·m, CMA 0.15) are data-suggested, not validated."),
  "#fbfbf6", ec = "#aaaaaa", fs = 8, tc = "#222222")
box(50, 10, 92, 7, paste0(
  "ABSOLUTE / FUTURE:  PC1 mean ≈ 369 Ω·m ≈ 103% predicted moisture; Otsu break ≈ 387 Ω·m ≈ 101%.  ",
  "Absolute cut premature — n=12 calibration 95% band ≈ 5–800 Ω·m.\n",
  "True incipient-decay vs wetwood, and per-species absolute thresholds, require destructive discs/cross-sections."),
  "#f0f0ea", ec = "#aaaaaa", fs = 7.6, tc = "#333333")

popViewport()
invisible(dev.off())
cat("saved flowchart\n")
