# Ported from montage.py -- exemplar paired SoT+ERT tomograms for each category.
# Run from the repo root: Rscript analysis/revision/scripts/montage.R
if (!l10n_info()[["UTF-8"]]) try(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)
suppressPackageStartupMessages(library(magick))
library(grid)
source("analysis/revision/scripts/revision_common.R")

ERT <- "images/main_ERT"; SOT <- "images/main_SoT"

# tree, category label, color, species, %dmg, mean res, CMA, PC1, one-line visual read
trees <- list(
  list(544, "I — No decay", "#26456e", "Q. rubra", 0, 672, 0.33, -2.15,
       "dry (red) throughout, no void → sound"),
  list(418, "II — Incipient", "#4d7f00", "A. rubrum", 0, 135, 0.41, 2.44,
       "wet (blue) interior but NO structural loss → incipient/wetwood"),
  list(217, "III — Active", "#7f4f00", "T. canadensis", 11, 128, 0.24, 5.01,
       "wet (blue) interior WITH structural loss → active decay"),
  list(880, "IV — Cavity", "#7f0000", "Q. rubra", 30, 408, 0.10, -0.42,
       "DRY (red) core + large SoT void → open/desiccated cavity"),
  list(420, "IV — Cavity", "#7f0000", "A. rubrum", 16, 337, 0.04, -1.08,
       "dry (red) central mass + crack → cavity"),
  list(893, "III vs IV (borderline)", "#555555", "T. canadensis", 12, 454, 0.40, -0.30,
       "composite labels \"cavity\", but ERT shows a WET blue core → really Active"),
  list(380, "III — Active (not cavity)", "#7f4f00", "A. rubrum", 34, 81, 0.05, 3.95,
       "cracked & 34% loss, but ERT uniformly WET → Active, not a cavity")
)

n <- length(trees); ncol <- 2; nrow_ <- ceiling(n / ncol)
W <- 15; H <- 4.7 * nrow_

draw_fit <- function(img, x, y, w, h, figw, figh) {
  info <- image_info(img)
  ar <- info$height / info$width
  wi <- w * figw; hi <- h * figh
  if (wi * ar <= hi) { dw <- wi; dh <- wi * ar } else { dh <- hi; dw <- hi / ar }
  grid.raster(as.raster(img), x = unit(x + w / 2, "npc"), y = unit(y + h / 2, "npc"),
              width = unit(dw, "in"), height = unit(dh, "in"))
}

out <- file.path(OUT_DIR, "CJFR-scan-montage.png")
png(out, width = W, height = H, units = "in", res = 115, type = "cairo")
grid.newpage()
grid.text(paste0("What the categories look like: paired SoT (structural loss) + ERT (moisture) tomograms\n",
                 "Cavity = DRY red core + void   ·   Active = WET blue interior + loss   ·   ",
                 "Incipient = WET interior, no loss   ·   No-decay = DRY, no loss"),
          y = 0.985, gp = gpar(fontsize = 12, lineheight = 1.25))
top <- 0.955
ch <- top / nrow_; cw <- 1 / ncol
for (i in seq_len(n)) {
  tr <- trees[[i]]
  tid <- tr[[1]]; cat_ <- tr[[2]]; col <- tr[[3]]; sp <- tr[[4]]
  dmg <- tr[[5]]; res <- tr[[6]]; cma <- tr[[7]]; pc1 <- tr[[8]]; read <- tr[[9]]
  r <- (i - 1) %/% ncol; k <- (i - 1) %% ncol
  x0 <- k * cw; ytop <- top - r * ch
  grid.text(sprintf("Tree %s  ·  %s\n%s  ·  SoT loss %s%%  ·  mean ρ %s Ω·m  ·  CMA %s  ·  PC1 %+.2f\n“%s”",
                    tid, cat_, sp, dmg, res, format(cma), pc1, read),
            x = x0 + cw / 2, y = ytop - 0.075 * ch,
            gp = gpar(fontsize = 9.2, col = col, fontface = "bold", lineheight = 1.4))
  panels <- list(list(find_scan(SOT, tid), "SoT (structure)"),
                 list(find_scan(ERT, tid), "ERT (moisture)"))
  for (j in seq_along(panels)) {
    px <- x0 + 0.035 * cw + (j - 1) * 0.475 * cw
    pw <- 0.455 * cw
    img <- image_read(panels[[j]][[1]])   # full frame, not cropped (as in montage.py)
    grid.text(panels[[j]][[2]], x = px + pw / 2, y = ytop - 0.185 * ch,
              gp = gpar(fontsize = 8, col = "#333333"))
    draw_fit(img, px, ytop - ch + 0.015 * ch, pw, 0.775 * ch, W, H)
  }
}
invisible(dev.off())
cat("saved montage with", n, "trees\n")
