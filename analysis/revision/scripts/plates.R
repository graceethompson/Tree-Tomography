# Fig. S18 plates: paired SoT + ERT tomograms for every tree, grouped by the
# manuscript's decay category (PC1 scheme): one plate per category, each fitting
# one portrait page (category I at 4 pairs per row, the others at 3). The 12 cored
# hemlocks (DBH scans, classified by the same thresholds) are appended to
# their categories, labelled "(cored)".
# Run from the repo root: Rscript analysis/revision/scripts/plates.R
if (!l10n_info()[["UTF-8"]]) try(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)
suppressPackageStartupMessages(library(magick)); library(grid)
source("analysis/revision/scripts/revision_common.R")
ERT <- "images/main_ERT"; SOT <- "images/main_SoT"; CROP <- "410x470+60+30"
m <- load_merged(); pc <- build_pc1(m)$pc1; m$struct <- m$percent_damaged > 1
cats <- four_cat(pc > mean(pc), m$struct)
sp_full <- function(s) sub("^([A-Z])\\.", "\\1. ", s)
recs <- data.frame(label = sprintf("%s  %s\n%g%% damage  |  %.0f Ω·m", m$tree, sp_full(m$sp), m$percent_damaged, m$mean),
                   sot = vapply(m$tree, function(t) find_scan(SOT, t), character(1)),
                   ert = vapply(m$tree, function(t) find_scan(ERT, t), character(1)),
                   cat = cats, cored = FALSE, stringsAsFactors = FALSE)
v <- read.csv("data/hemlock/validation_phases.csv")
recs <- rbind(recs, data.frame(label = sprintf("%s (cored)\n%g%% damage  |  %.0f Ω·m", v$tree, v$percent_damaged, v$mean),
                               sot = file.path("images/hemlock_SoT", paste0(v$tree, "_DBH.jpg")),
                               ert = file.path("images/hemlock_ERT", paste0(v$tree, "_DBH.jpg")),
                               cat = sub(":.*$", "", v$quadrant), cored = TRUE, stringsAsFactors = FALSE))
catinfo <- data.frame(code = c("I", "II", "III", "IV"), nm = c("No Decay", "Incipient", "Active", "Cavity"),
                      col = c("#3b6fb0", "#4e9a2c", "#d98a1f", "#b83232"))
PER_FOR <- c(I = 4, II = 3, III = 3, IV = 3); ROWH_FOR <- c(`4` = 1.0, `3` = 1.2)
draw_fit <- function(img, x, y, w, h, figw, figh) {
  info <- image_info(img); ar <- info$height / info$width; wi <- w * figw; hi <- h * figh
  if (wi * ar <= hi) { dw <- wi; dh <- wi * ar } else { dh <- hi; dw <- hi / ar }
  grid.raster(as.raster(img), x = unit(x + w / 2, "npc"), y = unit(y + h / 2, "npc"), width = unit(dw, "in"), height = unit(dh, "in"))
}
plates <- list(); part <- 0
for (ci in seq_len(nrow(catinfo))) {
  idx <- which(recs$cat == catinfo$code[ci]); idx <- idx[order(recs$cored[idx])]
  part <- part + 1
  nv <- sum(recs$cored[idx]); ns <- length(idx) - nv
  hdr <- sprintf("(%s)  Category %s — %s: %d study tree%s%s", letters[part], catinfo$code[ci], catinfo$nm[ci],
                 ns, if (ns == 1) "" else "s", if (nv > 0) sprintf(" + %d cored hemlock%s", nv, if (nv == 1) "" else "s") else "")
  plates[[part]] <- list(hdr = hdr, col = catinfo$col[ci], ids = idx, per = PER_FOR[[catinfo$code[ci]]])
}
for (pi in seq_along(plates)) {
  p <- plates[[pi]]; PER <- p$per; nrow_ <- ceiling(length(p$ids) / PER)
  W <- 6.5; rowh <- ROWH_FOR[[as.character(PER)]]; H <- 0.42 + nrow_ * rowh
  out <- file.path(OUT_DIR, sprintf("CJFR-plate-%s.png", letters[pi]))
  png(out, width = W, height = H, units = "in", res = 220, type = "cairo"); grid.newpage()
  hh <- 0.34 / H
  grid.rect(x = 0.5, y = 1 - hh / 2 - 0.01, width = 0.98, height = hh, gp = gpar(fill = adjustcolor(p$col, alpha.f = 0.16), col = NA))
  grid.text(p$hdr, x = 0.02, y = 1 - hh / 2 - 0.01, just = "left", gp = gpar(fontsize = 10.5, fontface = "bold", col = p$col))
  top <- 1 - hh - 0.02; rh <- rowh / H
  for (k in seq_along(p$ids)) {
    r <- recs[p$ids[k], ]; row <- (k - 1) %/% PER; col <- (k - 1) %% PER
    y0 <- top - (row + 1) * rh; base <- 0.01 + col * (0.98 / PER); cellw <- 0.98 / PER
    for (jj in 1:2) {
      img <- image_crop(image_read(c(r$sot, r$ert)[jj]), CROP)
      draw_fit(img, base + (jj - 1) * (cellw * 0.49), y0 + 0.02 * rh, cellw * 0.47, rh * 0.74, W, H)
    }
    grid.text(r$label, x = base + cellw * 0.49, y = y0 + rh * 0.885,
              gp = gpar(fontsize = 6.6, lineheight = 1.05, col = if (r$cored) "#555555" else "black", fontface = if (r$cored) "italic" else "plain"))
  }
  dev.off(); cat("saved", out, sprintf(" (%d pairs, %.1f in tall)\n", length(p$ids), H))
}
