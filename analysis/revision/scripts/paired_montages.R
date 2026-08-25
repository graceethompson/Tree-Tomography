# Ported from paired_montages.py -- paired SoT+ERT scan montages, grouped by
# classification scheme. Run from the repo root:
#   Rscript analysis/revision/scripts/paired_montages.R
if (!l10n_info()[["UTF-8"]]) try(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)
suppressPackageStartupMessages(library(magick))
library(grid)
source("analysis/revision/scripts/revision_common.R")

ERT <- "images/main_ERT"; SOT <- "images/main_SoT"
CROP <- "410x470+60+30"

m <- load_merged()
b <- build_pc1(m)          # species-normalized PCA, PC1 oriented higher = wetter
pc <- b$pc1
m$struct <- m$percent_damaged > 1
m$xdev <- xdev_species_median(m)

catinfo <- data.frame(code = c("I", "II", "III", "IV"),
                      nm = c("No decay", "Incipient", "Active", "Cavity"),
                      col = c("#3b6fb0", "#4e9a2c", "#d98a1f", "#b83232"))

schemes <- list(
  PC1 = list(title = "Bin by PC1 (published scheme)",
             cats = four_cat(pc > mean(pc), m$struct)),
  spmed = list(title = "Bin by within-species median split",
               cats = four_cat(m$xdev > 0, m$struct)),
  `6cell` = list(title = "6-cell: species-median +/-0.5 SD dead band (defaults to row baseline)",
                 cats = ifelse(!m$struct & (m$xdev > 0.5), "II",
                        ifelse(!m$struct & !(m$xdev > 0.5), "I",
                        ifelse(m$struct & (m$xdev < -0.5), "IV", "III"))))
)
PER <- 4

draw_fit <- function(img, x, y, w, h, figw, figh) {
  info <- image_info(img)
  ar <- info$height / info$width
  wi <- w * figw; hi <- h * figh
  if (wi * ar <= hi) { dw <- wi; dh <- wi * ar } else { dh <- hi; dw <- hi / ar }
  grid.raster(as.raster(img), x = unit(x + w / 2, "npc"), y = unit(y + h / 2, "npc"),
              width = unit(dw, "in"), height = unit(dh, "in"))
}

montage <- function(key) {
  title <- schemes[[key]]$title
  cats <- schemes[[key]]$cats
  # one record per displayed tree: label + image paths + category
  recs <- data.frame(
    label = sprintf("#%s %s SoT%s ρ%.0f", m$tree, substr(m$sp, 1, 4),
                    m$percent_damaged, m$mean),
    sot = vapply(m$tree, function(t) find_scan(SOT, t), character(1)),
    ert = vapply(m$tree, function(t) find_scan(ERT, t), character(1)),
    cat = cats, val = FALSE, stringsAsFactors = FALSE)
  # For the published scheme, append the 12 validation hemlocks (DBH scans),
  # classified by the identical thresholds (data/hemlock/validation_phases.csv,
  # written by code/final_phase_and_scans.R).
  vp <- "data/hemlock/validation_phases.csv"
  if (key == "PC1" && file.exists(vp)) {
    v <- read.csv(vp)
    recs <- rbind(recs, data.frame(
      label = sprintf("%s (val) SoT%s ρ%.0f", v$tree,
                      v$percent_damaged, v$mean),
      sot = file.path("images/hemlock_SoT", paste0(v$tree, "_DBH.jpg")),
      ert = file.path("images/hemlock_ERT", paste0(v$tree, "_DBH.jpg")),
      cat = sub(":.*$", "", v$quadrant), val = TRUE,
      stringsAsFactors = FALSE))
  }
  plan <- list()
  for (ci in seq_len(nrow(catinfo))) {
    idx <- which(recs$cat == catinfo$code[ci])
    idx <- idx[order(recs$val[idx])]   # study trees first, validation last
    nv <- sum(recs$val[idx])
    txt <- sprintf("%s — %s  (n=%d%s)", catinfo$code[ci], catinfo$nm[ci],
                   length(idx) - nv,
                   if (nv > 0) sprintf(" + %d validation", nv) else "")
    plan[[length(plan) + 1]] <- list(kind = "hdr", txt = txt,
                                     col = catinfo$col[ci])
    if (length(idx) == 0) {
      plan[[length(plan) + 1]] <- list(kind = "img", ids = integer(0))
    } else {
      for (s in seq(1, length(idx), by = PER))
        plan[[length(plan) + 1]] <- list(kind = "img",
                                         ids = idx[s:min(s + PER - 1, length(idx))])
    }
  }
  units <- vapply(plan, function(p) if (p$kind == "hdr") 0.32 else 1.0, numeric(1))
  tot <- sum(units)
  W <- PER * 3.5; H <- tot * 1.7
  out <- file.path(OUT_DIR, paste0("CJFR-paired-", key, ".png"))
  png(out, width = W, height = H, units = "in", res = 150, type = "cairo")
  grid.newpage()
  top <- 0.985
  ycur <- top
  for (pi in seq_along(plan)) {
    p <- plan[[pi]]
    h <- units[pi] / tot * top
    ycur <- ycur - h
    if (p$kind == "hdr") {
      grid.rect(x = 0.5, y = ycur + h / 2, width = 0.98, height = h - 0.003,
                gp = gpar(fill = adjustcolor(p$col, alpha.f = 0.16), col = NA))
      grid.text(p$txt, x = 0.02, y = ycur + h / 2, just = "left",
                gp = gpar(fontsize = 12, fontface = "bold", col = p$col))
    } else {
      for (k in seq_along(p$ids)) {
        r <- recs[p$ids[k], ]
        base <- 0.01 + (k - 1) * (0.98 / PER); cellw <- 0.98 / PER
        paths <- c(r$sot, r$ert)
        for (jj in seq_along(paths)) {
          img <- image_crop(image_read(paths[jj]), CROP)
          draw_fit(img, base + (jj - 1) * (cellw * 0.48), ycur + 0.03 * h,
                   cellw * 0.45, h * 0.74, W, H)
        }
        grid.text(r$label, x = base + cellw * 0.48, y = ycur + h * 0.88,
                  gp = gpar(fontsize = 6.6,
                            col = if (r$val) "#555555" else "black",
                            fontface = if (r$val) "italic" else "plain"))
      }
    }
  }
  dev.off()
  cat("saved", out, "\n")
}

# Publication run regenerates only the published (PC1) scheme montage;
# call montage("spmed") or montage("6cell") by hand if those variants are needed.
for (k in "PC1") montage(k)
