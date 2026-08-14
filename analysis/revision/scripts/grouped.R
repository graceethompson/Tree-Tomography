# Ported from grouped.py -- ERT scan montages grouped by category under four
# alternative binning schemes. Run from the repo root:
#   Rscript analysis/revision/scripts/grouped.R
if (!l10n_info()[["UTF-8"]]) try(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)
suppressPackageStartupMessages(library(magick))
library(grid)
source("analysis/revision/scripts/revision_common.R")

ERT <- "images/main_ERT"
CROP <- "410x470+60+30"

m <- load_merged()
b <- build_pc1(m)          # species-normalized PCA, PC1 oriented higher = wetter
pc <- b$pc1
m$struct <- m$percent_damaged > 1
m$res_dev <- ave(m$mean, m$sp, FUN = function(x) (x - median(x)) / sd(x))

schemes <- list(
  A_PC1    = list(title = "Bin by PC1 (original)", anom = pc > mean(pc)),
  B_CMA    = list(title = "Bin by CMA > 0.33", anom = m$cma > 0.33),
  C_absres = list(title = "Bin by absolute resistivity < sample median",
                  anom = m$mean < median(m$mean)),
  D_spmed  = list(title = "Bin by within-species median split (wetter than species median)",
                  anom = m$res_dev < 0)
)

catinfo <- data.frame(code = c("I", "II", "III", "IV"),
                      nm = c("No decay", "Incipient", "Active", "Cavity"),
                      col = c("#3b6fb0", "#4e9a2c", "#d98a1f", "#b83232"))
NC <- 8

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
  cats <- four_cat(schemes[[key]]$anom, m$struct)
  # row plan: header rows + image rows of up to NC scans
  plan <- list()
  for (ci in seq_len(nrow(catinfo))) {
    ids <- m$tree[cats == catinfo$code[ci]]
    plan[[length(plan) + 1]] <- list(kind = "hdr",
      txt = sprintf("%s — %s   (n=%d)", catinfo$code[ci], catinfo$nm[ci], length(ids)),
      col = catinfo$col[ci])
    if (length(ids) == 0) {
      plan[[length(plan) + 1]] <- list(kind = "img", ids = integer(0))
    } else {
      for (s in seq(1, length(ids), by = NC))
        plan[[length(plan) + 1]] <- list(kind = "img",
                                         ids = ids[s:min(s + NC - 1, length(ids))])
    }
  }
  units <- vapply(plan, function(p) if (p$kind == "hdr") 0.34 else 1.0, numeric(1))
  tot <- sum(units)
  W <- NC * 1.7; H <- tot * 1.55
  out <- file.path(OUT_DIR, paste0("CJFR-grouped-", key, ".png"))
  png(out, width = W, height = H, units = "in", res = 85, type = "cairo")
  grid.newpage()
  grid.text(paste0(title, "   —   ERT scans grouped by resulting category"),
            y = 0.988, gp = gpar(fontsize = 13))
  top <- 0.97
  ycur <- top
  for (pi in seq_along(plan)) {
    p <- plan[[pi]]
    h <- units[pi] / tot * top
    ycur <- ycur - h
    if (p$kind == "hdr") {
      grid.rect(x = 0.5, y = ycur + h / 2, width = 0.98, height = h - 0.004,
                gp = gpar(fill = adjustcolor(p$col, alpha.f = 0.16), col = NA))
      grid.text(p$txt, x = 0.02, y = ycur + h / 2, just = "left",
                gp = gpar(fontsize = 12, fontface = "bold", col = p$col))
    } else {
      for (k in seq_along(p$ids)) {
        t <- p$ids[k]
        x <- 0.01 + (k - 1) * (0.98 / NC); w <- 0.98 / NC * 0.94
        img <- image_crop(image_read(find_scan(ERT, t)), CROP)
        draw_fit(img, x, ycur + 0.01 * h, w, h * 0.72, W, H)
        row <- m[m$tree == t, ][1, ]
        grid.text(sprintf("#%s %s\nSoT%s ρ%.0f", t, substr(row$sp, 1, 4),
                          row$percent_damaged, row$mean),
                  x = x + w / 2, y = ycur + h * 0.755, just = c("centre", "bottom"),
                  gp = gpar(fontsize = 6.3, lineheight = 1.05))
      }
    }
  }
  dev.off()
  cat("saved", out, "\n")
}

for (k in names(schemes)) montage(k)
