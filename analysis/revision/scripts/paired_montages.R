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
  plan <- list()
  for (ci in seq_len(nrow(catinfo))) {
    ids <- m$tree[cats == catinfo$code[ci]]
    plan[[length(plan) + 1]] <- list(kind = "hdr",
      txt = sprintf("%s — %s  (n=%d)", catinfo$code[ci], catinfo$nm[ci], length(ids)),
      col = catinfo$col[ci])
    if (length(ids) == 0) {
      plan[[length(plan) + 1]] <- list(kind = "img", ids = integer(0))
    } else {
      for (s in seq(1, length(ids), by = PER))
        plan[[length(plan) + 1]] <- list(kind = "img",
                                         ids = ids[s:min(s + PER - 1, length(ids))])
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
        t <- p$ids[k]
        base <- 0.01 + (k - 1) * (0.98 / PER); cellw <- 0.98 / PER
        row <- m[m$tree == t, ][1, ]
        dirs <- c(SOT, ERT)
        for (jj in seq_along(dirs)) {
          img <- image_crop(image_read(find_scan(dirs[jj], t)), CROP)
          draw_fit(img, base + (jj - 1) * (cellw * 0.48), ycur + 0.03 * h,
                   cellw * 0.45, h * 0.74, W, H)
        }
        grid.text(sprintf("#%s %s SoT%s ρ%.0f", t, substr(row$sp, 1, 4),
                          row$percent_damaged, row$mean),
                  x = base + cellw * 0.48, y = ycur + h * 0.88,
                  gp = gpar(fontsize = 6.6))
      }
    }
  }
  dev.off()
  cat("saved", out, "\n")
}

# Publication run regenerates only the published (PC1) scheme montage;
# call montage("spmed") or montage("6cell") by hand if those variants are needed.
for (k in "PC1") montage(k)
