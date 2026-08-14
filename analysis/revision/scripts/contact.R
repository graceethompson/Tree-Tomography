# Ported from contact.py -- contact sheets of all 57 main-tree ERT tomograms.
# Run from the repo root: Rscript analysis/revision/scripts/contact.R
if (!l10n_info()[["UTF-8"]]) try(invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)
suppressPackageStartupMessages(library(magick))
library(grid)
source("analysis/revision/scripts/revision_common.R")

ERT <- "images/main_ERT"
CROP <- "410x470+60+30"   # PIL crop box (60,30)-(470,500)

df <- load_main()
rows <- df[order(df$tree), ]

# Draw a magick image inside the box (x, y, w, h) [npc, lower-left corner],
# preserving the image aspect ratio; figw/figh are the device size in inches.
draw_fit <- function(img, x, y, w, h, figw, figh) {
  info <- image_info(img)
  ar <- info$height / info$width
  wi <- w * figw; hi <- h * figh
  if (wi * ar <= hi) { dw <- wi; dh <- wi * ar } else { dh <- hi; dw <- hi / ar }
  grid.raster(as.raster(img), x = unit(x + w / 2, "npc"), y = unit(y + h / 2, "npc"),
              width = unit(dw, "in"), height = unit(dh, "in"))
}

sheet <- function(sub, fname, title) {
  n <- nrow(sub); ncol <- 6; nr <- ceiling(n / ncol)
  W <- ncol * 2.4; H <- nr * 2.6
  png(fname, width = W, height = H, units = "in", res = 95, type = "cairo")
  grid.newpage()
  grid.text(title, y = 0.988, gp = gpar(fontsize = 12))
  top <- 0.965
  ch <- top / nr; cw <- 1 / ncol
  for (i in seq_len(n)) {
    r <- (i - 1) %/% ncol; k <- (i - 1) %% ncol
    x0 <- k * cw; ytop <- top - r * ch
    img <- image_crop(image_read(find_scan(ERT, sub$tree[i])), CROP)
    bw <- cw * 0.94; bh <- ch * 0.72
    draw_fit(img, x0 + (cw - bw) / 2, ytop - ch + 0.02 * ch, bw, bh, W, H)
    grid.text(sprintf("#%s %s\nSoT %s%%", sub$tree[i], substr(sub$sp[i], 1, 4),
                      sub$percent_damaged[i]),
              x = x0 + cw / 2, y = ytop - ch * 0.115,
              gp = gpar(fontsize = 8, lineheight = 1.1))
  }
  dev.off()
  cat("saved", fname, n, "trees\n")
}

sheet(rows[1:30, ], file.path(OUT_DIR, "contact1.png"),
      "ERT tomograms 1/2  (blue=wet, red=dry)")
sheet(rows[31:nrow(rows), ], file.path(OUT_DIR, "contact2.png"),
      "ERT tomograms 2/2  (blue=wet, red=dry)")
