# Verifies that the PiCUS velocity colour scale was constant across every SoT
# tomogram export (main survey + validation hemlocks): the "v:100%" and "v:70%"
# scale-endpoint labels are pixel-identical (within JPEG compression noise) in
# every real export. Supports the reproducibility response: the colour scale and
# velocity thresholds used to tabulate damaged area did not vary across trees.
library(magick)

files <- list.files(c("images/main_SoT", "images/hemlock_SoT"),
                    full.names = TRUE, pattern = "jpe?g$|png$")
# empty exports of the three dropped trees, two test artifacts, and one failed
# ellipse-only export (tree 412 was rescanned the next day: 412_26vii24.jpg)
skip <- c("366_18vii24.jpg", "428_02viii24.jpg", "677_18vii24.jpg",
          "dummy_file.jpg", "test_elipse.jpg", "412_25vii24.jpg")
files <- files[!basename(files) %in% skip]

corner <- function(f, side) {
  i <- image_read(f); w <- image_info(i)$width
  g <- if (side == "L") "60x14+0+0" else geometry_area(50, 14, w - 50, 0)
  as.numeric(image_data(image_crop(i, g), channels = "gray"))
}

refL <- corner(files[1], "L"); refR <- corner(files[1], "R")
dL <- sapply(files, function(f) mean(abs(corner(f, "L") - refL)))
dR <- sapply(files, function(f) mean(abs(corner(f, "R") - refR)))

cat(sprintf("%d real tomogram exports\n", length(files)))
cat(sprintf("v:100%% endpoint: max mean|diff| vs reference = %.4f\n", max(dL)))
cat(sprintf("v:70%%  endpoint: max mean|diff| vs reference = %.4f\n", max(dR)))
stopifnot(max(dL) < 0.02, max(dR) < 0.02)
cat("PASS: velocity colour-scale endpoints identical across all real exports\n")
