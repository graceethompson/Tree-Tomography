# fig_ert_workflow.R
# 3-panel figure illustrating the ERT image analysis workflow
# Panel A: Raw ERT tomogram
# Panel B: Boundary delineation & calibration (cross-section + calibration curve)
# Panel C: Quantitative output (radial profile + metrics table)
#
# Uses tree 876 (Quercus rubra) with strong central moisture anomaly (CMA = 0.48)

library(jpeg)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)
library(sp)

# ============================================================================
# Processing functions (from ERT_App.R)
# ============================================================================

arc_length_param <- function(r, g, b) {
  n <- length(r)
  if (n < 2) return(rep(0, n))
  d <- sqrt(diff(r)^2 + diff(g)^2 + diff(b)^2)
  s <- c(0, cumsum(d))
  if (max(s) > 0) s / max(s) else s
}

extract_colorbar <- function(img) {
  W <- ncol(img)
  rows <- 3:8
  r <- colMeans(img[rows, , 1])
  g <- colMeans(img[rows, , 2])
  b <- colMeans(img[rows, , 3])
  trim <- 3
  idx <- (trim):(length(r) - trim + 1)
  r <- r[idx]; g <- g[idx]; b <- b[idx]
  smooth5 <- function(v) {
    n <- length(v); out <- numeric(n)
    for (i in 1:n) { lo <- max(1, i-2); hi <- min(n, i+2); out[i] <- mean(v[lo:hi]) }
    out
  }
  r <- smooth5(r); g <- smooth5(g); b <- smooth5(b)
  if (length(r) > 256) {
    idx <- round(seq(1, length(r), length.out = 256))
    r <- r[idx]; g <- g[idx]; b <- b[idx]
  }
  if ((b[1] - r[1]) + (r[length(r)] - b[length(r)]) < 0) {
    r <- rev(r); g <- rev(g); b <- rev(b)
  }
  t_param <- arc_length_param(r, g, b)
  list(r = r, g = g, b = b, t = t_param)
}

create_log_calibration <- function(positions, values) {
  ord <- order(positions)
  positions <- positions[ord]; values <- values[ord]
  f <- approxfun(positions, log10(values), rule = 2)
  function(t) 10^(f(t))
}

map_colors_to_values <- function(px_r, px_g, px_b, cb, calib_func) {
  cb_r <- cb$r; cb_g <- cb$g; cb_b <- cb$b; cb_t <- cb$t
  n_px <- length(px_r)
  values <- numeric(n_px)
  chunk <- 5000
  for (start in seq(1, n_px, by = chunk)) {
    end <- min(start + chunk - 1, n_px)
    idx <- start:end
    dr <- outer(px_r[idx], cb_r, "-")
    dg <- outer(px_g[idx], cb_g, "-")
    db <- outer(px_b[idx], cb_b, "-")
    dist2 <- dr^2 + dg^2 + db^2
    best <- apply(dist2, 1, which.min)
    values[idx] <- calib_func(cb_t[best])
  }
  values
}

# Auto-detect PiCUS blue boundary polygon (same logic as ERT_App.R)
auto_detect_polygon_jpeg <- function(img) {
  H <- nrow(img); W <- ncol(img)
  colorbar_cutoff <- 40

  # Detect blue boundary pixels below colorbar
  rows <- (colorbar_cutoff + 1):H
  r_vals <- as.vector(img[rows, , 1])
  g_vals <- as.vector(img[rows, , 2])
  b_vals <- as.vector(img[rows, , 3])
  grid_xy <- expand.grid(row = rows, col = 1:W)
  is_blue <- (b_vals > 0.4) & (r_vals < 0.3) & (g_vals < 0.4) & ((b_vals - g_vals) > 0.15)
  blue_pts <- data.frame(row = grid_xy$row[is_blue], col = grid_xy$col[is_blue])
  cat("Blue polygon pixels detected:", nrow(blue_pts), "\n")

  cx <- median(blue_pts$col); cy <- median(blue_pts$row)
  angles <- atan2(blue_pts$row - cy, blue_pts$col - cx)
  dists  <- sqrt((blue_pts$col - cx)^2 + (blue_pts$row - cy)^2)

  n_bins <- 180
  bin_edges <- seq(-pi, pi, length.out = n_bins + 1)
  poly_row <- poly_col <- numeric(n_bins)
  valid <- logical(n_bins)
  for (i in 1:n_bins) {
    in_bin <- angles >= bin_edges[i] & angles < bin_edges[i + 1]
    if (any(in_bin)) {
      max_idx <- which(in_bin)[which.max(dists[in_bin])]
      poly_col[i] <- blue_pts$col[max_idx]
      poly_row[i] <- blue_pts$row[max_idx]
      valid[i] <- TRUE
    }
  }
  poly_col <- poly_col[valid]; poly_row <- poly_row[valid]

  # Smooth (circular k=3)
  n <- length(poly_col)
  if (n > 10) {
    k <- 3
    pc_s <- pr_s <- numeric(n)
    for (i in 1:n) {
      idx <- ((i + (-k:k)) - 1) %% n + 1
      pc_s[i] <- mean(poly_col[idx]); pr_s[i] <- mean(poly_row[idx])
    }
    poly_col <- pc_s; poly_row <- pr_s
  }

  list(col = poly_col, row = poly_row)
}

# ============================================================================
# Load image and process
# ============================================================================

img_path <- "images/main_ERT/876_18vii24.jpg"
img <- readJPEG(img_path)
H <- nrow(img); W <- ncol(img)

# --- Colorbar extraction & calibration ---
cb <- extract_colorbar(img)
calib_pos <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
calib_val <- c(30, 61, 125, 254, 518, 1000)
calib_func <- create_log_calibration(calib_pos, calib_val)

# --- Blue-polygon boundary detection ---
poly <- auto_detect_polygon_jpeg(img)

# Create mask using point-in-polygon
# Note: expand.grid row order must match matrix fill order (row varies fastest)
grid_all <- expand.grid(row = 1:H, col = 1:W)
pip <- point.in.polygon(grid_all$col, grid_all$row, poly$col, poly$row)
mask_mat <- matrix(FALSE, nrow = H, ncol = W)
for (i in which(pip > 0)) {
  mask_mat[grid_all$row[i], grid_all$col[i]] <- TRUE
}

# Erode 3 pixels
for (e in 1:3) {
  er <- mask_mat
  er[1:(H-1), ] <- er[1:(H-1), ] & mask_mat[2:H, ]
  er[2:H, ]     <- er[2:H, ]     & mask_mat[1:(H-1), ]
  er[, 1:(W-1)] <- er[, 1:(W-1)] & mask_mat[, 2:W]
  er[, 2:W]     <- er[, 2:W]     & mask_mat[, 1:(W-1)]
  mask_mat <- er
}

# Extract masked pixel coordinates and colors
masked_idx <- which(mask_mat, arr.ind = TRUE)  # row, col
colnames(masked_idx) <- c("row", "col")
px_r <- img[cbind(masked_idx[,"row"], masked_idx[,"col"], 1)]
px_g <- img[cbind(masked_idx[,"row"], masked_idx[,"col"], 2)]
px_b <- img[cbind(masked_idx[,"row"], masked_idx[,"col"], 3)]

# Map to resistivity
cat("Mapping", length(px_r), "pixels to resistivity values...\n")
resistivity <- map_colors_to_values(px_r, px_g, px_b, cb, calib_func)

# Spatial data (flip y for plotting)
spatial <- data.frame(
  x = masked_idx[,"col"],
  y = H - masked_idx[,"row"] + 1,
  resistivity = resistivity
)

# --- Compute radial profile ---
cx <- mean(spatial$x); cy <- mean(spatial$y)
spatial$r <- sqrt((spatial$x - cx)^2 + (spatial$y - cy)^2)
spatial$r_norm <- spatial$r / max(spatial$r)

n_rings <- 8
ring_breaks <- seq(0, 1, length.out = n_rings + 1)
spatial$ring <- cut(spatial$r_norm, breaks = ring_breaks,
                    include.lowest = TRUE, labels = FALSE)

ring_profile <- spatial %>%
  group_by(ring) %>%
  summarise(
    mean_res = mean(resistivity, na.rm = TRUE),
    se_res   = sd(resistivity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    r_mid = (ring_breaks[-length(ring_breaks)] + ring_breaks[-1]) / 2
  )

# --- Summary metrics (from pre-computed results) ---
metrics_df <- data.frame(
  Metric = c("Mean (Ohm-m)", "CV", "Gini", "CMA", "Radial grad."),
  Value  = c("505.8", "0.450", "0.254", "0.483", "+504.7"),
  stringsAsFactors = FALSE
)

# ============================================================================
# Panel A -- Raw ERT tomogram
# ============================================================================

panel_a <- ggplot() +
  annotation_raster(img, xmin = 0, xmax = W, ymin = 0, ymax = H) +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  labs(title = "Raw ERT tomogram") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold",
                              margin = margin(b = 2)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 4, 2, 4)
  )

# ============================================================================
# Panel B -- Boundary delineation & calibration (two sub-panels stacked)
# ============================================================================

# ERT color scale
ert_colors <- c("#00008B", "#0000FF", "#00BFFF", "#00FFFF", "#00FF00",
                "#ADFF2F", "#FFFF00", "#FFA500", "#FF4500", "#FF0000", "#8B0000")

# B-top: calibrated cross-section with polygon boundary overlay
panel_b_top <- ggplot(spatial, aes(x = x, y = y, fill = resistivity)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradientn(
    colours = ert_colors,
    trans = "log10",
    limits = c(30, 1000),
    breaks = c(30, 100, 300, 1000),
    labels = c("30", "100", "300", "1000"),
    name = expression(rho ~ (Omega %.% m)),
    oob = scales::squish
  ) +
  coord_fixed(expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.margin = margin(t = 0, b = 0),
    plot.margin = margin(0, 2, 0, 2)
  )

# B-bottom: calibration curve (colorbar position vs resistivity)
calib_curve_df <- data.frame(
  position = seq(0, 1, length.out = 200)
)
calib_curve_df$resistivity <- calib_func(calib_curve_df$position)
calib_points <- data.frame(position = calib_pos, resistivity = calib_val)

panel_b_bottom <- ggplot(calib_curve_df, aes(x = position, y = resistivity)) +
  geom_line(color = "grey30", linewidth = 0.6) +
  geom_point(data = calib_points, aes(x = position, y = resistivity),
             color = "red", size = 1.5) +
  scale_y_log10(
    breaks = c(30, 100, 300, 1000),
    labels = c("30", "100", "300", "1000")
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  labs(
    x = "Colorbar position",
    y = expression(rho ~ (Omega %.% m))
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.title = element_text(size = 6),
    axis.text  = element_text(size = 5),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 6, 2, 6)
  )

# Stack B-top + B-bottom
panel_b <- plot_grid(
  panel_b_top,
  panel_b_bottom,
  ncol = 1, rel_heights = c(0.65, 0.45)
)

# Add title
panel_b_titled <- plot_grid(
  ggdraw() + draw_label("Delineation & calibration",
                         fontface = "bold", size = 9, x = 0.55),
  panel_b,
  ncol = 1, rel_heights = c(0.06, 1)
)

# ============================================================================
# Panel C -- Quantitative output: radial profile + metrics table
# ============================================================================

profile_plot <- ggplot(ring_profile, aes(x = r_mid, y = mean_res)) +
  geom_ribbon(aes(ymin = mean_res - se_res, ymax = mean_res + se_res),
              fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 2) +
  scale_x_continuous(
    breaks = ring_profile$r_mid,
    labels = 1:8,
    expand = expansion(mult = 0.05)
  ) +
  labs(
    x = "Ring (center to edge)",
    y = expression(Mean ~ rho ~ (Omega %.% m))
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 7),
    axis.text  = element_text(size = 6),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 4, 2, 4)
  )

tbl_theme <- ttheme_minimal(
  core    = list(fg_params = list(cex = 0.55, col = "grey20"),
                 bg_params = list(fill = c("grey95", "white")),
                 padding = unit(c(4, 3), "mm")),
  colhead = list(fg_params = list(cex = 0.6, fontface = "bold"),
                 bg_params = list(fill = "grey85"),
                 padding = unit(c(4, 3), "mm"))
)
table_grob <- tableGrob(metrics_df, rows = NULL, theme = tbl_theme)

panel_c <- plot_grid(
  profile_plot,
  ggdraw() + draw_grob(table_grob),
  ncol = 1, rel_heights = c(1, 0.8)
)

panel_c_titled <- plot_grid(
  ggdraw() + draw_label("Quantitative output", fontface = "bold", size = 10),
  panel_c,
  ncol = 1, rel_heights = c(0.06, 1)
)

# ============================================================================
# Assemble 3-panel figure with arrows
# ============================================================================

arrow_grob <- ggdraw() +
  draw_label(expression(symbol("\256")), size = 22, color = "grey40")

fig <- plot_grid(
  panel_a, arrow_grob, panel_b_titled, arrow_grob, panel_c_titled,
  nrow = 1,
  rel_widths = c(1, 0.1, 1, 0.1, 1),
  labels = c("A", "", "B", "", "C"),
  label_size = 12,
  label_fontface = "bold"
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Save
ggsave("output/figures/fig_ert_workflow.pdf", fig,
       width = 7.5, height = 3.8, units = "in")
ggsave("output/figures/fig_ert_workflow.png", fig,
       width = 7.5, height = 3.8, units = "in", dpi = 300, bg = "white")

cat("Figure saved to output/figures/fig_ert_workflow.pdf and .png\n")
