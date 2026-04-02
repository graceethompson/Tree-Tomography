# ERT Analysis Shiny App - Non-Linear Colorbar Mapping
# Save as: ert_app_nonlinear.R

library(shiny)
library(shinydashboard)
library(DT)
library(imager)
library(dplyr)
library(DescTools)
library(sp)

# ===== Processing Functions =====

# Helper: cumulative arc-length parameter (robust to uneven sampling)
arc_length_param <- function(col_df) {
  if (nrow(col_df) < 2) return(rep(0, nrow(col_df)))
  d <- sqrt(diff(col_df$r)^2 + diff(col_df$g)^2 + diff(col_df$b)^2)
  s <- c(0, cumsum(d))
  if (max(s) > 0) s / max(s) else s  # normalize 0..1
}

# Extract rows 3:8 only
extract_colorbar <- function(img) {
  W <- width(img)
  rows <- 3:8
  
  # Average rows 3-8 across RGB
  r_mat <- matrix(NA_real_, nrow = length(rows), ncol = W)
  g_mat <- matrix(NA_real_, nrow = length(rows), ncol = W)
  b_mat <- matrix(NA_real_, nrow = length(rows), ncol = W)
  
  for (i in seq_along(rows)) {
    r_mat[i, ] <- as.vector(img[, rows[i], 1, 1])
    g_mat[i, ] <- as.vector(img[, rows[i], 1, 2])
    b_mat[i, ] <- as.vector(img[, rows[i], 1, 3])
  }
  
  r <- colMeans(r_mat, na.rm = TRUE)
  g <- colMeans(g_mat, na.rm = TRUE)
  b <- colMeans(b_mat, na.rm = TRUE)
  
  cb <- data.frame(r = r, g = g, b = b)
  
  # Optional: trim a few pixels at each edge to avoid border artifacts
  if (nrow(cb) > 6) {
    cb <- cb[3:(nrow(cb)-2), ]
  }
  
  # Light smoothing to denoise but preserve ends
  if (nrow(cb) >= 5) {
    k <- 5
    smooth <- function(v) {
      n <- length(v)
      result <- numeric(n)
      for(i in 1:n) {
        start_idx <- max(1, i - 2)
        end_idx <- min(n, i + 2)
        result[i] <- mean(v[start_idx:end_idx])
      }
      result
    }
    cb$r <- smooth(cb$r)
    cb$g <- smooth(cb$g)
    cb$b <- smooth(cb$b)
  }
  
  # Downsample uniformly to at most 256 samples
  if (nrow(cb) > 256) {
    idx <- round(seq(1, nrow(cb), length.out = 256))
    cb <- cb[idx, , drop = FALSE]
  }
  
  # Auto direction: expect blue-ish at low end, red-ish at high end
  left_blue_score  <- cb$b[1] - cb$r[1]
  right_red_score  <- cb$r[nrow(cb)] - cb$b[nrow(cb)]
  if ( (left_blue_score + right_red_score) < 0 ) {
    cb <- cb[nrow(cb):1, , drop = FALSE]
  }
  
  # Store arc-length parameter
  attr(cb, "t") <- arc_length_param(cb)
  cb
}

# Create calibration function from user-specified points
create_calibration_function <- function(calib_positions, calib_values, transform = "log") {
  # Remove any NA entries
  valid <- !is.na(calib_positions) & !is.na(calib_values)
  if(sum(valid) < 2) {
    stop("Need at least 2 calibration points")
  }
  
  positions <- calib_positions[valid]
  values <- calib_values[valid]
  
  # Sort by position
  ord <- order(positions)
  positions <- positions[ord]
  values <- values[ord]
  
  # Create interpolation function based on transform
  if(transform == "log") {
    # Log transform for typical resistivity data
    log_values <- log10(values)
    
    # Fit polynomial or spline to log-transformed data
    if(length(positions) <= 3) {
      # Simple linear interpolation in log space
      interp_func <- approxfun(positions, log_values, rule = 2)
      calib_func <- function(x) {
        10^interp_func(x)
      }
    } else {
      # Smooth spline in log space
      fit <- smooth.spline(positions, log_values, spar = 0.5)
      calib_func <- function(x) {
        10^predict(fit, x)$y
      }
    }
  } else if(transform == "linear") {
    # Linear interpolation
    calib_func <- approxfun(positions, values, rule = 2)
  } else if(transform == "power") {
    # Power law fit: value = a * position^b
    if(length(positions) >= 2) {
      log_pos <- log10(positions[positions > 0])
      log_val <- log10(values[positions > 0])
      fit <- lm(log_val ~ log_pos)
      a <- 10^coef(fit)[1]
      b <- coef(fit)[2]
      calib_func <- function(x) {
        a * x^b
      }
    } else {
      calib_func <- approxfun(positions, values, rule = 2)
    }
  }
  
  return(calib_func)
}

# Map pixel colors using calibrated scale
map_colors_calibrated <- function(pixel_colors, colorbar_colors, calib_func) {
  if (is.null(colorbar_colors) || nrow(colorbar_colors) < 2) {
    stop("Colorbar not available or too short")
  }
  
  cb_mat <- as.matrix(colorbar_colors[, c("r","g","b")])
  t_cb   <- attr(colorbar_colors, "t")
  if (is.null(t_cb)) {
    t_cb <- seq(0, 1, length.out = nrow(colorbar_colors))
  }
  
  px <- as.matrix(pixel_colors)
  n  <- nrow(px)
  t_hat <- numeric(n)
  
  # Chunk to keep memory manageable
  chunk <- 2000L
  for (start in seq(1L, n, by = chunk)) {
    end <- min(n, start + chunk - 1L)
    x <- px[start:end, , drop = FALSE]
    
    # Euclidean NN in RGB
    for (i in 1:nrow(x)) {
      d <- (cb_mat[,1] - x[i,1])^2 + (cb_mat[,2] - x[i,2])^2 + (cb_mat[,3] - x[i,3])^2
      k <- which.min(d)
      t_hat[start + i - 1L] <- t_cb[k]
    }
  }
  
  # Apply calibration function
  vals <- calib_func(t_hat)
  vals
}

process_ert_with_mask <- function(img, mask_vector, image_name, 
                                  calib_positions = NULL, calib_values = NULL,
                                  transform = "log") {
  
  W <- width(img)
  H <- height(img)
  
  # Extract horizontal colorbar from top of image
  colorbar <- extract_colorbar(img)
  
  # Extract values from masked region
  img_df <- as.data.frame(img, wide = "c") %>%
    mutate(
      inside = mask_vector,
      r = c.1,
      g = c.2, 
      b = c.3
    ) %>%
    filter(inside == TRUE)
  
  if(nrow(img_df) == 0) {
    stop("No pixels found in mask")
  }
  
  # Get pixel colors
  pixel_colors <- as.matrix(img_df[, c("r", "g", "b")])
  
  # Map colors to values
  if(!is.null(colorbar) && nrow(colorbar) > 10) {
    if(!is.null(calib_positions) && !is.null(calib_values)) {
      # Use calibration curve
      calib_func <- create_calibration_function(calib_positions, calib_values, transform)
      values <- map_colors_calibrated(pixel_colors, colorbar, calib_func)
    } else {
      # Fallback to linear 0-1
      warning("No calibration provided, using 0-1 scale")
      cb_mat <- as.matrix(colorbar[, c("r","g","b")])
      t_cb <- attr(colorbar, "t")
      if(is.null(t_cb)) t_cb <- seq(0, 1, length.out = nrow(colorbar))
      
      px <- as.matrix(pixel_colors)
      values <- numeric(nrow(px))
      for(i in 1:nrow(px)) {
        d <- (cb_mat[,1] - px[i,1])^2 + (cb_mat[,2] - px[i,2])^2 + (cb_mat[,3] - px[i,3])^2
        k <- which.min(d)
        values[i] <- t_cb[k]
      }
    }
  } else {
    stop("Could not extract colorbar")
  }
  
  # Remove NA values
  valid_idx <- !is.na(values)
  values <- values[valid_idx]
  img_df <- img_df[valid_idx,]
  
  if(length(values) == 0) {
    stop("No valid values after processing")
  }
  
  # Calculate metrics
  cx <- mean(img_df$x)
  cy <- mean(img_df$y)
  r <- sqrt((img_df$x - cx)^2 + (img_df$y - cy)^2)
  r_norm <- r / max(r)
  
  # 8-ring radial profile
  n_rings <- 8
  ring_breaks <- seq(0, 1, length.out = n_rings + 1)
  ring_ids <- cut(r_norm, breaks = ring_breaks, include.lowest = TRUE, labels = FALSE)
  ring_means <- tapply(values, ring_ids, mean, na.rm = TRUE)
  ring_means[is.na(ring_means)] <- 0
  
  # CMA (Central Moisture Accumulation)
  q_low <- quantile(values, probs = 0.30, na.rm = TRUE)
  is_low <- values <= q_low
  is_inner <- r_norm <= 0.33
  
  if(sum(is_low) > 0) {
    CMA <- sum(is_low & is_inner) / sum(is_low)
  } else {
    CMA <- 0
  }
  
  # Center and edge means
  center_vals <- values[r_norm <= 0.33]
  edge_vals <- values[r_norm >= 0.67]
  
  center_mean <- if(length(center_vals) > 0) mean(center_vals) else mean(values)
  edge_mean <- if(length(edge_vals) > 0) mean(edge_vals) else mean(values)
  
  # Return comprehensive results
  list(
    file = image_name,
    mean = mean(values),
    median = median(values),
    sd = sd(values),
    cv = sd(values) / mean(values),
    gini = Gini(values),
    entropy = Entropy(values, method = "ML"),
    CMA = CMA,
    center_mean = center_mean,
    edge_mean = edge_mean,
    radial_gradient = edge_mean - center_mean,
    n_pixels = length(values),
    ring_profile = as.numeric(ring_means),
    values = values,
    spatial_data = img_df
  )
}

# Auto-detect PiCUS blue boundary polygon from image
auto_detect_polygon <- function(img) {
  W <- width(img)
  H <- height(img)
  colorbar_cutoff <- min(40, round(H * 0.08))

  # Extract pixel RGB below colorbar
  rows <- (colorbar_cutoff + 1):H
  r_vals <- as.vector(img[, rows, 1, 1])
  g_vals <- as.vector(img[, rows, 1, 2])
  b_vals <- as.vector(img[, rows, 1, 3])

  # Grid coordinates (imager: x across columns, y across rows)
  grid <- expand.grid(x = 1:W, row_idx = seq_along(rows))
  grid$y <- rows[grid$row_idx]
  grid$r <- r_vals
  grid$g <- g_vals
  grid$b <- b_vals

  # Detect blue boundary pixels: high blue, low red, low green
  grid$is_blue <- (grid$b > 0.4) & (grid$r < 0.3) & (grid$g < 0.4) &
                  ((grid$b - grid$g) > 0.15)

  blue_pts <- grid[grid$is_blue, c("x", "y")]

  if (nrow(blue_pts) < 20) {
    stop("Could not detect blue boundary polygon (only ", nrow(blue_pts),
         " blue pixels found). Use manual polygon instead.")
  }

  # Centroid of blue pixels
  cx <- median(blue_pts$x)
  cy <- median(blue_pts$y)

  # Bin by angle, take outermost blue pixel per bin
  angles <- atan2(blue_pts$y - cy, blue_pts$x - cx)
  dists  <- sqrt((blue_pts$x - cx)^2 + (blue_pts$y - cy)^2)

  n_bins <- 180
  bin_edges <- seq(-pi, pi, length.out = n_bins + 1)
  poly_x <- numeric(n_bins)
  poly_y <- numeric(n_bins)
  valid   <- logical(n_bins)

  for (i in 1:n_bins) {
    in_bin <- angles >= bin_edges[i] & angles < bin_edges[i + 1]
    if (any(in_bin)) {
      max_idx <- which(in_bin)[which.max(dists[in_bin])]
      poly_x[i] <- blue_pts$x[max_idx]
      poly_y[i] <- blue_pts$y[max_idx]
      valid[i]  <- TRUE
    }
  }

  poly_x <- poly_x[valid]
  poly_y <- poly_y[valid]

  # Smooth (circular moving average, k=3)
  n <- length(poly_x)
  if (n > 10) {
    k <- 3
    px_s <- numeric(n); py_s <- numeric(n)
    for (i in 1:n) {
      idx <- ((i + (-k:k)) - 1) %% n + 1
      px_s[i] <- mean(poly_x[idx])
      py_s[i] <- mean(poly_y[idx])
    }
    poly_x <- px_s; poly_y <- py_s
  }

  # Return as Nx2 matrix (same format as manual polygon_points)
  cbind(poly_x, poly_y)
}

# ===== Shiny UI =====
ui <- dashboardPage(
  dashboardHeader(title = "ERT Analysis - Non-Linear Scale"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Process Images", tabName = "process", icon = icon("image")),
      menuItem("Results Table", tabName = "results", icon = icon("table")),
      menuItem("Help & Guide", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #f4f4f4; }
        .calib-input { margin-bottom: 10px; }
      "))
    ),
    
    tabItems(
      # Main Processing Tab
      tabItem(tabName = "process",
              fluidRow(
                # Left Panel - Controls (narrow)
                column(width = 3,
                       box(
                         title = "Step 1: Upload Images",
                         status = "primary",
                         solidHeader = TRUE,
                         width = NULL,

                         fileInput("images", "Select ERT Images:",
                                   multiple = TRUE,
                                   accept = c(".jpg", ".jpeg", ".png")),

                         verbatimTextOutput("file_info")
                       ),

                       box(
                         title = "Step 2: Colorbar Calibration",
                         status = "warning",
                         solidHeader = TRUE,
                         width = NULL,
                         collapsible = TRUE,
                         collapsed = TRUE,
                         
                         p("Enter scale values at specific positions (0=left, 1=right):"),
                         
                         # Calibration points
                         div(class = "calib-input",
                             fluidRow(
                               column(6, numericInput("pos1", "Position 1:", value = 0, min = 0, max = 1, step = 0.01)),
                               column(6, numericInput("val1", "Value 1:", value = 30))
                             )
                         ),
                         div(class = "calib-input",
                             fluidRow(
                               column(6, numericInput("pos2", "Position 2:", value = 0.2, min = 0, max = 1, step = 0.01)),
                               column(6, numericInput("val2", "Value 2:", value = 61))
                             )
                         ),
                         div(class = "calib-input",
                             fluidRow(
                               column(6, numericInput("pos3", "Position 3:", value = 0.4, min = 0, max = 1, step = 0.01)),
                               column(6, numericInput("val3", "Value 3:", value = 125))
                             )
                         ),
                         div(class = "calib-input",
                             fluidRow(
                               column(6, numericInput("pos4", "Position 4:", value = 0.6, min = 0, max = 1, step = 0.01)),
                               column(6, numericInput("val4", "Value 4:", value = 254))
                             )
                         ),
                         div(class = "calib-input",
                             fluidRow(
                               column(6, numericInput("pos5", "Position 5:", value = 0.8, min = 0, max = 1, step = 0.01)),
                               column(6, numericInput("val5", "Value 5:", value = 518))
                             )
                         ),
                         div(class = "calib-input",
                             fluidRow(
                               column(6, numericInput("pos6", "Position 6:", value = 1, min = 0, max = 1, step = 0.01)),
                               column(6, numericInput("val6", "Value 6:", value = 1000))
                             )
                         ),
                         
                         radioButtons("transform", "Transform:",
                                      choices = c("Logarithmic" = "log",
                                                  "Linear" = "linear",
                                                  "Power Law" = "power"),
                                      selected = "log"),
                         
                         actionButton("preview_calibration", "Preview Calibration", width = "100%")
                       ),
                       
                       box(
                         title = "Step 3: Define Boundary",
                         status = "info",
                         solidHeader = TRUE,
                         width = NULL,

                         actionButton("auto_polygon", "Auto-Detect Boundary",
                                      icon = icon("magic"),
                                      style = "width: 100%; margin-bottom: 10px;",
                                      class = "btn-info"),
                         helpText("Detects PiCUS blue boundary polygon automatically."),
                         hr(),
                         p("Or click points on image manually:"),
                         actionButton("clear_polygon", "Clear Points", width = "45%"),
                         actionButton("finish_polygon", "Finish Polygon", width = "45%"),
                         br(), br(),
                         verbatimTextOutput("polygon_status"),
                         sliderInput("erosion", "Edge Erosion (pixels):", 
                                     min = 0, max = 10, value = 2, step = 1)
                       ),
                       
                       box(
                         title = "Step 4: Process",
                         status = "success",
                         solidHeader = TRUE,
                         width = NULL,
                         
                         actionButton("process_btn", "Process This Image", 
                                      style = "width: 100%; font-size: 16px;"),
                         br(), br(),
                         verbatimTextOutput("process_status")
                       )
                ),
                
                # Right Panel - Image Display
                column(width = 9,
                       conditionalPanel(
                         condition = "input.images",
                         box(
                           width = NULL,
                           fluidRow(
                             column(4, actionButton("prev_img", "← Previous")),
                             column(4, h4(textOutput("image_counter"), align = "center")),
                             column(4, actionButton("next_img", "Next →"))
                           )
                         )
                       ),
                       
                       box(
                         title = "Image Display - Click to add polygon points",
                         status = "primary",
                         solidHeader = TRUE,
                         width = NULL,
                         
                         plotOutput("image_display", 
                                    height = "500px",
                                    click = "plot_click")
                       ),
                       
                       fluidRow(
                         box(
                           title = "Statistics",
                           width = 6,
                           tableOutput("current_stats")
                         ),
                         
                         box(
                           title = "Radial Profile",
                           width = 6,
                           plotOutput("radial_plot", height = "250px")
                         )
                       )
                )
              )
      ),
      
      # Results Table Tab
      tabItem(tabName = "results",
              box(
                title = "Results Table",
                width = 12,
                
                downloadButton("download_csv", "Download CSV"),
                actionButton("clear_results", "Clear All"),
                br(), br(),
                DTOutput("results_table")
              )
      ),
      
      # Help Tab
      tabItem(tabName = "help",
              fluidRow(
                column(width = 12,
                       box(
                         title = "ERT Image Analysis - User Guide",
                         width = NULL,
                         status = "info",
                         solidHeader = TRUE,
                         
                         h3("Overview"),
                         p("This application analyzes Electrical Resistivity Tomography (ERT) images of tree cross-sections to extract quantitative metrics about moisture distribution and heterogeneity. The tool handles non-linear colorbars commonly found in ERT data and can automatically detect the blue measurement polygon from PiCUS ERT exports."),
                         hr(),
                         
                         h3("Step-by-Step Instructions"),
                         
                         h4("Step 1: Upload ERT Images"),
                         tags$ol(
                           tags$li("Click 'Browse' to select one or more ERT images (JPG, PNG format)"),
                           tags$li("Images should have a colorbar at the top showing the resistivity scale"),
                           tags$li("Use the Previous/Next buttons to navigate between multiple images")
                         ),
                         
                         h4("Step 2: Calibrate the Colorbar"),
                         tags$ol(
                           tags$li("Look at the colorbar at the top of your image"),
                           tags$li("Note the labeled values (e.g., 30, 61, 125, 254, 518, 1000 Ohm.m)"),
                           tags$li("For each labeled value, estimate its position along the colorbar:"),
                           tags$ul(
                             tags$li("Position 0 = leftmost edge of colorbar"),
                             tags$li("Position 1 = rightmost edge of colorbar"),
                             tags$li("Position 0.5 = middle of colorbar")
                           ),
                           tags$li("Enter the position-value pairs in the calibration section"),
                           tags$li("Select the appropriate transform:"),
                           tags$ul(
                             tags$li(strong("Logarithmic:"), " Best for resistivity data spanning orders of magnitude (recommended)"),
                             tags$li(strong("Linear:"), " For uniformly spaced scales"),
                             tags$li(strong("Power Law:"), " For power-law relationships")
                           ),
                           tags$li("Click 'Preview Calibration' to verify your calibration")
                         ),
                         
                         h4("Step 3: Define the Tree Boundary"),
                         p(strong("Option A: Auto-Detect (recommended for PiCUS images)")),
                         tags$ol(
                           tags$li("Click 'Auto-Detect Boundary' to automatically detect the blue measurement polygon drawn by the PiCUS software"),
                           tags$li("The app identifies blue boundary pixels, traces the outermost boundary, and creates a smoothed polygon"),
                           tags$li("The polygon status will show the number of detected vertices (typically 150-180 points)"),
                           tags$li("If auto-detection fails (e.g., no blue polygon in the image), use manual mode instead")
                         ),
                         p(strong("Option B: Manual polygon")),
                         tags$ol(
                           tags$li("Click points around the tree cross-section boundary on the image"),
                           tags$li("Points will appear as red dots connected by lines"),
                           tags$li("Click 'Finish Polygon' when you've outlined the complete boundary (minimum 3 points)"),
                           tags$li("Click 'Clear Points' to start over if needed")
                         ),
                         p("Adjust 'Edge Erosion' (default: 2 pixels) to exclude pixels near the boundary and reduce edge artifacts."),
                         
                         h4("Step 4: Process the Image"),
                         tags$ol(
                           tags$li("Click 'Process This Image' to analyze the current image"),
                           tags$li("The Statistics box will display key metrics (Mean, CV, Gini, CMA) for the current image"),
                           tags$li("The Radial Profile plot will show mean resistivity across 8 concentric rings from center to edge"),
                           tags$li("Results are automatically appended to the cumulative Results Table (one row per processed image)"),
                           tags$li("The polygon is cleared after processing so you can proceed to the next image"),
                           tags$li("Use the Previous/Next buttons to navigate between uploaded images and repeat Steps 3-4 for each")
                         ),

                         h4("Step 5: Review and Export Results"),
                         tags$ol(
                           tags$li("Click the 'Results Table' tab in the left sidebar to see all processed images"),
                           tags$li("The table contains one row per image with columns: Filename, Mean, Median, SD, CV, Gini, Entropy, CMA, RadialGradient, NPixels"),
                           tags$li("Click 'Download CSV' to export the full table as a CSV file"),
                           tags$li("Click 'Clear All' to reset the results table (this cannot be undone)")
                         ),
                         p(strong("Important:"), " Results are stored only for the current session. If you close or refresh the browser, all results are lost. Export your CSV before closing the app."),
                         
                         hr(),
                         
                         h3("Understanding the Metrics"),
                         
                         h4("Basic Statistics"),
                         tags$ul(
                           tags$li(strong("Mean:"), " Average resistivity value within the tree boundary. Higher values indicate drier conditions."),
                           tags$li(strong("Median:"), " Middle value when all pixels are sorted. Less sensitive to outliers than mean."),
                           tags$li(strong("SD (Standard Deviation):"), " Measure of variability in resistivity values. Higher SD indicates more heterogeneous moisture distribution.")
                         ),
                         
                         h4("Heterogeneity Metrics"),
                         tags$ul(
                           tags$li(strong("CV (Coefficient of Variation):"), " SD divided by mean. Normalized measure of relative variability (0 = uniform, higher = more variable). Values typically range from 0.1 to 1.0."),
                           tags$li(strong("Gini Coefficient:"), " Measure of inequality in the distribution (0 = perfect equality, 1 = maximum inequality). Higher values indicate concentrated moisture in specific areas."),
                           tags$li(strong("Entropy:"), " Information-theoretic measure of disorder. Higher entropy indicates more random/dispersed moisture patterns.")
                         ),
                         
                         h4("Spatial Metrics"),
                         tags$ul(
                           tags$li(strong("CMA (Central Moisture Accumulation):"), " Proportion of the wettest pixels (lowest 30% resistivity) that are located in the central third of the tree. Values > 0.33 indicate moisture concentration in the center."),
                           tags$li(strong("Radial Gradient:"), " Difference between edge mean and center mean resistivity. Positive values indicate drier edges (healthy pattern), negative values indicate drier center (possible heartwood or decay)."),
                           tags$li(strong("Radial Profile:"), " Mean resistivity values in 8 concentric rings from center to edge. Visualizes the radial moisture distribution pattern.")
                         ),
                         
                         h4("Quality Metrics"),
                         tags$ul(
                           tags$li(strong("NPixels:"), " Number of pixels analyzed within the tree boundary. Indicates image resolution and analysis coverage.")
                         ),
                         
                         hr(),
                         
                         h3("Interpreting Results"),
                         
                         h4("Healthy Trees Typically Show:"),
                         tags$ul(
                           tags$li("Moderate mean resistivity values (depends on species and season)"),
                           tags$li("Lower CV and Gini values (more uniform moisture)"),
                           tags$li("Positive radial gradient (wetter sapwood, drier heartwood)"),
                           tags$li("Smooth radial profile increasing from center to edge")
                         ),
                         
                         h4("Signs of Potential Issues:"),
                         tags$ul(
                           tags$li("Very high CV or Gini (> 0.7): Indicates patchy moisture, possible decay pockets"),
                           tags$li("High CMA (> 0.5): Moisture accumulation in center, possible hollowing"),
                           tags$li("Negative radial gradient: Dry center, could indicate advanced heartwood decay"),
                           tags$li("Irregular radial profile: Non-uniform moisture distribution, structural issues")
                         ),
                         
                         hr(),
                         
                         h3("Troubleshooting"),
                         
                         h4("Colorbar Detection Issues"),
                         tags$ul(
                           tags$li("Ensure the colorbar is horizontal and at the top of the image"),
                           tags$li("The colorbar should be continuous without gaps"),
                           tags$li("Try adjusting calibration positions if mapping seems incorrect")
                         ),
                         
                         h4("Auto-Detect Boundary Issues"),
                         tags$ul(
                           tags$li("'Could not detect blue boundary polygon': The image may not have a PiCUS-style blue measurement polygon. Use manual polygon mode instead."),
                           tags$li("Auto-detected boundary looks wrong: The algorithm looks for blue pixels (high blue channel, low red and green). Non-standard color schemes may confuse it. Switch to manual mode."),
                           tags$li("Boundary includes artifacts: Increase the Edge Erosion slider to shrink the boundary inward by more pixels.")
                         ),

                         h4("Processing Errors"),
                         tags$ul(
                           tags$li("'No pixels in mask': Reduce edge erosion or redraw polygon"),
                           tags$li("'Need at least 2 calibration points': Fill in at least 2 position-value pairs"),
                           tags$li("Unexpected values: Check calibration and ensure colorbar goes from low (blue) to high (red)")
                         ),

                         hr(),

                         h3("About"),
                         p("This application was developed for extracting quantitative metrics from PiCUS ERT tomogram images, enabling reproducible, threshold-based decay classification. It is open-source and available at:"),
                         tags$a(href = "https://github.com/graceethompson/Tree-Tomography",
                                "github.com/graceethompson/Tree-Tomography",
                                target = "_blank")
                       )
                )
              )
      )
    )
  )
)

# ===== Shiny Server =====
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    current_idx = 1,
    polygon_points = NULL,
    polygon_finished = FALSE,
    results_table = data.frame(),
    current_result = NULL
  )
  
  # Current image
  current_image <- reactive({
    req(input$images)
    load.image(input$images$datapath[values$current_idx])
  })
  
  # Get calibration data
  get_calibration <- reactive({
    positions <- c(input$pos1, input$pos2, input$pos3, input$pos4, input$pos5, input$pos6)
    values <- c(input$val1, input$val2, input$val3, input$val4, input$val5, input$val6)
    list(positions = positions, values = values)
  })
  
  # File info
  output$file_info <- renderText({
    if(is.null(input$images)) {
      "No files uploaded"
    } else {
      paste(nrow(input$images), "files uploaded")
    }
  })
  
  # Image display
  output$image_display <- renderPlot({
    if(is.null(input$images)) {
      plot.new()
      text(0.5, 0.5, "Upload images to begin", cex = 2)
      return()
    }
    
    img <- current_image()
    plot(img, axes = FALSE)
    title(input$images$name[values$current_idx])
    
    # Draw polygon if points exist
    if(!is.null(values$polygon_points)) {
      points(values$polygon_points[,1], values$polygon_points[,2], 
             col = "red", pch = 19, cex = 1.5)
      
      if(nrow(values$polygon_points) > 1) {
        lines(values$polygon_points[,1], values$polygon_points[,2], 
              col = "red", lwd = 2)
      }
      
      if(values$polygon_finished && nrow(values$polygon_points) >= 3) {
        polygon(values$polygon_points[,1], values$polygon_points[,2], 
                border = "red", lwd = 2, col = rgb(1,0,0,0.2))
      }
    }
  })
  
  # Handle clicks
  observeEvent(input$plot_click, {
    if(!values$polygon_finished) {
      new_point <- c(input$plot_click$x, input$plot_click$y)
      if(is.null(values$polygon_points)) {
        values$polygon_points <- matrix(new_point, ncol = 2)
      } else {
        values$polygon_points <- rbind(values$polygon_points, new_point)
      }
    }
  })
  
  # Clear polygon
  observeEvent(input$clear_polygon, {
    values$polygon_points <- NULL
    values$polygon_finished <- FALSE
  })
  
  # Finish polygon
  observeEvent(input$finish_polygon, {
    if(!is.null(values$polygon_points) && nrow(values$polygon_points) >= 3) {
      values$polygon_finished <- TRUE
      showNotification("Polygon ready!", type = "message")
    } else {
      showNotification("Need at least 3 points", type = "warning")
    }
  })

  # Auto-detect blue boundary polygon
  observeEvent(input$auto_polygon, {
    req(input$images)
    tryCatch({
      img <- current_image()
      poly <- auto_detect_polygon(img)
      values$polygon_points <- poly
      values$polygon_finished <- TRUE
      showNotification(
        paste("Auto-detected boundary:", nrow(poly), "vertices"),
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Auto-detect failed:", e$message), type = "error")
    })
  })

  # Polygon status
  output$polygon_status <- renderText({
    if(is.null(values$polygon_points)) {
      "No points"
    } else if(values$polygon_finished) {
      paste("Ready:", nrow(values$polygon_points), "points")
    } else {
      paste(nrow(values$polygon_points), "points")
    }
  })
  
  # Preview calibration
  observeEvent(input$preview_calibration, {
    req(input$images)
    
    img <- current_image()
    colorbar <- extract_colorbar(img)
    
    if(!is.null(colorbar) && nrow(colorbar) > 10) {
      calib <- get_calibration()
      
      # Create calibration function
      calib_func <- create_calibration_function(calib$positions, calib$values, input$transform)
      
      # Test the calibration
      test_positions <- seq(0, 1, length.out = 100)
      test_values <- calib_func(test_positions)
      
      showModal(modalDialog(
        title = "Colorbar Calibration Preview",
        size = "l",
        renderPlot({
          par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))
          
          # Show the extracted colorbar
          plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
               xlab = "Position (0=left, 1=right)", ylab = "", 
               main = "Extracted Colorbar", yaxt = "n")
          
          t_cb <- attr(colorbar, "t")
          if(is.null(t_cb)) t_cb <- seq(0, 1, length.out = nrow(colorbar))
          
          for(i in 1:(nrow(colorbar)-1)) {
            rect(t_cb[i], 0, t_cb[i+1], 1, 
                 col = rgb(colorbar$r[i], colorbar$g[i], colorbar$b[i]), 
                 border = NA)
          }
          
          # Add calibration points
          valid <- !is.na(calib$positions) & !is.na(calib$values)
          if(sum(valid) > 0) {
            points(calib$positions[valid], rep(0.5, sum(valid)), 
                   pch = 19, cex = 2, col = "white")
            text(calib$positions[valid], rep(0.5, sum(valid)), 
                 calib$values[valid], col = "black", cex = 0.8)
          }
          
          # Show calibration curve
          plot(test_positions, test_values, type = "l", lwd = 2,
               main = paste("Calibration Curve (", input$transform, ")", sep=""),
               xlab = "Colorbar Position", ylab = "Calibrated Value",
               log = if(input$transform == "log") "y" else "")
          
          points(calib$positions[valid], calib$values[valid], 
                 pch = 19, cex = 1.5, col = "red")
          grid()
          
          # Show RGB components
          plot(t_cb, colorbar$r, type = "l", col = "red", 
               ylim = c(0, 1), xlab = "Position", ylab = "RGB Value",
               main = "RGB Components", lwd = 2)
          lines(t_cb, colorbar$g, col = "green", lwd = 2)
          lines(t_cb, colorbar$b, col = "blue", lwd = 2)
          legend("topright", c("Red", "Green", "Blue"), 
                 col = c("red", "green", "blue"), lty = 1, lwd = 2)
          grid()
          
        }, height = 600),
        p(strong("Calibration ready!")),
        p(paste("Transform type:", input$transform)),
        p("The calibration points are shown as white dots with values."),
        footer = modalButton("Close")
      ))
    } else {
      showNotification("Could not detect colorbar", type = "warning", duration = 5)
    }
  })
  
  # Process image
  observeEvent(input$process_btn, {
    req(input$images)
    
    if(!values$polygon_finished || is.null(values$polygon_points)) {
      showNotification("Draw and finish polygon first", type = "error")
      return()
    }
    
    withProgress(message = 'Processing...', {
      tryCatch({
        img <- current_image()
        W <- width(img)
        H <- height(img)
        
        # Create mask
        grid <- expand.grid(x = 1:W, y = 1:H)
        mask <- point.in.polygon(grid$x, grid$y,
                                 values$polygon_points[,1],
                                 values$polygon_points[,2])
        mask <- mask > 0
        
        # Apply erosion
        if(input$erosion > 0) {
          mask_array <- array(mask, dim = c(W, H, 1, 1))
          mask_cimg <- as.cimg(mask_array)
          
          for(i in 1:input$erosion) {
            mask_cimg <- erode_square(mask_cimg, size = 3)
          }
          
          mask <- as.logical(as.vector(mask_cimg))
        }
        
        if(sum(mask) == 0) {
          stop("No pixels in mask after erosion")
        }
        
        # Get calibration
        calib <- get_calibration()
        
        # Process
        result <- process_ert_with_mask(
          img = img,
          mask_vector = mask,
          image_name = input$images$name[values$current_idx],
          calib_positions = calib$positions,
          calib_values = calib$values,
          transform = input$transform
        )
        
        values$current_result <- result
        
        # Add to table
        new_row <- data.frame(
          Filename = result$file,
          Mean = round(result$mean, 2),
          Median = round(result$median, 2),
          SD = round(result$sd, 2),
          CV = round(result$cv, 4),
          Gini = round(result$gini, 4),
          Entropy = round(result$entropy, 4),
          CMA = round(result$CMA, 4),
          RadialGradient = round(result$radial_gradient, 2),
          NPixels = result$n_pixels,
          stringsAsFactors = FALSE
        )
        
        values$results_table <- rbind(values$results_table, new_row)
        
        showNotification("Processing complete!", type = "message")
        
        # Clear polygon
        values$polygon_points <- NULL
        values$polygon_finished <- FALSE
        
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  
  # Navigation
  observeEvent(input$prev_img, {
    if(values$current_idx > 1) {
      values$current_idx <- values$current_idx - 1
      values$polygon_points <- NULL
      values$polygon_finished <- FALSE
      values$current_result <- NULL
    }
  })
  
  observeEvent(input$next_img, {
    req(input$images)
    if(values$current_idx < nrow(input$images)) {
      values$current_idx <- values$current_idx + 1
      values$polygon_points <- NULL
      values$polygon_finished <- FALSE
      values$current_result <- NULL
    }
  })
  
  # Outputs
  output$image_counter <- renderText({
    req(input$images)
    paste(values$current_idx, "/", nrow(input$images))
  })
  
  output$process_status <- renderText({
    paste(nrow(values$results_table), "processed")
  })
  
  output$current_stats <- renderTable({
    req(values$current_result)
    data.frame(
      Metric = c("Mean", "CV", "Gini", "CMA"),
      Value = c(
        sprintf("%.1f", values$current_result$mean),
        sprintf("%.3f", values$current_result$cv),
        sprintf("%.3f", values$current_result$gini),
        sprintf("%.3f", values$current_result$CMA)
      )
    )
  })
  
  output$radial_plot <- renderPlot({
    req(values$current_result)
    plot(values$current_result$ring_profile,
         type = "b", pch = 19, col = "blue",
         xlab = "Ring", ylab = "Mean Value",
         main = "Radial Profile (8 rings)")
    grid()
  })
  
  output$results_table <- renderDT({
    datatable(values$results_table, options = list(pageLength = 10))
  })
  
  # Download CSV
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("ERT_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$results_table, file, row.names = FALSE)
    }
  )
  
  # Clear results
  observeEvent(input$clear_results, {
    values$results_table <- data.frame()
    values$current_result <- NULL
    showNotification("Results cleared", type = "message")
  })
}

# Run the application
shinyApp(ui = ui, server = server)