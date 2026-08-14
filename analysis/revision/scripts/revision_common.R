# Shared helpers for the revision sensitivity analyses (cjfr-2026-0202).
# Source from any script in analysis/revision/scripts/; run from the repo root.
#
# Reconstructs the manuscript's ERT PC1 (species-normalized z-scores, PCA on the
# 57 main trees, oriented so higher = wetter / lower resistivity) and reproduces
# the published values (PC1 47.1% / PC2 30.0% variance; classification 30/12/12/3),
# consistent with code/final_phase_and_scans.R.

SPECIES_MAP <- c(rm = "A.rubrum", bg = "N.sylvatica", ro = "Q.rubra", hem = "T.canadensis")
ERT_METRICS <- c("mean", "median", "sd", "cv", "gini", "entropy", "cma", "radialgradient")

SPECIES_COLS <- c(A.rubrum = "#2c7fb8", T.canadensis = "#41b6c4",
                  N.sylvatica = "#7fcdbb", Q.rubra = "#d95f0e")
CAT_COLS <- c(I = "#3b6fb0", II = "#4e9a2c", III = "#d98a1f", IV = "#b83232")

load_main <- function() {
  df <- read.csv("data/Tree_ID_info.csv")
  df$site <- ifelse(df$plot == "BGS", "BGS", "EMS")
  df$sp <- unname(SPECIES_MAP[df$species])
  df
}

# Main 57 trees merged with ERT metrics (inner join, left order preserved)
load_merged <- function() {
  df <- load_main()
  ert <- read.csv("data/ERT_application_results.csv")
  m <- merge(df, ert, by = "tree", sort = FALSE)
  m <- m[order(match(m$tree, df$tree)), ]
  rownames(m) <- NULL
  m
}

# Hemlock validation set; mean resistivity from conductance
load_validation <- function() {
  v <- read.csv("data/hemlock/validation_summary.csv")
  v$mean_res <- 1000 / v$Conductance
  v
}

species_z <- function(m, cols = ERT_METRICS) {
  Z <- m[cols]
  for (cc in cols) {
    Z[[cc]] <- ave(m[[cc]], m$sp, FUN = function(x) (x - mean(x)) / sd(x))
  }
  Z
}

# PCA on z-scored metrics via SVD; returns scores, variance fractions, loadings.
# scheme = "species" (within-species z, published) or "pooled".
# PC1 oriented so that higher = lower mean resistivity (wetter).
build_pc1 <- function(m, scheme = c("species", "pooled"), cols = ERT_METRICS) {
  scheme <- match.arg(scheme)
  Z <- if (scheme == "species") species_z(m, cols) else as.data.frame(scale(m[cols]))
  Zc <- scale(as.matrix(Z), center = TRUE, scale = FALSE)
  sv <- svd(Zc)
  scores <- sv$u %*% diag(sv$d)
  ve <- sv$d^2 / sum(sv$d^2)
  pc1 <- scores[, 1]
  load1 <- sv$v[, 1]
  if (cor(pc1, m$mean) > 0) {
    pc1 <- -pc1
    load1 <- -load1
  }
  names(load1) <- cols
  list(pc1 = pc1, ve = ve, loadings = load1, scores = scores, Z = Z)
}

# Species-standardized resistivity deviation, oriented so higher = wetter
xdev_species_median <- function(m) {
  -ave(m$mean, m$sp, FUN = function(x) (x - median(x)) / sd(x))
}

# Four-quadrant classification. anom, struct: logical vectors.
four_cat <- function(anom, struct) {
  ifelse(!struct & !anom, "I",
    ifelse(!struct & anom, "II",
      ifelse(struct & anom, "III", "IV")))
}

# Otsu threshold on a numeric vector (histogram-based, as in the Python scripts)
otsu_threshold <- function(x, nbins = 256) {
  h <- hist(x, breaks = nbins, plot = FALSE)
  # replicate numpy.histogram: exactly nbins equal-width bins over the range
  edges <- seq(min(x), max(x), length.out = nbins + 1)
  counts <- as.numeric(table(cut(x, edges, include.lowest = TRUE)))
  centers <- (edges[-1] + edges[-length(edges)]) / 2
  w <- counts / sum(counts)
  cw <- cumsum(w); cm <- cumsum(w * centers); gm <- cm[length(cm)]
  best_t <- NA_real_; best_v <- -1
  for (i in seq_len(nbins - 1) + 1) {   # i = 2..nbins, matching python range(1, nbins)
    w0 <- cw[i]; w1 <- 1 - w0
    if (w0 <= 0 || w1 <= 0) next
    v <- w0 * w1 * ((cm[i] / w0) - ((gm - cm[i]) / w1))^2
    if (v > best_v) { best_v <- v; best_t <- centers[i] }
  }
  best_t
}

# 1-D 2-component Gaussian mixture via mclust; returns BICs and the density
# crossover between the two weighted component densities (between the means).
gmm_1d <- function(x) {
  suppressMessages(library(mclust, quietly = TRUE))
  x <- as.numeric(x)
  g1 <- mclust::Mclust(x, G = 1, modelNames = "V", verbose = FALSE)
  g2 <- mclust::Mclust(x, G = 2, modelNames = "V", verbose = FALSE)
  # mclust BIC = 2*loglik - k*log(n)  (higher better); convert to the
  # sklearn convention BIC = -2*loglik + k*log(n) (lower better)
  bic1 <- -g1$bic; bic2 <- -g2$bic
  mns <- g2$parameters$mean
  o <- order(mns)
  xs <- seq(mns[o[1]], mns[o[2]], length.out = 2000)
  p <- g2$parameters
  lp1 <- log(p$pro[o[1]]) + dnorm(xs, p$mean[o[1]], sqrt(p$variance$sigmasq[min(o[1], length(p$variance$sigmasq))]), log = TRUE)
  lp2 <- log(p$pro[o[2]]) + dnorm(xs, p$mean[o[2]], sqrt(p$variance$sigmasq[min(o[2], length(p$variance$sigmasq))]), log = TRUE)
  d <- lp1 - lp2
  idx <- which(diff(sign(d)) != 0)
  crossover <- if (length(idx)) xs[idx[1]] else NA_real_
  list(bic1 = bic1, bic2 = bic2, crossover = crossover, fit2 = g2)
}

# Find a tree's tomogram jpg, excluding duplicates/test shots
find_scan <- function(dir, tree_id,
                      exclude = c("labtest", "_2.jpg", "_jon", "30ii", "dummy", "test")) {
  g <- sort(Sys.glob(file.path(dir, paste0(tree_id, "_*.jpg"))))
  keep <- !Reduce(`|`, lapply(exclude, function(s) grepl(s, g, fixed = TRUE)))
  g <- g[keep]
  if (!length(g)) stop("no scan found for tree ", tree_id, " in ", dir)
  g[1]
}

OUT_DIR <- "analysis/revision/output"
