library(tidyverse)

dat <- read_csv("data/Tree_ID_info.csv", show_col_types=FALSE) %>% rename_all(tolower)
ert <- read_csv("data/ERT_application_results.csv", show_col_types=FALSE) %>% rename_all(tolower)
dat <- dat %>%
  left_join(ert %>% group_by(tree) %>%
    summarise(across(c(mean,median,sd,cv,gini,entropy,cma,radialgradient,npixels),
      \(x) mean(x,na.rm=TRUE)), .groups="drop"), by="tree") %>%
  mutate(structural_loss=percent_damaged, dataset="training",
         tree=as.character(tree),
         abs_cma=abs(cma), abs_radgrad=abs(radialgradient),
         neg_mean=-mean, neg_median=-median)

hem_val_ert <- read_csv("data/hemlock/validation_summary.csv", show_col_types=FALSE)
hem_sot <- read_csv("data/hemlock/SOT_results.csv", show_col_types=FALSE)
hem_sot_dbh <- hem_sot %>%
  mutate(tree_id=str_extract(Filename,"HF_[0-9]+")) %>%
  filter(str_detect(Filename,"DBH")) %>%
  group_by(tree_id) %>% summarise(percent_damaged=mean(pct_damaged,na.rm=TRUE),.groups="drop")

hem_val <- hem_val_ert %>%
  rename(mean=Mean,median=Median,sd=SD,cv=CV,gini=Gini,entropy=Entropy,
         cma=CMA,radialgradient=RadialGradient) %>%
  left_join(hem_sot_dbh, by="tree_id") %>%
  mutate(tree=tree_id, species="hem", site="HF", dataset="validation",
         structural_loss=percent_damaged,
         abs_cma=abs(cma), abs_radgrad=abs(radialgradient),
         neg_mean=-mean, neg_median=-median)

dat <- bind_rows(dat, hem_val)

znorm_from_training <- function(dat, metric_col, group_vars) {
  train_stats <- dat %>% filter(dataset=="training") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(mu=mean(.data[[metric_col]],na.rm=TRUE),
              sigma=sd(.data[[metric_col]],na.rm=TRUE), .groups="drop")
  dat %>% left_join(train_stats, by=group_vars) %>%
    mutate(z=if_else(!is.na(sigma)&sigma>0, (.data[[metric_col]]-mu)/sigma, .data[[metric_col]]-mu)) %>%
    pull(z)
}

sot_threshold <- 1

# CMA species x site znorm
dat$cma_z <- znorm_from_training(dat, "cma", c("species","site"))
na_mask <- is.na(dat$cma_z)
if(any(na_mask)) dat$cma_z[na_mask] <- znorm_from_training(dat,"cma","species")[na_mask]
ert_thresh_cma <- mean(dat$cma_z[dat$dataset=="training"], na.rm=TRUE)

# neg_median species x site znorm
dat$negmed_z <- znorm_from_training(dat, "neg_median", c("species","site"))
na_mask <- is.na(dat$negmed_z)
if(any(na_mask)) dat$negmed_z[na_mask] <- znorm_from_training(dat,"neg_median","species")[na_mask]
ert_thresh_negmed <- mean(dat$negmed_z[dat$dataset=="training"], na.rm=TRUE)

# CV species x site znorm
dat$cv_z <- znorm_from_training(dat, "cv", c("species","site"))
na_mask <- is.na(dat$cv_z)
if(any(na_mask)) dat$cv_z[na_mask] <- znorm_from_training(dat,"cv","species")[na_mask]
find_otsu <- function(values) {
  values <- values[!is.na(values)]
  best <- -Inf; best_t <- median(values)
  for (t in quantile(values, seq(0.1, 0.9, 0.02))) {
    w0 <- mean(values <= t); w1 <- 1 - w0
    if (w0 > 0 & w1 > 0) {
      bv <- w0 * w1 * (mean(values[values <= t]) - mean(values[values > t]))^2
      if (bv > best) { best <- bv; best_t <- t }
    }
  }
  best_t
}
ert_thresh_cv <- find_otsu(dat$cv_z[dat$dataset=="training"])

# PCA
pca_metrics <- c("mean","median","sd","cv","gini","entropy","abs_cma","abs_radgrad")
train_spp_stats <- dat %>% filter(dataset=="training") %>% group_by(species) %>%
  summarise(across(all_of(pca_metrics), list(mu=~mean(.,na.rm=TRUE), sigma=~sd(.,na.rm=TRUE))), .groups="drop")
pca_normed <- dat %>% select(tree, species, dataset, all_of(pca_metrics)) %>%
  left_join(train_spp_stats, by="species")
for(m in pca_metrics) {
  pca_normed[[m]] <- (pca_normed[[m]] - pca_normed[[paste0(m,"_mu")]]) / pca_normed[[paste0(m,"_sigma")]]
}
pca_mat <- pca_normed %>% select(all_of(pca_metrics)) %>% as.matrix()
pca_mat[is.nan(pca_mat)] <- 0
train_rows <- which(dat$dataset=="training")
pca_fit <- prcomp(pca_mat[train_rows,], center=FALSE, scale.=FALSE)
all_scores <- pca_mat %*% pca_fit$rotation
dat$pc1 <- all_scores[,1]; dat$pc2 <- all_scores[,2]
if(pca_fit$rotation["mean",1] > 0) dat$pc1 <- -dat$pc1
ert_thresh_pc1 <- mean(dat$pc1[dat$dataset=="training"], na.rm=TRUE)
ert_thresh_pc2 <- mean(dat$pc2[dat$dataset=="training"], na.rm=TRUE)

assign_quad <- function(sot, ert, sot_t, ert_t) {
  case_when(sot<=sot_t & ert<=ert_t ~ "I",
            sot<=sot_t & ert>ert_t ~ "II",
            sot>sot_t & ert>ert_t ~ "III",
            sot>sot_t & ert<=ert_t ~ "IV")
}

dat$q_cma <- assign_quad(dat$structural_loss, dat$cma_z, sot_threshold, ert_thresh_cma)
dat$q_negmed <- assign_quad(dat$structural_loss, dat$negmed_z, sot_threshold, ert_thresh_negmed)
dat$q_cv <- assign_quad(dat$structural_loss, dat$cv_z, sot_threshold, ert_thresh_cv)
dat$q_pc1 <- assign_quad(dat$structural_loss, dat$pc1, sot_threshold, ert_thresh_pc1)
dat$q_pc2 <- assign_quad(dat$structural_loss, dat$pc2, sot_threshold, ert_thresh_pc2)

cat("\n=== QUADRANT COUNTS (training) ===\n")
axes <- c("q_cma","q_negmed","q_cv","q_pc1","q_pc2")
labels <- c("CMA spp×site","neg_median spp×site","CV spp×site","PC1 (spp)","PC2 (spp)")
for(k in seq_along(axes)) {
  tb <- table(dat[[axes[k]]][dat$dataset=="training"])
  cat(sprintf("%-22s  I=%2d  II=%2d  III=%2d  IV=%2d\n", labels[k],
    tb["I"], tb["II"], tb["III"], tb["IV"]))
}

cat("\n=== PAIRWISE AGREEMENT (training, n=57) ===\n")
train <- dat %>% filter(dataset=="training")
for(i in 1:(length(axes)-1)) {
  for(j in (i+1):length(axes)) {
    agree <- mean(train[[axes[i]]] == train[[axes[j]]], na.rm=TRUE)
    cat(sprintf("%-22s vs %-22s : %4.1f%%\n", labels[i], labels[j], agree*100))
  }
}

cat("\n=== CROSS-TAB: CMA vs PC1 (training) ===\n")
print(table(CMA=train$q_cma, PC1=train$q_pc1))

cat("\n=== CROSS-TAB: CMA vs neg_median (training) ===\n")
print(table(CMA=train$q_cma, negMed=train$q_negmed))

cat("\n=== CROSS-TAB: PC1 vs PC2 (training) ===\n")
print(table(PC1=train$q_pc1, PC2=train$q_pc2))

cat("\n=== VALIDATION HEMLOCK ASSIGNMENTS ===\n")
val <- dat %>% filter(dataset=="validation") %>%
  select(tree, structural_loss, q_cma, q_negmed, q_cv, q_pc1, q_pc2)
print(as.data.frame(val), row.names=FALSE)

cat("\n=== VALIDATION AGREEMENT ===\n")
val_dat <- dat %>% filter(dataset=="validation")
for(i in 1:(length(axes)-1)) {
  for(j in (i+1):length(axes)) {
    agree <- mean(val_dat[[axes[i]]] == val_dat[[axes[j]]], na.rm=TRUE)
    cat(sprintf("%-22s vs %-22s : %4.1f%%\n", labels[i], labels[j], agree*100))
  }
}

cat("\n=== PC LOADINGS (reminder) ===\n")
print(round(pca_fit$rotation[,1:3], 3))

cat("\n=== THRESHOLDS USED ===\n")
cat("SoT:", sot_threshold, "\n")
cat("CMA z:", round(ert_thresh_cma,3), "\n")
cat("neg_median z:", round(ert_thresh_negmed,3), "\n")
cat("CV z (Otsu):", round(ert_thresh_cv,3), "\n")
cat("PC1:", round(ert_thresh_pc1,3), "\n")
cat("PC2:", round(ert_thresh_pc2,3), "\n")
