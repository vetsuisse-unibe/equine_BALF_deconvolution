############################################################
# mRNA-content correction: leave-one-sample-out (LOSO) sensitivity
# Drop each donor from the reference, recompute m_k from the other 10,
# re-apply the correction, and track m_k and accuracy across folds.
# Resampling is at the donor level, n = 11 (Supplementary Table S6).
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

data_dir    <- "data"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

## ===== Load inputs =====
umi <- read.csv(file.path(data_dir, "per_cell_umi.csv"))
val <- read.csv(file.path(data_dir, "validation_real_bulk_matched.csv"))
umi$Sample <- sub("^A", "", as.character(umi$SampleID))

theta <- val %>% select(Sample, CellType, BayesPrism) %>%
  pivot_wider(names_from = CellType, values_from = BayesPrism) %>%
  column_to_rownames("Sample") %>% as.matrix()
truth <- val %>% select(Sample, CellType, Truth) %>%
  pivot_wider(names_from = CellType, values_from = Truth) %>%
  column_to_rownames("Sample") %>% as.matrix()
samples   <- rownames(theta)
celltypes <- colnames(theta)

correct <- function(theta, m) {
  out <- sweep(theta, 2, m[colnames(theta)], "/")
  out / rowSums(out)
}
score <- function(est) {
  obs <- truth[rownames(est), colnames(est)]
  c(Pearson_r = cor(c(est), c(obs)), MAE = mean(abs(est - obs)))
}

## ===== Full-reference baseline =====
m_full   <- tapply(umi$UMI, umi$CellType, mean)
ref_corr <- correct(theta, m_full)
base     <- score(ref_corr)
message(sprintf("Full reference: r = %.4f, MAE = %.4f", base["Pearson_r"], base["MAE"]))

## ===== LOSO folds (one per donor) =====
mk_list  <- list()
per_fold <- list()
for (s in samples) {
  keep <- umi$Sample != s
  m    <- tapply(umi$UMI[keep], umi$CellType[keep], mean)
  est  <- correct(theta, m)
  sc   <- score(est)
  mk_list[[s]]  <- m[celltypes]
  per_fold[[s]] <- data.frame(HeldOut = s,
                              Pearson_r = sc["Pearson_r"],
                              MAE       = sc["MAE"],
                              Delta_heldout = mean(abs(est[s, ] - ref_corr[s, ])))
}
per_fold <- do.call(rbind, per_fold); rownames(per_fold) <- NULL
mk       <- do.call(rbind, mk_list)

## ===== m_k spread across folds =====
mk_range <- data.frame(
  CellType  = celltypes,
  m_full    = round(m_full[celltypes], 1),
  LOSO_min  = round(apply(mk, 2, min), 1),
  LOSO_max  = round(apply(mk, 2, max), 1),
  range_pct = round(100 * (apply(mk, 2, max) - apply(mk, 2, min)) / m_full[celltypes], 1),
  CV_pct    = round(100 * apply(mk, 2, sd) / colMeans(mk), 1)
)

write.csv(per_fold, file.path(results_dir, "sensitivity_loso_per_fold.csv"), row.names = FALSE)
write.csv(mk_range, file.path(results_dir, "sensitivity_loso_mk_range.csv"), row.names = FALSE)
write.csv(mk,       file.path(results_dir, "sensitivity_loso_mk_matrix.csv"))

print(mk_range)
message(sprintf("LOSO Pearson r: %.4f - %.4f", min(per_fold$Pearson_r), max(per_fold$Pearson_r)))
