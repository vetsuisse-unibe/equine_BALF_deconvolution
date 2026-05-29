############################################################
# mRNA-content correction: sensitivity to disease state
# Recompute m_k separately for CTL and SEA, re-apply by condition,
# and compare with the pooled baseline (Supplementary Table S6)
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

# val is long-format: Sample, CellType, BayesPrism (theta), Truth (cell-count prop)
theta <- val %>% select(Sample, CellType, BayesPrism) %>%
  pivot_wider(names_from = CellType, values_from = BayesPrism) %>%
  column_to_rownames("Sample") %>% as.matrix()
truth <- val %>% select(Sample, CellType, Truth) %>%
  pivot_wider(names_from = CellType, values_from = Truth) %>%
  column_to_rownames("Sample") %>% as.matrix()

umi$Sample <- sub("^A", "", as.character(umi$SampleID))  # "A30" -> "30"
samples   <- rownames(theta)
celltypes <- colnames(theta)

# divide theta by m_k and renormalise to sum 1
correct <- function(theta, m) {
  out <- sweep(theta, 2, m[colnames(theta)], "/")
  out / rowSums(out)
}
score <- function(est) {
  obs <- truth[rownames(est), colnames(est)]
  data.frame(Pearson_r    = cor(c(est), c(obs)),
             Spearman_rho = cor(c(est), c(obs), method = "spearman"),
             MAE          = mean(abs(est - obs)))
}

## ===== Pooled vs condition-specific m_k =====
m_pooled <- tapply(umi$UMI, umi$CellType, mean)
m_ctl    <- tapply(umi$UMI[umi$Condition == "CTL"], umi$CellType[umi$Condition == "CTL"], mean)
m_sea    <- tapply(umi$UMI[umi$Condition == "SEA"], umi$CellType[umi$Condition == "SEA"], mean)

cond <- umi %>% distinct(Sample, Condition) %>% filter(Sample %in% samples)
theta_cond <- theta
for (s in samples) {
  m <- if (cond$Condition[cond$Sample == s] == "CTL") m_ctl else m_sea
  v <- theta[s, ] / m[celltypes]
  theta_cond[s, ] <- v / sum(v)
}

## ===== Outputs =====
mk_tbl <- data.frame(CellType = names(m_ctl),
                     m_CTL = as.numeric(m_ctl),
                     m_SEA = as.numeric(m_sea[names(m_ctl)]))
mk_tbl$SEA_over_CTL <- round(mk_tbl$m_SEA / mk_tbl$m_CTL, 2)

overall <- rbind(cbind(Variant = "pooled (baseline)",  score(correct(theta, m_pooled))),
                 cbind(Variant = "condition_specific", score(theta_cond)))

write.csv(mk_tbl,  file.path(results_dir, "mRNA_scaling_condition_specific.csv"), row.names = FALSE)
write.csv(overall, file.path(results_dir, "sensitivity_condition_overall.csv"),   row.names = FALSE)

print(mk_tbl)
print(overall)
