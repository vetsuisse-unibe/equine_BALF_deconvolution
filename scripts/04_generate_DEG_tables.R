## Generate filtered DEG_tables from already-computed DEG_analysis results
## This is the standalone version of section 7b from 02_DEG_comparison.R

library(tidyverse)

# Working directory should contain DEG_analysis/ output from 02_DEG_comparison.R

deg_analysis_dir <- "DEG_analysis"
deg_tables_dir <- "DEG_tables"
dir.create(deg_tables_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: filter significant DEGs from full DESeq2 output
filter_sig <- function(df, lfc_col = "log2FoldChange", padj_col = "padj") {
  df[!is.na(df[[padj_col]]) & df[[padj_col]] < 0.05 & abs(df[[lfc_col]]) > 1, ]
}

# Cell types with both pseudobulk and deconvolution results
cell_types <- c("Neutrophils", "T_cells", "Mast_cells", "Mo_Ma", "Dendritic_cells")

for (ct in cell_types) {
  ct_clean <- gsub("_", "", ct)  # Match naming: T_cells -> Tcells

  # --- Pseudobulk: filter from FULL pseudobulk DESeq2 results ---
  pb_full <- read.csv(file.path(deg_analysis_dir, paste0("DEGs_pseudobulk_", ct, ".csv")))
  pb_sig <- filter_sig(pb_full)
  pb_sig <- pb_sig[order(pb_sig$padj), ]

  pb_out <- data.frame(gene = pb_sig$gene,
                       log2FC_pseudobulk = pb_sig$log2FoldChange,
                       padj_pseudobulk = pb_sig$padj)
  write.csv(pb_out, file.path(deg_tables_dir, paste0("DEGs_pseudobulk_", ct_clean, ".csv")), row.names = FALSE)

  # --- Deconvolution: filter from FULL deconvoluted DESeq2 results ---
  dc_full <- read.csv(file.path(deg_analysis_dir, paste0("DEGs_deconvoluted_", ct, ".csv")))
  dc_sig <- filter_sig(dc_full)
  dc_sig <- dc_sig[order(dc_sig$padj), ]

  dc_out <- data.frame(gene = dc_sig$gene,
                       log2FC_deconvoluted = dc_sig$log2FoldChange,
                       padj_deconvoluted = dc_sig$padj)
  write.csv(dc_out, file.path(deg_tables_dir, paste0("DEGs_deconvolution_", ct_clean, ".csv")), row.names = FALSE)

  # --- Overlap / *_only: use comparison table (common genes between methods) ---
  comp <- read.csv(file.path(deg_analysis_dir, paste0("DEGs_comparison_", ct, ".csv")))

  overlap_genes <- comp[comp$sig_pseudobulk == TRUE & comp$sig_deconvoluted == TRUE, ]
  pb_only_genes <- comp[comp$sig_pseudobulk == TRUE & comp$sig_deconvoluted == FALSE, ]
  dc_only_genes <- comp[comp$sig_deconvoluted == TRUE & comp$sig_pseudobulk == FALSE, ]

  # Overlap (significant in both)
  if (nrow(overlap_genes) > 0) {
    overlap_genes <- overlap_genes[order(overlap_genes$padj_pseudobulk), ]
    write.csv(
      data.frame(gene = overlap_genes$gene,
                 log2FC_pseudobulk = overlap_genes$log2FC_pseudobulk,
                 padj_pseudobulk = overlap_genes$padj_pseudobulk,
                 log2FC_deconvoluted = overlap_genes$log2FC_deconvoluted,
                 padj_deconvoluted = overlap_genes$padj_deconvoluted),
      file.path(deg_tables_dir, paste0("DEGs_overlap_", ct_clean, ".csv")),
      row.names = FALSE)
  } else {
    write.csv(
      data.frame(gene = character(0), log2FC_pseudobulk = numeric(0), padj_pseudobulk = numeric(0),
                 log2FC_deconvoluted = numeric(0), padj_deconvoluted = numeric(0)),
      file.path(deg_tables_dir, paste0("DEGs_overlap_", ct_clean, ".csv")),
      row.names = FALSE)
  }

  # Pseudobulk only (significant in pseudobulk, not deconvolution — within common genes)
  if (nrow(pb_only_genes) > 0) {
    pb_only_genes <- pb_only_genes[order(pb_only_genes$padj_pseudobulk), ]
    write.csv(
      data.frame(gene = pb_only_genes$gene,
                 log2FC_pseudobulk = pb_only_genes$log2FC_pseudobulk,
                 padj_pseudobulk = pb_only_genes$padj_pseudobulk),
      file.path(deg_tables_dir, paste0("DEGs_pseudobulk_only_", ct_clean, ".csv")),
      row.names = FALSE)
  } else {
    write.csv(
      data.frame(gene = character(0), log2FC_pseudobulk = numeric(0), padj_pseudobulk = numeric(0)),
      file.path(deg_tables_dir, paste0("DEGs_pseudobulk_only_", ct_clean, ".csv")),
      row.names = FALSE)
  }

  # Deconvolution only (significant in deconvolution, not pseudobulk — within common genes)
  if (nrow(dc_only_genes) > 0) {
    dc_only_genes <- dc_only_genes[order(dc_only_genes$padj_deconvoluted), ]
    write.csv(
      data.frame(gene = dc_only_genes$gene,
                 log2FC_deconvoluted = dc_only_genes$log2FC_deconvoluted,
                 padj_deconvoluted = dc_only_genes$padj_deconvoluted),
      file.path(deg_tables_dir, paste0("DEGs_deconvolution_only_", ct_clean, ".csv")),
      row.names = FALSE)
  } else {
    write.csv(
      data.frame(gene = character(0), log2FC_deconvoluted = numeric(0), padj_deconvoluted = numeric(0)),
      file.path(deg_tables_dir, paste0("DEGs_deconvolution_only_", ct_clean, ".csv")),
      row.names = FALSE)
  }

  cat(ct, ": PB=", nrow(pb_sig), " (was ", nrow(comp[comp$sig_pseudobulk == TRUE, ]), " from comparison)",
          " DC=", nrow(dc_sig), " (was ", nrow(comp[comp$sig_deconvoluted == TRUE, ]), " from comparison)",
          " Overlap=", nrow(overlap_genes),
          " PB-only=", nrow(pb_only_genes),
          " DC-only=", nrow(dc_only_genes))
}

# B cells: only pseudobulk (no deconvolution comparison available)
pb_full <- read.csv(file.path(deg_analysis_dir, "DEGs_pseudobulk_B_cells.csv"))
pb_sig <- filter_sig(pb_full)
pb_sig <- pb_sig[order(pb_sig$padj), ]
write.csv(
  data.frame(gene = pb_sig$gene,
             log2FC_pseudobulk = pb_sig$log2FoldChange,
             padj_pseudobulk = pb_sig$padj),
  file.path(deg_tables_dir, "DEGs_pseudobulk_Bcells.csv"),
  row.names = FALSE
)
cat("B cells: ", nrow(pb_sig), " pseudobulk DEGs (no deconvolution)")

## ===== Also regenerate DEG_comparison_summary.csv with corrected counts =====
message("\nRegenerating DEG_comparison_summary.csv with full DEG counts")

cell_types_summary <- c("Neutrophils", "T_cells", "Mast_cells", "Mo_Ma", "Dendritic_cells")
comparison_stats <- data.frame()

for (ct in cell_types_summary) {
  pb_full <- read.csv(file.path(deg_analysis_dir, paste0("DEGs_pseudobulk_", ct, ".csv")))
  dc_full <- read.csv(file.path(deg_analysis_dir, paste0("DEGs_deconvoluted_", ct, ".csv")))
  comp <- read.csv(file.path(deg_analysis_dir, paste0("DEGs_comparison_", ct, ".csv")))

  # Full counts
  pb_sig_full <- filter_sig(pb_full)
  dc_sig_full <- filter_sig(dc_full)

  # Common-gene counts
  pb_sig_common <- comp$gene[comp$sig_pseudobulk == TRUE]
  dc_sig_common <- comp$gene[comp$sig_deconvoluted == TRUE]

  # Overlap
  overlap <- length(intersect(pb_sig_common, dc_sig_common))
  union_size <- length(union(pb_sig_common, dc_sig_common))
  jaccard <- ifelse(union_size > 0, overlap / union_size, NA)

  # Log2FC correlation (all common genes)
  valid <- !is.na(comp$log2FC_pseudobulk) & !is.na(comp$log2FC_deconvoluted)
  lfc_pearson <- cor(comp$log2FC_pseudobulk[valid], comp$log2FC_deconvoluted[valid], method = "pearson")
  lfc_spearman <- cor(comp$log2FC_pseudobulk[valid], comp$log2FC_deconvoluted[valid], method = "spearman")

  # Direction concordance
  sig_either <- union(pb_sig_common, dc_sig_common)
  if (length(sig_either) > 0) {
    comp_sig <- comp[comp$gene %in% sig_either, ]
    pb_dir <- sign(comp_sig$log2FC_pseudobulk)
    dc_dir <- sign(comp_sig$log2FC_deconvoluted)
    valid_dir <- !is.na(pb_dir) & !is.na(dc_dir)
    direction_concordance <- mean(pb_dir[valid_dir] == dc_dir[valid_dir])
  } else {
    direction_concordance <- NA
  }

  # Cell type display name
  ct_display <- gsub("_", " ", ct)

  comparison_stats <- rbind(comparison_stats, data.frame(
    CellType = ct_display,
    N_Genes_Common = nrow(comp),
    N_DEG_Pseudobulk_Total = nrow(pb_sig_full),
    N_DEG_Deconvoluted_Total = nrow(dc_sig_full),
    N_DEG_Pseudobulk_Common = length(pb_sig_common),
    N_DEG_Deconvoluted_Common = length(dc_sig_common),
    N_DEG_Overlap = overlap,
    Jaccard = round(jaccard, 3),
    Log2FC_Pearson = round(lfc_pearson, 3),
    Log2FC_Spearman = round(lfc_spearman, 3),
    Direction_Concordance = round(direction_concordance, 3)
  ))
}

print(comparison_stats)
write.csv(comparison_stats, file.path(deg_analysis_dir, "DEG_comparison_summary.csv"), row.names = FALSE)
message("  Saved: DEG_analysis/DEG_comparison_summary.csv")

