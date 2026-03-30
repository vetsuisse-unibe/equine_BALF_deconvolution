############################################################
# DEG Comparison: scRNA Ground Truth vs BayesPrism Deconvolution
# Compare DEGs of SEA vs CTL between pseudo-bulk and deconvoluted expression
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(gridExtra)
  library(cowplot)
  library(ggbreak)
  library(BayesPrism)
})

## ===== Setup =====
data_dir <- "data"  # path to input files
set.seed(12345)

out_dir <- "DEG_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


theme_publication <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

## ===== Load Data =====
message("Loading data")

# Load scRNA reference
scRNA_ref <- readRDS(file.path(data_dir, "sc22_all_seed.rds"))
DefaultAssay(scRNA_ref) <- "RNA"

# Load BayesPrism results
load("bp.res.real_bulk_matched.RData")

# Create sample metadata from scRNA (keep original CTL/SEA labels)
sample_meta <- unique(scRNA_ref@meta.data[, c("sample_id", "disease_state")])
colnames(sample_meta)[2] <- "condition"
sample_meta$sample_id <- as.character(sample_meta$sample_id)
rownames(sample_meta) <- sample_meta$sample_id

print(table(sample_meta$condition))

# Get cell types
cell_types <- unique(scRNA_ref@meta.data$MajCellType)
print(paste("  Cell types: ", paste(cell_types, collapse = ", ")))

## ===== Create Pseudo-bulk from scRNA-seq (Ground Truth) =====
message("Creating pseudo-bulk from scRNA-seq")

pseudobulk_list <- list()

for (ct in cell_types) {
  message("  Processing: ", ct)
  
  cells_ct <- WhichCells(scRNA_ref, expression = MajCellType == ct)
  
  if (length(cells_ct) < 10) {
    next
  }
  
  seurat_ct <- subset(scRNA_ref, cells = cells_ct)
  
  agg_expr <- AggregateExpression(seurat_ct, 
                                  group.by = "sample_id",
                                  assays = "RNA",
                                  slot = "counts",
                                  return.seurat = FALSE)$RNA
  
  # Fix: Remove 'g' prefix added by AggregateExpression
  colnames(agg_expr) <- gsub("^g", "", colnames(agg_expr))
  
  samples_in_agg <- colnames(agg_expr)
  conditions_in_agg <- sample_meta$condition[match(samples_in_agg, sample_meta$sample_id)]
  
  if (length(unique(conditions_in_agg)) < 2) {
    next
  }
  
  pseudobulk_list[[ct]] <- list(
    counts = as.matrix(agg_expr),
    samples = samples_in_agg,
    conditions = conditions_in_agg
  )
}


## ===== Extract Cell-Type-Specific Expression from BayesPrism =====
message("Extracting deconvoluted expression from BayesPrism")

# Get Z matrix: samples x genes x cell_types
Z_matrix <- get.exp(bp = bp.res, state.or.type = "type")


## ===== B Cell Detection Comparison Table =====
message("Creating B cell detection comparison table")

# scRNA-seq B cell counts per sample
scrna_bcells <- as.data.frame(table(
  scRNA_ref$sample_id[scRNA_ref$MajCellType == "B cells"],
  scRNA_ref$disease_state[scRNA_ref$MajCellType == "B cells"]
))
colnames(scrna_bcells) <- c("sample_id", "condition", "scRNA_Bcell_count")
scrna_bcells <- scrna_bcells[scrna_bcells$scRNA_Bcell_count > 0, ]
scrna_bcells$sample_id <- as.character(scrna_bcells$sample_id)

# Deconvoluted B cell counts per sample
z_bcells <- t(Z_matrix[, , "B cells"])
b_cell_libsize <- colSums(round(z_bcells))

deconv_bcells <- data.frame(
  sample_id = names(b_cell_libsize),
  Deconv_Bcell_counts = as.numeric(b_cell_libsize)
)

# Merge tables
bcell_comparison <- merge(
  data.frame(
    sample_id = sample_meta$sample_id,
    condition = sample_meta$condition
  ),
  scrna_bcells[, c("sample_id", "scRNA_Bcell_count")],
  by = "sample_id",
  all.x = TRUE
)
bcell_comparison <- merge(bcell_comparison, deconv_bcells, by = "sample_id", all.x = TRUE)
bcell_comparison$scRNA_Bcell_count[is.na(bcell_comparison$scRNA_Bcell_count)] <- 0

# Order by condition and sample
bcell_comparison <- bcell_comparison[order(bcell_comparison$condition, bcell_comparison$sample_id), ]

# Add detection status
bcell_comparison$Deconv_detected <- ifelse(bcell_comparison$Deconv_Bcell_counts > 0, "Yes", "No")

message("\n  B Cell Detection Comparison:")
print(bcell_comparison)

# Save table
write.csv(bcell_comparison, file.path(out_dir, "Table_Bcell_detection_comparison.csv"), row.names = FALSE)

# Summary statistics
cat("\n  Summary:")
cat("    scRNA-seq B cells - CTL samples: ", 
        sum(bcell_comparison$scRNA_Bcell_count[bcell_comparison$condition == "CTL"]), " cells (",
        paste(bcell_comparison$scRNA_Bcell_count[bcell_comparison$condition == "CTL"], collapse = ", "), ")")
cat("    scRNA-seq B cells - SEA samples: ", 
        sum(bcell_comparison$scRNA_Bcell_count[bcell_comparison$condition == "SEA"]), " cells (",
        paste(bcell_comparison$scRNA_Bcell_count[bcell_comparison$condition == "SEA"], collapse = ", "), ")")
cat("    Deconvolution detected B cells in CTL: ", 
        sum(bcell_comparison$Deconv_Bcell_counts[bcell_comparison$condition == "CTL"] > 0), "/",
        sum(bcell_comparison$condition == "CTL"), " samples\n")
cat("    Deconvolution detected B cells in SEA: ", 
        sum(bcell_comparison$Deconv_Bcell_counts[bcell_comparison$condition == "SEA"] > 0), "/",
        sum(bcell_comparison$condition == "SEA"), " samples\n")

## ===== Run DESeq2 for Each Cell Type =====
message("Running differential expression analysis")

run_deseq2 <- function(count_matrix, col_data, design_formula = ~ condition,
                       contrast = c("condition", "SEA", "CTL")) {
  
  count_matrix <- count_matrix[rowSums(count_matrix) >= 10, ]
  count_matrix <- round(count_matrix)
  count_matrix[count_matrix < 0] <- 0
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = design_formula
  )
  
  # poscounts works for all data (handles zeros)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, quiet = TRUE)
  
  res <- results(dds, contrast = contrast)
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  
  return(list(dds = dds, results = res))
}

deg_results <- list()

for (ct in names(pseudobulk_list)) {
  message("  ", ct)
  
  # Find matching cell type in Z matrix
  ct_clean <- ct
  if (!ct_clean %in% dimnames(Z_matrix)[[3]]) {
    ct_matches <- grep(ct, dimnames(Z_matrix)[[3]], value = TRUE, ignore.case = TRUE)
    if (length(ct_matches) == 0) {
      message("    Cell type not found in Z matrix, skipping")
      next
    }
    ct_clean <- ct_matches[1]
  }
  
  # --- Ground Truth DEGs (Pseudo-bulk) ---
  pb_data <- pseudobulk_list[[ct]]
  pb_coldata <- data.frame(
    row.names = pb_data$samples,
    condition = factor(pb_data$conditions, levels = c("CTL", "SEA"))
  )
  
  tryCatch({
    pb_res <- run_deseq2(pb_data$counts, pb_coldata)
    
    # --- Check if B cells - skip deconvolution comparison ---
    if (ct == "B cells") {
      message("    Skipping deconvolution DEG - CTL samples have zero B cell expression")
      
      deg_results[[ct]] <- list(
        pseudobulk = pb_res$results,
        deconvoluted = NULL,
        note = "Deconvolution comparison skipped: BayesPrism estimated zero B cell expression in all CTL samples"
      )
      
      message("    Pseudo-bulk DEGs: ", sum(pb_res$results$padj < 0.05, na.rm = TRUE), " (padj < 0.05)")
      next
    }
    
    # --- Deconvoluted DEGs (Z matrix) ---
    
    # Extract: Z is samples x genes x cell_types -> transpose to genes x samples
    z_expr <- t(Z_matrix[, , ct_clean])  # Now genes x samples
    
    # Check for zero library samples
    lib_sizes <- colSums(round(z_expr))
    ctl_has_data <- sum(lib_sizes[sample_meta$condition[match(names(lib_sizes), sample_meta$sample_id)] == "CTL"] > 0) > 0
    sea_has_data <- sum(lib_sizes[sample_meta$condition[match(names(lib_sizes), sample_meta$sample_id)] == "SEA"] > 0) > 0
    
    if (!ctl_has_data || !sea_has_data) {
      print("    Skipping deconvolution DEG - one condition has no expression data")
      deg_results[[ct]] <- list(
        pseudobulk = pb_res$results,
        deconvoluted = NULL,
        note = "Deconvolution comparison skipped: one condition has zero expression"
      )
      message("    Pseudo-bulk DEGs: ", sum(pb_res$results$padj < 0.05, na.rm = TRUE), " (padj < 0.05)")
      next
    }
    
    # Match samples
    common_samples <- intersect(colnames(z_expr), pb_data$samples)
    
    if (length(common_samples) < 4) {
      message("    Too few common samples, skipping")
      next
    }
    
    z_expr <- z_expr[, common_samples]
    z_coldata <- data.frame(
      row.names = common_samples,
      condition = factor(sample_meta$condition[match(common_samples, sample_meta$sample_id)],
                         levels = c("CTL", "SEA"))
    )
    
    z_res <- run_deseq2(z_expr, z_coldata)
    
    deg_results[[ct]] <- list(
      pseudobulk = pb_res$results,
      deconvoluted = z_res$results
    )
    
    
  }, error = function(e) {
  })
}

## ===== Compare DEGs =====
message("Comparing DEGs between methods")

comparison_stats <- data.frame()

for (ct in names(deg_results)) {
  message("  ", ct)
  
  pb_degs <- deg_results[[ct]]$pseudobulk
  dc_degs <- deg_results[[ct]]$deconvoluted
  
  # Skip if no deconvoluted results (e.g., B cells)
  if (is.null(dc_degs)) {
    next
  }
  
  # Common genes
  common_genes <- intersect(pb_degs$gene, dc_degs$gene)
  
  pb_common <- pb_degs[pb_degs$gene %in% common_genes, ]
  dc_common <- dc_degs[dc_degs$gene %in% common_genes, ]
  
  # Align by gene
  rownames(pb_common) <- pb_common$gene
  rownames(dc_common) <- dc_common$gene
  pb_common <- pb_common[common_genes, ]
  dc_common <- dc_common[common_genes, ]
  
  # Significant DEGs (padj < 0.05, |log2FC| > 1)
  # Full counts (all genes tested per method)
  pb_sig_full <- pb_degs$gene[!is.na(pb_degs$padj) & pb_degs$padj < 0.05 & abs(pb_degs$log2FoldChange) > 1]
  dc_sig_full <- dc_degs$gene[!is.na(dc_degs$padj) & dc_degs$padj < 0.05 & abs(dc_degs$log2FoldChange) > 1]
  # Counts within common genes (for overlap/Jaccard/concordance)
  pb_sig <- pb_common$gene[!is.na(pb_common$padj) & pb_common$padj < 0.05 & abs(pb_common$log2FoldChange) > 1]
  dc_sig <- dc_common$gene[!is.na(dc_common$padj) & dc_common$padj < 0.05 & abs(dc_common$log2FoldChange) > 1]
  
  # Overlap metrics
  intersection <- length(intersect(pb_sig, dc_sig))
  union_size <- length(union(pb_sig, dc_sig))
  jaccard <- ifelse(union_size > 0, intersection / union_size, NA)
  
  # Correlation of log2FC
  valid_idx <- !is.na(pb_common$log2FoldChange) & !is.na(dc_common$log2FoldChange)
  lfc_cor <- cor(pb_common$log2FoldChange[valid_idx], 
                 dc_common$log2FoldChange[valid_idx], 
                 method = "pearson")
  lfc_spearman <- cor(pb_common$log2FoldChange[valid_idx], 
                      dc_common$log2FoldChange[valid_idx], 
                      method = "spearman")
  
  # Direction concordance
  sig_either <- union(pb_sig, dc_sig)
  if (length(sig_either) > 0) {
    pb_direction <- sign(pb_common[sig_either, "log2FoldChange"])
    dc_direction <- sign(dc_common[sig_either, "log2FoldChange"])
    valid_dir <- !is.na(pb_direction) & !is.na(dc_direction)
    direction_concordance <- mean(pb_direction[valid_dir] == dc_direction[valid_dir])
  } else {
    direction_concordance <- NA
  }
  
  stats_row <- data.frame(
    CellType = ct,
    N_Genes_Common = length(common_genes),
    N_DEG_Pseudobulk_Total = length(pb_sig_full),
    N_DEG_Deconvoluted_Total = length(dc_sig_full),
    N_DEG_Pseudobulk_Common = length(pb_sig),
    N_DEG_Deconvoluted_Common = length(dc_sig),
    N_DEG_Overlap = intersection,
    Jaccard = round(jaccard, 3),
    Log2FC_Pearson = round(lfc_cor, 3),
    Log2FC_Spearman = round(lfc_spearman, 3),
    Direction_Concordance = round(direction_concordance, 3)
  )
  
  comparison_stats <- rbind(comparison_stats, stats_row)
  
  # Store for plotting
  deg_results[[ct]]$comparison <- data.frame(
    gene = common_genes,
    log2FC_pseudobulk = pb_common$log2FoldChange,
    log2FC_deconvoluted = dc_common$log2FoldChange,
    padj_pseudobulk = pb_common$padj,
    padj_deconvoluted = dc_common$padj,
    sig_pseudobulk = common_genes %in% pb_sig,
    sig_deconvoluted = common_genes %in% dc_sig
  )
}

print(comparison_stats)
write.csv(comparison_stats, file.path(out_dir, "DEG_comparison_summary.csv"), row.names = FALSE)

## ===== Generate Figures =====
message("Generating figures")

# --- Figure 5: Log2FC Correlation Scatter Plots (Single Shared Legend) ---
scatter_plots <- list()

for (ct in names(deg_results)) {
  if (is.null(deg_results[[ct]]$comparison)) next

  df <- deg_results[[ct]]$comparison
  df <- df[!is.na(df$log2FC_pseudobulk) & !is.na(df$log2FC_deconvoluted), ]

  df$Category <- case_when(
    df$sig_pseudobulk & df$sig_deconvoluted ~ "Both",
    df$sig_pseudobulk ~ "Experimental scRNA-seq only",
    df$sig_deconvoluted ~ "Deconvoluted only",
    TRUE ~ "Not significant"
  )
  df$Category <- factor(df$Category, levels = c("Both", "Experimental scRNA-seq only", "Deconvoluted only", "Not significant"))

  r <- cor(df$log2FC_pseudobulk, df$log2FC_deconvoluted, use = "complete.obs")

  p <- ggplot(df, aes(x = log2FC_pseudobulk, y = log2FC_deconvoluted, color = Category)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8,
                data = df, aes(x = log2FC_pseudobulk, y = log2FC_deconvoluted), inherit.aes = FALSE) +
    scale_color_manual(values = c("Both" = "#E64B35", "Experimental scRNA-seq only" = "#4DBBD5",
                                  "Deconvoluted only" = "#00A087", "Not significant" = "gray70")) +
    annotate("text", x = min(df$log2FC_pseudobulk, na.rm = TRUE),
             y = max(df$log2FC_deconvoluted, na.rm = TRUE),
             label = paste0("r = ", round(r, 3)), hjust = 0, vjust = 1, size = 4, fontface = "bold") +
    theme_publication +
    labs(title = gsub("Mo/Ma", "MoMa", ct), x = "log2FC (Experimental scRNA-seq)", y = "log2FC (Deconvoluted)", color = "DEG Status") +
    theme(legend.position = "none")  # Remove individual legends

  scatter_plots[[ct]] <- p
}

if (length(scatter_plots) > 0) {
  n_plots <- length(scatter_plots)
  ncol <- min(3, n_plots)
  nrow <- ceiling(n_plots / ncol)

  # Extract legend from one plot (temporarily add legend back)
  legend_plot <- scatter_plots[[1]] +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1, title = "DEG Status"))

  # Extract the legend using cowplot
  shared_legend <- get_legend(legend_plot)

  # Combine plots without legends
  plots_combined <- plot_grid(plotlist = scatter_plots, ncol = ncol, align = "hv")

  # Add shared legend at bottom
  fig5 <- plot_grid(plots_combined, shared_legend, ncol = 1, rel_heights = c(1, 0.08))

  ggsave(file.path(out_dir, "Figure5_DEG_log2FC_correlation.pdf"), fig5, width = 5 * ncol, height = 5 * nrow + 0.5, device = cairo_pdf)
  ggsave(file.path(out_dir, "Figure5_DEG_log2FC_correlation.png"), fig5, width = 5 * ncol, height = 5 * nrow + 0.5, dpi = 300)
  message("  Saved: Figure5_DEG_log2FC_correlation.pdf/png")
}

# --- Figure 2: DEG Comparison Bar Plot ---
bar_data <- data.frame()

for (ct in names(deg_results)) {
  df <- deg_results[[ct]]$comparison
  if (is.null(df)) next
  
  pb_sig <- df$gene[df$sig_pseudobulk]
  dc_sig <- df$gene[df$sig_deconvoluted]
  
  overlap <- length(intersect(pb_sig, dc_sig))
  pb_only <- length(setdiff(pb_sig, dc_sig))
  dc_only <- length(setdiff(dc_sig, pb_sig))
  
  bar_data <- rbind(bar_data, data.frame(
    CellType = ct,
    Category = c("Pseudo-bulk only", "Overlap", "Deconvoluted only"),
    Count = c(pb_only, overlap, dc_only)
  ))
}

if (nrow(bar_data) > 0) {
  # Set factor order
  bar_data$Category <- factor(bar_data$Category, 
                              levels = c("Pseudo-bulk only", "Overlap", "Deconvoluted only"))
  
  # Publication colors
  bar_colors <- c("Pseudo-bulk only" = "#005CAB", 
                  "Overlap" = "#F46A25", 
                  "Deconvoluted only" = "#3D3D3D")
  
  # Plot
  fig_bar <- ggplot(bar_data, aes(x = CellType, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = Count), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = bar_colors) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Comparison of DEGs: Pseudo-bulk vs Deconvoluted",
      subtitle = "Significant DEGs (padj < 0.05, |log2FC| > 1)",
      x = "Cell Type",
      y = "Number of DEGs"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  ggsave(file.path(out_dir, "Figure_DEG_comparison_barplot.pdf"), fig_bar, width = 10, height = 6, device = cairo_pdf)
  ggsave(file.path(out_dir, "Figure_DEG_comparison_barplot.png"), fig_bar, width = 10, height = 6, dpi = 300)
  message("  Saved: Figure_DEG_comparison_barplot.pdf/png")
}

# --- Figure 3: Correlation Metrics Heatmap ---
if (nrow(comparison_stats) > 0) {
  metrics_mat <- comparison_stats %>%
    select(CellType, Jaccard, Log2FC_Pearson, Log2FC_Spearman, Direction_Concordance) %>%
    column_to_rownames("CellType") %>%
    as.matrix()
  
  cairo_pdf(file.path(out_dir, "Figure_DEG_metrics_heatmap.pdf"), width = 8, height = 6)
  pheatmap(metrics_mat, 
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.2f",
           color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
           main = "DEG Comparison Metrics",
           fontsize = 12)
  dev.off()
  message("  Saved: Figure_DEG_metrics_heatmap.pdf")
}

# --- Figure S: B Cell Detection Comparison ---
bcell_plot_df <- bcell_comparison %>%
  pivot_longer(cols = c(scRNA_Bcell_count, Deconv_Bcell_counts),
               names_to = "Method", values_to = "Count") %>%
  mutate(Method = ifelse(Method == "scRNA_Bcell_count", "scRNA-seq", "Deconvoluted"))


fig_bcell <- ggplot(bcell_plot_df, aes(x = sample_id, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  facet_wrap(~condition, scales = "free_x") +
  scale_fill_manual(values = c("scRNA-seq" = "#4DBBD5", "Deconvoluted" = "#E64B35")) +
  scale_y_break(c(300, 5000), scales = 0.5) +  # Break between 300 and 5000
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "B Cell Detection: scRNA-seq vs Deconvolution",
       subtitle = "Deconvolution failed to detect B cells in CTL samples due to low abundance",
       x = "Sample ID", y = "B Cell Count / Expression", fill = "Method")


ggsave(file.path(out_dir, "Figure_Bcell_detection_comparison.pdf"), fig_bcell, width = 10, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_Bcell_detection_comparison.png"), fig_bcell, width = 10, height = 6, dpi = 300)

## ===== Export DEG Lists =====
message("Exporting DEG lists")

for (ct in names(deg_results)) {
  # Clean cell type name for filename
  ct_filename <- gsub("[/ ]", "_", ct)  # Replace both / and space with _
  
  # Pseudo-bulk DEGs (always available)
  pb_degs <- deg_results[[ct]]$pseudobulk
  pb_degs <- pb_degs[order(pb_degs$padj), ]
  write.csv(pb_degs, file.path(out_dir, paste0("DEGs_pseudobulk_", ct_filename, ".csv")), row.names = FALSE)
  
  # Deconvoluted DEGs (may be NULL for B cells)
  dc_degs <- deg_results[[ct]]$deconvoluted
  if (!is.null(dc_degs)) {
    dc_degs <- dc_degs[order(dc_degs$padj), ]
    write.csv(dc_degs, file.path(out_dir, paste0("DEGs_deconvoluted_", ct_filename, ".csv")), row.names = FALSE)
  }
  
  # Comparison table (may be NULL for B cells)
  if (!is.null(deg_results[[ct]]$comparison)) {
    comp <- deg_results[[ct]]$comparison
    comp <- comp[order(comp$padj_pseudobulk), ]
    write.csv(comp, file.path(out_dir, paste0("DEGs_comparison_", ct_filename, ".csv")), row.names = FALSE)
  }
}

## ===== Generate Filtered DEG Tables =====
## NOTE: Pseudobulk and deconvolution DEGs are filtered from their own full
##       DESeq2 results (not from the comparison table, which only contains
##       genes common to both methods). Overlap / *_only files use the
##       comparison table as intended.
message("Generating filtered DEG tables")

deg_tables_dir <- "DEG_tables"
dir.create(deg_tables_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: filter significant DEGs (padj < 0.05 and |log2FC| > 1)
filter_sig <- function(df) {
  df[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1, ]
}

for (ct in names(deg_results)) {
  ct_clean <- gsub("[/ ]", "", ct)

  pb_full <- deg_results[[ct]]$pseudobulk
  dc_full <- deg_results[[ct]]$deconvoluted
  comp    <- deg_results[[ct]]$comparison

  # --- Pseudobulk: filter from FULL pseudobulk DESeq2 results ---
  if (!is.null(pb_full)) {
    pb_sig <- filter_sig(pb_full)
    pb_sig <- pb_sig[order(pb_sig$padj), ]
    write.csv(
      data.frame(gene = pb_sig$gene,
                 log2FC_pseudobulk = pb_sig$log2FoldChange,
                 padj_pseudobulk = pb_sig$padj),
      file.path(deg_tables_dir, paste0("DEGs_pseudobulk_", ct_clean, ".csv")),
      row.names = FALSE
    )
  }

  # --- Deconvolution: filter from FULL deconvoluted DESeq2 results ---
  if (!is.null(dc_full)) {
    dc_sig <- filter_sig(dc_full)
    dc_sig <- dc_sig[order(dc_sig$padj), ]
    write.csv(
      data.frame(gene = dc_sig$gene,
                 log2FC_deconvoluted = dc_sig$log2FoldChange,
                 padj_deconvoluted = dc_sig$padj),
      file.path(deg_tables_dir, paste0("DEGs_deconvolution_", ct_clean, ".csv")),
      row.names = FALSE
    )
  }

  # --- Overlap / *_only: use comparison table (common genes between methods) ---
  if (!is.null(comp)) {
    overlap_genes <- comp[comp$sig_pseudobulk == TRUE & comp$sig_deconvoluted == TRUE, ]
    pb_only_genes <- comp[comp$sig_pseudobulk == TRUE & comp$sig_deconvoluted == FALSE, ]
    dc_only_genes <- comp[comp$sig_deconvoluted == TRUE & comp$sig_pseudobulk == FALSE, ]

    # Overlap
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

    # Pseudobulk only (within common genes)
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

    # Deconvolution only (within common genes)
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

    n_pb <- if (!is.null(pb_full)) nrow(filter_sig(pb_full)) else 0
    n_dc <- if (!is.null(dc_full)) nrow(filter_sig(dc_full)) else 0
    message("  ", ct, ": PB=", n_pb, " DC=", n_dc,
            " Overlap=", nrow(overlap_genes), " PB-only=", nrow(pb_only_genes),
            " DC-only=", nrow(dc_only_genes))
  } else {
    n_pb <- if (!is.null(pb_full)) nrow(filter_sig(pb_full)) else 0
    message("  ", ct, ": PB=", n_pb, " (no deconvolution)")
  }
}

message("  Filtered DEG tables saved to: ", deg_tables_dir)

## ===== Summary =====
print("\n=== Summary ===")
print("Figures generated:")
print("\nTables generated:")

print("\n=== Key Metrics ===")
print(comparison_stats)

message("\n=== B Cell Note ===")
message("B cells were excluded from DEG comparison between pseudo-bulk and deconvoluted data.")
message("Reason: BayesPrism estimated zero B cell expression in all CTL samples (n=5),")
message("        despite scRNA-seq detecting 12-27 B cells per CTL sample.")
message("This reflects a sensitivity limitation of deconvolution for low-abundance cell types.")
message("See Table_Bcell_detection_comparison.csv for details.")

## ===== Bulk RNA-seq DEG Analysis (CTL vs SEA) =====
message("Running DESeq2 on bulk RNA-seq (CTL vs SEA)")

# Load bulk RNA-seq data
bulk_counts <- read.table(file.path(data_dir, "bulkRNAseq/counts_extended/merged_count_matrix_renamed.txt"),
                          header = TRUE, row.names = 1, check.names = FALSE)

# Remove spike-ins and metadata rows
spike_rows <- grep("^(ERCC|SIRV|__|Assigned|Unassigned|No_feature|Ambiguous)", rownames(bulk_counts))
if (length(spike_rows) > 0) bulk_counts <- bulk_counts[-spike_rows, ]

# Filter low-count genes
bulk_counts <- bulk_counts[rowSums(bulk_counts) >= 10, ]

# Harmonize gene IDs and sample IDs
rownames(bulk_counts) <- sub("^gene-", "", rownames(bulk_counts))
colnames(bulk_counts) <- sub("^A", "", colnames(bulk_counts))


# Match samples with metadata
bulk_samples <- colnames(bulk_counts)
matched_bulk_samples <- intersect(bulk_samples, sample_meta$sample_id)

if (length(matched_bulk_samples) < 4) {
  stop("Too few matched samples between bulk RNA-seq and metadata")
}

bulk_counts_matched <- bulk_counts[, matched_bulk_samples]
bulk_coldata <- data.frame(
  row.names = matched_bulk_samples,
  condition = factor(sample_meta$condition[match(matched_bulk_samples, sample_meta$sample_id)],
                     levels = c("CTL", "SEA"))
)

print(table(bulk_coldata$condition))

# Run DESeq2 on bulk RNA-seq
message("  Running DESeq2")
bulk_dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(bulk_counts_matched),
  colData = bulk_coldata,
  design = ~ condition
)

bulk_dds <- estimateSizeFactors(bulk_dds, type = "poscounts")
bulk_dds <- DESeq(bulk_dds, quiet = TRUE)

bulk_res <- results(bulk_dds, contrast = c("condition", "SEA", "CTL"))
bulk_res_df <- as.data.frame(bulk_res)
bulk_res_df$gene <- rownames(bulk_res_df)

# Significant bulk DEGs
bulk_sig <- bulk_res_df$gene[!is.na(bulk_res_df$padj) &
                              bulk_res_df$padj < 0.05 &
                              abs(bulk_res_df$log2FoldChange) > 1]


# Save bulk DEGs
write.csv(bulk_res_df[order(bulk_res_df$padj), ], file.path(out_dir, "DEGs_bulk_RNAseq.csv"), row.names = FALSE)

## ===== Three-Way DEG Overlap (Bulk vs Deconvolution vs scRNA-seq) =====
message("Computing three-way DEG overlap")

# Combine all scRNA-seq pseudo-bulk DEGs across cell types
all_pseudobulk_degs <- c()
for (ct in names(deg_results)) {
  pb_degs <- deg_results[[ct]]$pseudobulk
  pb_sig <- pb_degs$gene[!is.na(pb_degs$padj) &
                          pb_degs$padj < 0.05 &
                          abs(pb_degs$log2FoldChange) > 1]
  all_pseudobulk_degs <- union(all_pseudobulk_degs, pb_sig)
}

# Combine all deconvolution DEGs across cell types
all_deconv_degs <- c()
for (ct in names(deg_results)) {
  dc_degs <- deg_results[[ct]]$deconvoluted
  if (!is.null(dc_degs)) {
    dc_sig <- dc_degs$gene[!is.na(dc_degs$padj) &
                            dc_degs$padj < 0.05 &
                            abs(dc_degs$log2FoldChange) > 1]
    all_deconv_degs <- union(all_deconv_degs, dc_sig)
  }
}

# Bulk DEGs (already computed)

# Compute overlaps
overlap_bulk_scrna <- intersect(bulk_sig, all_pseudobulk_degs)
overlap_bulk_deconv <- intersect(bulk_sig, all_deconv_degs)
overlap_scrna_deconv <- intersect(all_pseudobulk_degs, all_deconv_degs)
overlap_all_three <- Reduce(intersect, list(bulk_sig, all_pseudobulk_degs, all_deconv_degs))

print("\n  === Three-Way Overlap Summary ===")
print("  Bulk ∩ scRNA-seq: ", length(overlap_bulk_scrna))
print("  Bulk ∩ Deconvolution: ", length(overlap_bulk_deconv))
print("  scRNA-seq ∩ Deconvolution: ", length(overlap_scrna_deconv))
print("  All three methods: ", length(overlap_all_three))

# Unique to each method
only_bulk <- setdiff(setdiff(bulk_sig, all_pseudobulk_degs), all_deconv_degs)
only_scrna <- setdiff(setdiff(all_pseudobulk_degs, bulk_sig), all_deconv_degs)
only_deconv <- setdiff(setdiff(all_deconv_degs, bulk_sig), all_pseudobulk_degs)

print("\n  Unique to each method:")
print("  Bulk only: ", length(only_bulk))
print("  scRNA-seq only: ", length(only_scrna))
print("  Deconvolution only: ", length(only_deconv))

# Create summary table
threeway_summary <- data.frame(
  Category = c("Bulk RNA-seq", "scRNA-seq (pseudo-bulk)", "Deconvolution",
               "Bulk ∩ scRNA-seq", "Bulk ∩ Deconvolution", "scRNA-seq ∩ Deconvolution",
               "All three methods",
               "Bulk only", "scRNA-seq only", "Deconvolution only"),
  Count = c(length(bulk_sig), length(all_pseudobulk_degs), length(all_deconv_degs),
            length(overlap_bulk_scrna), length(overlap_bulk_deconv), length(overlap_scrna_deconv),
            length(overlap_all_three),
            length(only_bulk), length(only_scrna), length(only_deconv))
)

print(threeway_summary)
write.csv(threeway_summary, file.path(out_dir, "DEG_threeway_overlap_summary.csv"), row.names = FALSE)

# Save gene lists
write.csv(data.frame(gene = overlap_all_three), file.path(out_dir, "DEG_overlap_all_three_methods.csv"), row.names = FALSE)
write.csv(data.frame(gene = overlap_bulk_scrna), file.path(out_dir, "DEG_overlap_bulk_scrna.csv"), row.names = FALSE)
write.csv(data.frame(gene = overlap_bulk_deconv), file.path(out_dir, "DEG_overlap_bulk_deconv.csv"), row.names = FALSE)

## ===== Bulk vs Cell-Type-Specific DEG Comparison =====
message("Comparing bulk DEGs to each cell type")

bulk_vs_celltype <- data.frame()

for (ct in names(deg_results)) {
  # Pseudo-bulk DEGs for this cell type
  pb_degs <- deg_results[[ct]]$pseudobulk
  pb_sig <- pb_degs$gene[!is.na(pb_degs$padj) &
                          pb_degs$padj < 0.05 &
                          abs(pb_degs$log2FoldChange) > 1]

  # Deconvolution DEGs for this cell type
  dc_degs <- deg_results[[ct]]$deconvoluted
  if (!is.null(dc_degs)) {
    dc_sig <- dc_degs$gene[!is.na(dc_degs$padj) &
                            dc_degs$padj < 0.05 &
                            abs(dc_degs$log2FoldChange) > 1]
  } else {
    dc_sig <- character(0)
  }

  # Overlaps with bulk
  overlap_bulk_pb <- length(intersect(bulk_sig, pb_sig))
  overlap_bulk_dc <- length(intersect(bulk_sig, dc_sig))

  # Jaccard indices
  jaccard_bulk_pb <- ifelse(length(union(bulk_sig, pb_sig)) > 0,
                            overlap_bulk_pb / length(union(bulk_sig, pb_sig)), NA)
  jaccard_bulk_dc <- ifelse(length(union(bulk_sig, dc_sig)) > 0,
                            overlap_bulk_dc / length(union(bulk_sig, dc_sig)), NA)

  bulk_vs_celltype <- rbind(bulk_vs_celltype, data.frame(
    CellType = ct,
    N_scRNA_DEGs = length(pb_sig),
    N_Deconv_DEGs = length(dc_sig),
    Overlap_Bulk_scRNA = overlap_bulk_pb,
    Overlap_Bulk_Deconv = overlap_bulk_dc,
    Jaccard_Bulk_scRNA = round(jaccard_bulk_pb, 3),
    Jaccard_Bulk_Deconv = round(jaccard_bulk_dc, 3)
  ))
}

# Add bulk total row
bulk_vs_celltype <- rbind(bulk_vs_celltype, data.frame(
  CellType = "BULK TOTAL",
  N_scRNA_DEGs = NA,
  N_Deconv_DEGs = NA,
  Overlap_Bulk_scRNA = NA,
  Overlap_Bulk_Deconv = NA,
  Jaccard_Bulk_scRNA = NA,
  Jaccard_Bulk_Deconv = NA
))
bulk_vs_celltype$N_Bulk_DEGs <- length(bulk_sig)

print(bulk_vs_celltype)
write.csv(bulk_vs_celltype, file.path(out_dir, "DEG_bulk_vs_celltype_comparison.csv"), row.names = FALSE)

## ===== Updated Summary =====
print("\n=== Updated Summary ===")
print(paste0("\nOutput directory: ", out_dir))
print("\nFiles generated:")
print("  Filtered DEG tables (DEG_tables/):")

print("\n=== Three-Way DEG Overlap ===")
print(paste0("Bulk RNA-seq DEGs: ", length(bulk_sig)))
print(paste0("scRNA-seq DEGs (all cell types combined): ", length(all_pseudobulk_degs)))
print(paste0("Deconvolution DEGs (all cell types combined): ", length(all_deconv_degs)))
print(paste0("Overlap (all three): ", length(overlap_all_three)))

if (length(overlap_all_three) > 0) {
  print("\nGenes found in all three methods:")
  print(head(overlap_all_three, 20))
  if (length(overlap_all_three) > 20) {
    print(paste0("  ... and ", length(overlap_all_three) - 20, " more"))
  }
}
