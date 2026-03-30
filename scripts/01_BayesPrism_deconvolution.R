############################################################
# BayesPrism Deconvolution of Real Bulk RNA-seq
# Compare BayesPrism estimates to experimental scRNA-seq  (matched samples)
############################################################

suppressPackageStartupMessages({
  library(BayesPrism)
  library(Seurat)
  library(tidyverse)
})

## ===== Setup =====
data_dir <- "data"  # path to input files
set.seed(12345)

## ===== Load scRNA reference =====
scRNA_ref <- readRDS(file.path(data_dir, "sc22_all_seed.rds"))
DefaultAssay(scRNA_ref) <- "RNA"

sc.dat <- t(as.matrix(GetAssayData(scRNA_ref, assay = "RNA", layer = "counts")))

cell.type <- as.character(scRNA_ref@meta.data$MajCellType)
names(cell.type) <- rownames(scRNA_ref@meta.data)

sample_id <- as.character(scRNA_ref@meta.data$sample_id)
names(sample_id) <- rownames(scRNA_ref@meta.data)

merged_meta <- read.csv(file.path(data_dir, "merged_metadata.csv"))
cell.state <- paste(merged_meta$MajorCellType, merged_meta$cellSubtype)
names(cell.state) <- merged_meta$CellID
cell.state <- as.character(cell.state[rownames(sc.dat)])


## ===== Load real bulk RNA-seq =====
message("Loading real bulk RNA-seq")
bulk_counts <- read.table(file.path(data_dir, "bulkRNAseq/counts_extended/merged_count_matrix_renamed.txt"),
                          header = TRUE, row.names = 1, check.names = FALSE)

# Remove spike-ins and metadata rows
spike_rows <- grep("^(ERCC|SIRV|__|Assigned|Unassigned|No_feature|Ambiguous)", rownames(bulk_counts))
if (length(spike_rows) > 0) bulk_counts <- bulk_counts[-spike_rows, ]

# Filter low-count genes
bulk_counts <- bulk_counts[rowSums(bulk_counts) >= 10, ]

# Harmonize gene IDs: strip "gene-" prefix from bulk row names
rownames(bulk_counts) <- sub("^gene-", "", rownames(bulk_counts))

# Harmonize sample IDs: strip "A" prefix from bulk column names
colnames(bulk_counts) <- sub("^A", "", colnames(bulk_counts))


## ===== Bulk RNA-seq Quality Control =====
message("Bulk RNA-seq quality control")

# Calculate QC metrics
bulk_qc <- data.frame(
  Sample = colnames(bulk_counts),
  LibrarySize = colSums(bulk_counts),
  GenesDetected = colSums(bulk_counts > 0),
  ZeroInflation = colSums(bulk_counts == 0) / nrow(bulk_counts) * 100
)

# Proportion of reads in top 10 genes per sample
top10_prop <- sapply(colnames(bulk_counts), function(s) {
  counts <- bulk_counts[, s]
  sum(sort(counts, decreasing = TRUE)[1:10]) / sum(counts) * 100
})
bulk_qc$Top10GenesPct <- top10_prop

        round(max(bulk_qc$LibrarySize)/1e6, 2), " million reads")
        round(max(bulk_qc$ZeroInflation), 1), "%")
        round(max(bulk_qc$Top10GenesPct), 1), "%")

write.csv(bulk_qc, "bulk_RNA_QC_metrics.csv", row.names = FALSE)

library(ggplot2)
library(gridExtra)

sample_palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", 
                "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
                "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")

n_samples <- length(unique(bulk_qc$Sample))
sample_colors <- sample_palette[1:n_samples]
names(sample_colors) <- unique(bulk_qc$Sample)

theme_qc <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# A) Library Size
qc_a <- ggplot(bulk_qc, aes(x = Sample, y = LibrarySize/1e6, fill = Sample)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = sample_colors) +
  theme_qc +
  labs(title = "A) Library Size", x = "Sample ID", y = "Total Reads (millions)")

# B) Genes Detected
qc_b <- ggplot(bulk_qc, aes(x = Sample, y = GenesDetected, fill = Sample)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = sample_colors) +
  theme_qc +
  labs(title = "B) Genes Detected", x = "Sample ID", y = "Number of Genes")

# C) Zero Inflation
qc_c <- ggplot(bulk_qc, aes(x = Sample, y = ZeroInflation, fill = Sample)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = sample_colors) +
  theme_qc +
  labs(title = "C) Zero Inflation", x = "Sample ID", y = "Percentage (%)")

# D) Top 10 Genes
qc_d <- ggplot(bulk_qc, aes(x = Sample, y = Top10GenesPct, fill = Sample)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = sample_colors) +
  theme_qc +
  labs(title = "D) Top 10 Genes", x = "Sample ID", y = "Percentage (%)")

# Combine into 2x2 grid
qc_plot <- grid.arrange(qc_a, qc_b, qc_c, qc_d, ncol = 2,
                        top = grid::textGrob("Bulk RNA-seq Quality Metrics", 
                                             gp = grid::gpar(fontface = "bold", fontsize = 14)))

ggsave("FigureS1_bulk_QC.pdf", qc_plot, width = 10, height = 8, device = cairo_pdf)
ggsave("FigureS1_bulk_QC.png", qc_plot, width = 10, height = 8, dpi = 300)


## ===== Find matched samples =====
message("Matching samples between bulk and scRNA")

bulk_samples <- colnames(bulk_counts)
scrna_samples <- unique(sample_id)
matched_samples <- intersect(bulk_samples, scrna_samples)


if (length(matched_samples) == 0) {

  stop("No matching sample IDs found between bulk and scRNA!")
}

# Subset bulk to matched samples
bulk_counts <- bulk_counts[, matched_samples, drop = FALSE]

## ===== Compute scRNA ground truth proportions =====
message("Calculating experimenatl scRNA-seq  proportions")

scrna_meta <- data.frame(
  CellID = rownames(scRNA_ref@meta.data),
  CellType = cell.type,
  SampleID = sample_id
)

truth_df <- scrna_meta %>%
  filter(SampleID %in% matched_samples) %>%
  group_by(SampleID, CellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(Truth = n / sum(n)) %>%
  select(Sample = SampleID, CellType, Truth)

## ===== BayesPrism deconvolution =====
message("Running BayesPrism deconvolution")

bk.dat <- t(as.matrix(bulk_counts))  # samples x genes

# Filter and select markers
sc.dat.filt <- sc.dat[, colSums(sc.dat > 0) > 3]
diff.exp.stat <- get.exp.stat(sc.dat = sc.dat.filt, cell.type.labels = cell.type,
                               cell.state.labels = cell.state, pseudo.count = 0.1,
                               cell.count.cutoff = 10, n.cores = 1)
sc.dat.markers <- select.marker(sc.dat = sc.dat, stat = diff.exp.stat,
                                 pval.max = 0.01, lfc.min = 0.1)

# Check marker uniqueness across cell types
message("Marker uniqueness")
marker_lists <- lapply(names(diff.exp.stat), function(ct) {
  stats <- diff.exp.stat[[ct]]
  genes <- rownames(stats)[stats$pval.up.min < 0.01 & stats$min.lfc > 0.1]
  return(genes)
})
names(marker_lists) <- names(diff.exp.stat)

# Count how many cell types each marker belongs to
all_markers <- unique(unlist(marker_lists))
marker_counts <- sapply(all_markers, function(g) {
  sum(sapply(marker_lists, function(ml) g %in% ml))
})

# Summary
unique_markers <- names(marker_counts)[marker_counts == 1]
shared_markers <- names(marker_counts)[marker_counts > 1]

print(paste("  Total markers: ", length(all_markers)))
print(paste("  Unique to one cell type: ", length(unique_markers), " (", 
        round(100 * length(unique_markers) / length(all_markers), 1), "%)"))

# Per cell type uniqueness
marker_summary <- data.frame(
  CellType = names(marker_lists),
  TotalMarkers = sapply(marker_lists, length),
  UniqueMarkers = sapply(marker_lists, function(ml) sum(marker_counts[ml] == 1)),
  SharedMarkers = sapply(marker_lists, function(ml) sum(marker_counts[ml] > 1))
)
marker_summary$PctUnique <- round(100 * marker_summary$UniqueMarkers / marker_summary$TotalMarkers, 1)
print(marker_summary)
write.csv(marker_summary, "marker_uniqueness_real_bulk.csv", row.names = FALSE)

# FIX: Use full filtered reference  per BayesPrism tutorial
# The tutorial passes the full protein-coding filtered reference to new.prism(),
# Using the full reference gives the Z matrix expression estimates for all genes.
common_genes <- intersect(colnames(bk.dat), colnames(sc.dat.filt))
print(paste("Common genes : ", length(common_genes)))

bk.dat.final <- bk.dat[, common_genes]
sc.dat.final <- sc.dat.filt[, common_genes]

# Build prism and run (full reference)
myPrism <- new.prism(reference = sc.dat.final, mixture = bk.dat.final,
                      input.type = "count.matrix", cell.type.labels = cell.type,
                      cell.state.labels = cell.state, key = NULL,
                      outlier.cut = 0.01, outlier.fraction = 0.1)

bp.res <- run.prism(prism = myPrism, n.cores = 4)
save(bp.res, file = "bp.res.real_bulk_matched.RData")

# Extract estimated fractions
theta <- get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")

## ===== Compare BayesPrism vs scRNA Ground Truth =====
message("Comparing BayesPrism estimates to Experimenatal scRNA-seq")

theta_long <- as.data.frame(theta) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "BayesPrism")

# Merge with ground truth (only matched samples)
validation_df <- inner_join(theta_long, truth_df, by = c("Sample", "CellType"))
validation_df <- validation_df %>%
  mutate(Difference = BayesPrism - Truth,
         AbsDiff = abs(Difference))

overall_cor <- cor(validation_df$BayesPrism, validation_df$Truth)
overall_rho <- cor(validation_df$BayesPrism, validation_df$Truth, method = "spearman")
overall_mae <- mean(validation_df$AbsDiff)
overall_rmse <- sqrt(mean(validation_df$Difference^2))

print("Matched samples: ", paste(matched_samples, collapse = ", "))
print("Overall Pearson r: ", round(overall_cor, 3))
print("Overall Spearman rho: ", round(overall_rho, 3))
print("Overall MAE: ", round(overall_mae, 4))
print("Overall RMSE: ", round(overall_rmse, 4))

# Per cell type
ct_summary <- validation_df %>%
  group_by(CellType) %>%
  summarise(Pearson_r = cor(BayesPrism, Truth, use = "complete.obs"),
            Spearman_rho = cor(BayesPrism, Truth, method = "spearman", use = "complete.obs"),
            MAE = mean(AbsDiff),
            RMSE = sqrt(mean(Difference^2)), .groups = "drop")
print(ct_summary)

# Per sample
sample_summary <- validation_df %>%
  group_by(Sample) %>%
  summarise(Pearson_r = cor(BayesPrism, Truth),
            Spearman_rho = cor(BayesPrism, Truth, method = "spearman"),
            MAE = mean(AbsDiff), .groups = "drop")
print(sample_summary)

write.csv(validation_df, "validation_real_bulk_matched.csv", row.names = FALSE)
write.csv(ct_summary, "validation_real_bulk_matched_celltype.csv", row.names = FALSE)
write.csv(sample_summary, "validation_real_bulk_matched_sample.csv", row.names = FALSE)

## ===== mRNA Content Correction =====
message("mRNA content correction")

# Calculate mean mRNA content (total UMI) per cell type from scRNA
cell_mRNA <- rowSums(sc.dat)  # total counts per cell
mRNA_per_celltype <- tapply(cell_mRNA, cell.type, mean)
message("Mean mRNA per cell type:")
print(round(mRNA_per_celltype))

Correct BayesPrism estimates to approximate cell fractions
# Divide by mRNA content and renormalize
theta_corrected <- sweep(theta, 2, mRNA_per_celltype[colnames(theta)], "/")
theta_corrected <- theta_corrected / rowSums(theta_corrected)

# Reshape corrected estimates
theta_corrected_long <- as.data.frame(theta_corrected) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "BayesPrism_Corrected")

Create mRNA-weighted ground truth (for comparison)
scrna_meta_mRNA <- data.frame(
  CellID = names(cell.type),
  CellType = cell.type,
  SampleID = sample_id,
  mRNA = cell_mRNA
)

truth_mRNA_weighted <- scrna_meta_mRNA %>%
  filter(SampleID %in% matched_samples) %>%
  group_by(SampleID, CellType) %>%
  summarise(total_mRNA = sum(mRNA), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(Truth_mRNA = total_mRNA / sum(total_mRNA)) %>%
  select(Sample = SampleID, CellType, Truth_mRNA)

validation_full <- validation_df %>%
  left_join(theta_corrected_long, by = c("Sample", "CellType")) %>%
  left_join(truth_mRNA_weighted, by = c("Sample", "CellType")) %>%
  mutate(
    Diff_Corrected = BayesPrism_Corrected - Truth,
    Diff_vs_mRNA = BayesPrism - Truth_mRNA
  )

# Summary: Original BayesPrism vs Cell Count Truth (what you had)
ct_summary_original <- validation_full %>%
  group_by(CellType) %>%
  summarise(Pearson_r = cor(BayesPrism, Truth),
            Spearman_rho = cor(BayesPrism, Truth, method = "spearman"),
            MAE = mean(abs(BayesPrism - Truth)), .groups = "drop")
print(ct_summary_original)

# Summary: Corrected BayesPrism vs Cell Count Truth
ct_summary_corrected <- validation_full %>%
  group_by(CellType) %>%
  summarise(Pearson_r = cor(BayesPrism_Corrected, Truth),
            Spearman_rho = cor(BayesPrism_Corrected, Truth, method = "spearman"),
            MAE = mean(abs(Diff_Corrected)), .groups = "drop")
print(ct_summary_corrected)

# Summary: Original BayesPrism vs mRNA-weighted Truth
ct_summary_mRNA <- validation_full %>%
  group_by(CellType) %>%
  summarise(Pearson_r = cor(BayesPrism, Truth_mRNA),
            Spearman_rho = cor(BayesPrism, Truth_mRNA, method = "spearman"),
            MAE = mean(abs(Diff_vs_mRNA)), .groups = "drop")
print(ct_summary_mRNA)

# Overall correlations
cat("Original BayesPrism vs Cell Truth: r = ",
        round(cor(validation_full$BayesPrism, validation_full$Truth), 3),
        ", rho = ", round(cor(validation_full$BayesPrism, validation_full$Truth, method = "spearman"), 3), "\n")
cat("Corrected BayesPrism vs Cell Truth: r = ",
        round(cor(validation_full$BayesPrism_Corrected, validation_full$Truth), 3),
        ", rho = ", round(cor(validation_full$BayesPrism_Corrected, validation_full$Truth, method = "spearman"), 3), "\n")
cat("Original BayesPrism vs mRNA Truth: r = ",
        round(cor(validation_full$BayesPrism, validation_full$Truth_mRNA), 3),
        ", rho = ", round(cor(validation_full$BayesPrism, validation_full$Truth_mRNA, method = "spearman"), 3), "\n")

# Save full results
write.csv(validation_full, "validation_real_bulk_all_methods.csv", row.names = FALSE)

## ===== Visualization =====
library(ggplot2)
library(gridExtra)

cell_colors <- c(
  "B cells"         = "#4477AA",
  "T cells"         = "#EE6677",
  "Mo/Ma"           = "#228833",
  "Neutrophils"     = "#CCBB44",
  "Dendritic cells" = "#66CCEE",
  "Mast cells"      = "#AA3377"
)

theme_publication <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 11)
  )

# ---- Figure 1: mRNA Content per Cell Type (Combined with Table) ----
mRNA_df <- data.frame(
  CellType = names(mRNA_per_celltype),
  MeanUMI = as.numeric(mRNA_per_celltype)
)
mRNA_df <- mRNA_df[order(mRNA_df$MeanUMI, decreasing = TRUE), ]
mRNA_df$FoldVsMin <- round(mRNA_df$MeanUMI / min(mRNA_df$MeanUMI), 2)
mRNA_df$CellType <- factor(mRNA_df$CellType, levels = mRNA_df$CellType)

# Bar plot with value labels
fig1 <- ggplot(mRNA_df, aes(x = CellType, y = MeanUMI, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(format(round(MeanUMI), big.mark = ","), "\n(", FoldVsMin, "x)")), 
            vjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = cell_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "", y = "Mean UMI Counts per Cell")

ggsave("Figure1_mRNA_content.pdf", fig1, width = 7, height = 5, device = cairo_pdf)
ggsave("Figure1_mRNA_content.png", fig1, width = 7, height = 5, dpi = 300)

# ---- Figure 2: Three-panel Deconvolution Comparison ----
# Calculate statistics for annotations
stats_2a <- data.frame(
  r = round(cor(validation_full$BayesPrism, validation_full$Truth), 3),
  rho = round(cor(validation_full$BayesPrism, validation_full$Truth, method = "spearman"), 3),
  mse = round(mean((validation_full$BayesPrism - validation_full$Truth)^2), 4)
)
stats_2b <- data.frame(
  r = round(cor(validation_full$BayesPrism_Corrected, validation_full$Truth), 3),
  rho = round(cor(validation_full$BayesPrism_Corrected, validation_full$Truth, method = "spearman"), 3),
  mse = round(mean((validation_full$BayesPrism_Corrected - validation_full$Truth)^2), 4)
)
stats_2c <- data.frame(
  r = round(cor(validation_full$BayesPrism, validation_full$Truth_mRNA), 3),
  rho = round(cor(validation_full$BayesPrism, validation_full$Truth_mRNA, method = "spearman"), 3),
  mse = round(mean((validation_full$BayesPrism - validation_full$Truth_mRNA)^2), 4)
)

# Plot 2A: Original BayesPrism vs Cell Count Truth
p2a <- ggplot(validation_full, aes(x = Truth, y = BayesPrism, color = CellType)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  scale_color_manual(values = cell_colors) +
  annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0, size = 3.5,
           label = paste0("r = ", stats_2a$r, "\nρ = ", stats_2a$rho, "\nMSE = ", stats_2a$mse)) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 6, fontface = "bold", label = "A") +
  theme_publication +
  labs(x = "Experimental scRNA-seq (Cell Counts)", y = "Raw BayesPrism Estimate") +
  xlim(0, 1) + ylim(0, 1) +
  theme(legend.position = "none")

# Plot 2B: Corrected BayesPrism vs Cell Count Truth
p2b <- ggplot(validation_full, aes(x = Truth, y = BayesPrism_Corrected, color = CellType)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  scale_color_manual(values = cell_colors) +
  annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0, size = 3.5,
           label = paste0("r = ", stats_2b$r, "\nρ = ", stats_2b$rho, "\nMSE = ", stats_2b$mse)) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 6, fontface = "bold", label = "B") +
  theme_publication +
  labs(x = "Experimental scRNA-seq (Cell Counts)", y = "mRNA-corrected BayesPrism Estimate") +
  xlim(0, 1) + ylim(0, 1) +
  theme(legend.position = "none")

# Plot 2C: Original BayesPrism vs mRNA-weighted Truth
p2c <- ggplot(validation_full, aes(x = Truth_mRNA, y = BayesPrism, color = CellType)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  scale_color_manual(values = cell_colors) +
  annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0, size = 3.5,
           label = paste0("r = ", stats_2c$r, "\nρ = ", stats_2c$rho, "\nMSE = ", stats_2c$mse)) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 6, fontface = "bold", label = "C") +
  theme_publication +
  labs(x = "Experimental scRNA-seq (mRNA-weighted)", y = "Raw BayesPrism Estimate") +
  xlim(0, 1) + ylim(0, 1) +
  theme(legend.position = "none")

# Create shared legend
legend_plot <- ggplot(validation_full, aes(x = Truth, y = BayesPrism, color = CellType)) +
  geom_point(size = 4) +
  scale_color_manual(values = cell_colors, name = "Cell Type") +
  theme_publication +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"))
legend <- cowplot::get_legend(legend_plot)

# Combine plots with legend using cowplot
library(cowplot)
fig2_plots <- plot_grid(p2a, p2b, p2c, ncol = 3, align = "h")
fig2_combined <- plot_grid(fig2_plots, legend, ncol = 1, rel_heights = c(1, 0.08))

ggsave("Figure2_deconvolution_comparison.pdf", fig2_combined, width = 12, height = 5.5, device = cairo_pdf)
ggsave("Figure2_deconvolution_comparison.png", fig2_combined, width = 12, height = 5.5, dpi = 300)

# ---- Figure 3: Per-Cell Type Scatter Plots (Faceted) ----
# Original version (without condition split)
fig3 <- ggplot(validation_full, aes(x = Truth_mRNA, y = BayesPrism)) +
  geom_point(aes(color = CellType), size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  facet_wrap(~CellType, scales = "free", ncol = 3) +
  scale_color_manual(values = cell_colors) +
  theme_publication +
  theme(legend.position = "none") +
  labs(title = "BayesPrism Estimates vs mRNA-weighted Experimental scRNA-seq by Cell Type",
       x = "Experimental scRNA-seq (mRNA-weighted)", y = "Raw BayesPrism Estimate")

ggsave("Figure3_per_celltype_scatter.pdf", fig3, width = 10, height = 7, device = cairo_pdf)
ggsave("Figure3_per_celltype_scatter.png", fig3, width = 10, height = 7, dpi = 300)

# ---- Figure 3B: Per-Cell Type with Condition as Color/Shape (Option C) ----
message("Figure 3B: Per-cell type scatter with condition")

# Add condition to validation_full (needed before section 9)
validation_full_cond <- validation_full %>%
  left_join(
    scRNA_ref@meta.data %>%
      select(sample_id, disease_state) %>%
      distinct() %>%
      mutate(Sample = as.character(sample_id), Condition = disease_state) %>%
      select(Sample, Condition),
    by = "Sample"
  )

# Calculate per-cell-type, per-condition correlations for annotation
cor_annotations <- validation_full_cond %>%
  group_by(CellType, Condition) %>%
  summarise(
    r = round(cor(BayesPrism, Truth_mRNA, use = "complete.obs"), 2),
    rho = round(cor(BayesPrism, Truth_mRNA, method = "spearman", use = "complete.obs"), 2),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Condition, values_from = c(r, rho), names_sep = "_")

# Create label for each cell type
cor_labels <- validation_full_cond %>%
  group_by(CellType) %>%
  summarise(
    x_pos = max(Truth_mRNA, na.rm = TRUE) * 0.95,
    y_pos = min(BayesPrism, na.rm = TRUE) + (max(BayesPrism, na.rm = TRUE) - min(BayesPrism, na.rm = TRUE)) * 0.15,
    .groups = "drop"
  ) %>%
  left_join(cor_annotations, by = "CellType") %>%
  mutate(label = paste0("CTL: r=", r_CTL, ", \u03C1=", rho_CTL, "\nSEA: r=", r_SEA, ", \u03C1=", rho_SEA))

# Define condition colors and shapes
condition_colors <- c("CTL" = "#4DBBD5", "SEA" = "#E64B35")
condition_shapes <- c("CTL" = 16, "SEA" = 17)  # circle and triangle

# Create Figure 3B
fig3b <- ggplot(validation_full_cond, aes(x = Truth_mRNA, y = BayesPrism)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  geom_point(aes(color = Condition, shape = Condition), size = 3, alpha = 0.8) +
  geom_smooth(aes(color = Condition, linetype = Condition), method = "lm", se = FALSE, linewidth = 0.8) +
  geom_text(data = cor_labels, aes(x = x_pos, y = y_pos, label = label),
            hjust = 1, vjust = 0, size = 3, fontface = "bold") +
  facet_wrap(~CellType, scales = "free", ncol = 3) +
  scale_color_manual(values = condition_colors, name = "Condition") +
  scale_shape_manual(values = condition_shapes, name = "Condition") +
  scale_linetype_manual(values = c("CTL" = "solid", "SEA" = "solid"), guide = "none") +
  theme_publication +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "BayesPrism vs mRNA-weighted Experimental scRNA-seq by Cell Type and Condition",
    x = "Experimental scRNA-seq (mRNA-weighted)",
    y = "Raw BayesPrism Estimate"
  )

ggsave("Figure3B_per_celltype_by_condition.pdf", fig3b, width = 10, height = 7, device = cairo_pdf)
ggsave("Figure3B_per_celltype_by_condition.png", fig3b, width = 10, height = 7, dpi = 300)

# ---- Figure 4: Cell Type Composition Comparison (Stacked Bar) ----
# Prepare long-format data for composition plots
composition_df <- validation_full %>%
  select(Sample, CellType, Truth, BayesPrism, BayesPrism_Corrected, Truth_mRNA) %>%
  pivot_longer(cols = c(Truth, BayesPrism, BayesPrism_Corrected, Truth_mRNA),
               names_to = "Method", values_to = "Proportion") %>%
  mutate(Method = factor(Method, levels = c("Truth", "BayesPrism", "BayesPrism_Corrected", "Truth_mRNA"),
                         labels = c("Cell Count Truth", "BayesPrism Original", 
                                   "BayesPrism Corrected", "mRNA-weighted Truth")))

# Overall mean composition across samples for each method (main Figure 4)
composition_summary <- composition_df %>%
  group_by(Method, CellType) %>%
  summarise(Proportion = mean(Proportion), .groups = "drop")

fig4 <- ggplot(composition_summary, aes(x = Method, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom") +
  labs(title = "Mean Cell Type Composition Across Samples",
       x = "Method", y = "Mean Proportion")

ggsave("Figure4_composition_barplot.pdf", fig4, width = 7, height = 5, device = cairo_pdf)
ggsave("Figure4_composition_barplot.png", fig4, width = 7, height = 5, dpi = 300)

# ---- Figure S2: Cell Count vs mRNA-weighted Composition (Stacked Bar) ----
truth_composition_df <- composition_df %>%
  filter(Method %in% c("Cell Count Truth", "mRNA-weighted Truth"))

figS2 <- ggplot(truth_composition_df, aes(x = Method, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  facet_wrap(~Sample, ncol = 1) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom") +
  labs(title = "Cell Count vs mRNA-weighted Composition per Sample",
       x = "", y = "Proportion")

ggsave("FigureS2_truth_vs_mRNA_composition.pdf", figS2, width = 7, height = 18, device = cairo_pdf)
ggsave("FigureS2_truth_vs_mRNA_composition.png", figS2, width = 7, height = 18, dpi = 300)

# ---- Figure S3: Full Composition per Sample Across Methods (Stacked Bar) ----
figS3 <- ggplot(composition_df, aes(x = Method, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  facet_wrap(~Sample, ncol = 1) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom") +
  labs(title = "Cell Type Composition per Sample Across Methods",
       x = "Method", y = "Proportion")

ggsave("FigureS3_composition_per_sample.pdf", figS3, width = 7, height = 18, device = cairo_pdf)
ggsave("FigureS3_composition_per_sample.png", figS3, width = 7, height = 18, dpi = 300)

# ---- Export Summary Tables ----

# cell count and samples. 
# Total cells
total_cells <- ncol(scRNA_ref)
print("Total cells: ", total_cells)

# Total samples (unique sample IDs)
total_samples <- length(unique(scRNA_ref@meta.data$sample_id))
print("Total samples: ", total_samples)

# Or combined
print("scRNA-seq: ", total_cells, " cells across ", total_samples, " samples")

# Table 1: Marker genes per cell type
write.csv(marker_summary, "Table1_marker_summary.csv", row.names = FALSE)

# Table 2: mRNA content per cell type
mRNA_table <- data.frame(
  CellType = names(mRNA_per_celltype),
  MeanUMI = round(as.numeric(mRNA_per_celltype), 0),
  RelativeToMin = round(as.numeric(mRNA_per_celltype) / min(mRNA_per_celltype), 2)
)
mRNA_table <- mRNA_table[order(mRNA_table$MeanUMI, decreasing = TRUE), ]
write.csv(mRNA_table, "Table2_mRNA_content.csv", row.names = FALSE)

# Table 3: Deconvolution performance summary (all three comparisons)
performance_table <- data.frame(
  CellType = ct_summary_original$CellType,
  Original_r = round(ct_summary_original$Pearson_r, 3),
  Original_rho = round(ct_summary_original$Spearman_rho, 3),
  Original_MAE = round(ct_summary_original$MAE, 4),
  Corrected_r = round(ct_summary_corrected$Pearson_r, 3),
  Corrected_rho = round(ct_summary_corrected$Spearman_rho, 3),
  Corrected_MAE = round(ct_summary_corrected$MAE, 4),
  mRNATruth_r = round(ct_summary_mRNA$Pearson_r, 3),
  mRNATruth_rho = round(ct_summary_mRNA$Spearman_rho, 3),
  mRNATruth_MAE = round(ct_summary_mRNA$MAE, 4)
)
# Add overall row
performance_table <- rbind(performance_table, data.frame(
  CellType = "OVERALL",
  Original_r = round(cor(validation_full$BayesPrism, validation_full$Truth), 3),
  Original_rho = round(cor(validation_full$BayesPrism, validation_full$Truth, method = "spearman"), 3),
  Original_MAE = round(mean(abs(validation_full$BayesPrism - validation_full$Truth)), 4),
  Corrected_r = round(cor(validation_full$BayesPrism_Corrected, validation_full$Truth), 3),
  Corrected_rho = round(cor(validation_full$BayesPrism_Corrected, validation_full$Truth, method = "spearman"), 3),
  Corrected_MAE = round(mean(abs(validation_full$Diff_Corrected)), 4),
  mRNATruth_r = round(cor(validation_full$BayesPrism, validation_full$Truth_mRNA), 3),
  mRNATruth_rho = round(cor(validation_full$BayesPrism, validation_full$Truth_mRNA, method = "spearman"), 3),
  mRNATruth_MAE = round(mean(abs(validation_full$Diff_vs_mRNA)), 4)
))
write.csv(performance_table, "Table2_deconvolution_performance.csv", row.names = FALSE)

## ===== Condition-Specific Analysis (SEA vs CTL) =====
message("Condition-specific analysis (SEA vs CTL)")

# Get condition for each sample from scRNA metadata
sample_condition <- scRNA_ref@meta.data %>%
  select(sample_id, disease_state) %>%
  distinct() %>%
  mutate(Sample = as.character(sample_id),
         Condition = disease_state) %>%
  select(Sample, Condition)

# Add condition to validation data
validation_by_condition <- validation_full %>%
  left_join(sample_condition, by = "Sample")

# Condition table
print("Samples per condition:")
print(table(unique(validation_by_condition[, c("Sample", "Condition")])$Condition))

# ---- Calculate statistics by condition ----

# Function to calculate stats for a subset
calc_stats <- function(df, bp_col, truth_col) {
  data.frame(
    r = round(cor(df[[bp_col]], df[[truth_col]], use = "complete.obs"), 3),
    rho = round(cor(df[[bp_col]], df[[truth_col]], method = "spearman", use = "complete.obs"), 3),
    MAE = round(mean(abs(df[[bp_col]] - df[[truth_col]]), na.rm = TRUE), 4)
  )
}

# Overall stats by condition
condition_stats <- validation_by_condition %>%
  group_by(Condition) %>%
  summarise(
    n_samples = n_distinct(Sample),
    # Original BayesPrism vs Cell Count Truth
    Original_vs_CellCount_r = round(cor(BayesPrism, Truth), 3),
    Original_vs_CellCount_rho = round(cor(BayesPrism, Truth, method = "spearman"), 3),
    # Corrected BayesPrism vs Cell Count Truth
    Corrected_vs_CellCount_r = round(cor(BayesPrism_Corrected, Truth), 3),
    Corrected_vs_CellCount_rho = round(cor(BayesPrism_Corrected, Truth, method = "spearman"), 3),
    # Original BayesPrism vs mRNA-weighted Truth
    Original_vs_mRNA_r = round(cor(BayesPrism, Truth_mRNA), 3),
    Original_vs_mRNA_rho = round(cor(BayesPrism, Truth_mRNA, method = "spearman"), 3),
    .groups = "drop"
  )
print("\nOverall correlation by condition:")
print(condition_stats)

# Per cell type stats by condition
celltype_condition_stats <- validation_by_condition %>%
  group_by(Condition, CellType) %>%
  summarise(
    n = n(),
    Original_vs_CellCount_r = round(cor(BayesPrism, Truth, use = "complete.obs"), 3),
    Original_vs_CellCount_rho = round(cor(BayesPrism, Truth, method = "spearman", use = "complete.obs"), 3),
    Corrected_vs_CellCount_r = round(cor(BayesPrism_Corrected, Truth, use = "complete.obs"), 3),
    Corrected_vs_CellCount_rho = round(cor(BayesPrism_Corrected, Truth, method = "spearman", use = "complete.obs"), 3),
    Original_vs_mRNA_r = round(cor(BayesPrism, Truth_mRNA, use = "complete.obs"), 3),
    Original_vs_mRNA_rho = round(cor(BayesPrism, Truth_mRNA, method = "spearman", use = "complete.obs"), 3),
    Original_vs_CellCount_MAE = round(mean(abs(BayesPrism - Truth)), 4),
    Corrected_vs_CellCount_MAE = round(mean(abs(BayesPrism_Corrected - Truth)), 4),
    .groups = "drop"
  )

print("\nPer cell type correlation by condition:")
print(celltype_condition_stats %>% arrange(CellType, Condition))

# Save condition-specific stats
write.csv(condition_stats, "Table_condition_overall_stats.csv", row.names = FALSE)
write.csv(celltype_condition_stats, "Table_condition_celltype_stats.csv", row.names = FALSE)

# ---- Figure: Correlation plots by condition ----
message("\n  Creating condition-specific correlation plots")

# Define condition colors
condition_colors <- c("CTL" = "#4DBBD5", "SEA" = "#E64B35")

# Function to create scatter plot for a condition
create_condition_plot <- function(data, condition_name, x_var, y_var, x_label, y_label, panel_label) {
  df_subset <- data %>% filter(Condition == condition_name)

  r_val <- round(cor(df_subset[[x_var]], df_subset[[y_var]], use = "complete.obs"), 3)
  rho_val <- round(cor(df_subset[[x_var]], df_subset[[y_var]], method = "spearman", use = "complete.obs"), 3)
  n_samples <- n_distinct(df_subset$Sample)

  ggplot(df_subset, aes_string(x = x_var, y = y_var, color = "CellType")) +
    geom_point(size = 3, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
    scale_color_manual(values = cell_colors) +
    annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0, size = 3.5,
             label = paste0("r = ", r_val, "\n\u03C1 = ", rho_val, "\nn = ", n_samples, " samples")) +
    annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 5, fontface = "bold",
             label = panel_label) +
    theme_publication +
    labs(title = condition_name, x = x_label, y = y_label) +
    xlim(0, 1) + ylim(0, 1) +
    theme(legend.position = "none")
}

# ---- Figure 2B: Original BayesPrism vs Cell Count Truth by Condition ----
p_ctl_original <- create_condition_plot(
  validation_by_condition, "CTL", "Truth", "BayesPrism",
  "Ground Truth (Cell Counts)", "BayesPrism Estimate", "A"
)
p_sea_original <- create_condition_plot(
  validation_by_condition, "SEA", "Truth", "BayesPrism",
  "Ground Truth (Cell Counts)", "BayesPrism Estimate", "B"
)

# ---- Figure 2C: Corrected BayesPrism vs Cell Count Truth by Condition ----
p_ctl_corrected <- create_condition_plot(
  validation_by_condition, "CTL", "Truth", "BayesPrism_Corrected",
  "Ground Truth (Cell Counts)", "BayesPrism Corrected", "C"
)
p_sea_corrected <- create_condition_plot(
  validation_by_condition, "SEA", "Truth", "BayesPrism_Corrected",
  "Ground Truth (Cell Counts)", "BayesPrism Corrected", "D"
)

# ---- Figure 2D: Original BayesPrism vs mRNA Truth by Condition ----
p_ctl_mRNA <- create_condition_plot(
  validation_by_condition, "CTL", "Truth_mRNA", "BayesPrism",
  "Ground Truth (mRNA-weighted)", "BayesPrism Estimate", "E"
)
p_sea_mRNA <- create_condition_plot(
  validation_by_condition, "SEA", "Truth_mRNA", "BayesPrism",
  "Ground Truth (mRNA-weighted)", "BayesPrism Estimate", "F"
)

# Create shared legend
legend_plot_cond <- ggplot(validation_by_condition, aes(x = Truth, y = BayesPrism, color = CellType)) +
  geom_point(size = 4) +
  scale_color_manual(values = cell_colors, name = "Cell Type") +
  theme_publication +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"))
legend_cond <- cowplot::get_legend(legend_plot_cond)

# Combine: 3 rows (Original, Corrected, mRNA) x 2 columns (CTL, SEA)
fig_by_condition <- plot_grid(
  # Row labels
  plot_grid(
    ggdraw() + draw_label("Original vs\nCell Counts", fontface = "bold", size = 10),
    p_ctl_original, p_sea_original, ncol = 3, rel_widths = c(0.15, 1, 1)
  ),
  plot_grid(
    ggdraw() + draw_label("Corrected vs\nCell Counts", fontface = "bold", size = 10),
    p_ctl_corrected, p_sea_corrected, ncol = 3, rel_widths = c(0.15, 1, 1)
  ),
  plot_grid(
    ggdraw() + draw_label("Original vs\nmRNA Truth", fontface = "bold", size = 10),
    p_ctl_mRNA, p_sea_mRNA, ncol = 3, rel_widths = c(0.15, 1, 1)
  ),
  ncol = 1, rel_heights = c(1, 1, 1)
)

fig_by_condition_final <- plot_grid(fig_by_condition, legend_cond, ncol = 1, rel_heights = c(1, 0.05))

ggsave("Figure_correlation_by_condition.pdf", fig_by_condition_final, width = 10, height = 12, device = cairo_pdf)
ggsave("Figure_correlation_by_condition.png", fig_by_condition_final, width = 10, height = 12, dpi = 300)


# ---- Simpler 2-panel figure: CTL vs SEA for mRNA-weighted comparison ----
fig_mRNA_by_condition <- plot_grid(p_ctl_mRNA, p_sea_mRNA, ncol = 2, align = "h")
fig_mRNA_by_condition_final <- plot_grid(fig_mRNA_by_condition, legend_cond, ncol = 1, rel_heights = c(1, 0.1))

ggsave("Figure_mRNA_correlation_CTL_vs_SEA.pdf", fig_mRNA_by_condition_final, width = 10, height = 5.5, device = cairo_pdf)
ggsave("Figure_mRNA_correlation_CTL_vs_SEA.png", fig_mRNA_by_condition_final, width = 10, height = 5.5, dpi = 300)


# ---- Summary Table: Condition-specific performance ----
condition_summary_table <- validation_by_condition %>%
  group_by(Condition) %>%
  summarise(
    n_samples = n_distinct(Sample),
    `BayesPrism vs Cell Count (r)` = round(cor(BayesPrism, Truth), 3),
    `BayesPrism vs Cell Count (rho)` = round(cor(BayesPrism, Truth, method = "spearman"), 3),
    `Corrected vs Cell Count (r)` = round(cor(BayesPrism_Corrected, Truth), 3),
    `Corrected vs Cell Count (rho)` = round(cor(BayesPrism_Corrected, Truth, method = "spearman"), 3),
    `BayesPrism vs mRNA Truth (r)` = round(cor(BayesPrism, Truth_mRNA), 3),
    `BayesPrism vs mRNA Truth (rho)` = round(cor(BayesPrism, Truth_mRNA, method = "spearman"), 3),
    `BayesPrism vs Cell Count (MAE)` = round(mean(abs(BayesPrism - Truth)), 4),
    `Corrected vs Cell Count (MAE)` = round(mean(abs(BayesPrism_Corrected - Truth)), 4),
    .groups = "drop"
  )

print(condition_summary_table)

write.csv(condition_summary_table, "Table_condition_performance_summary.csv", row.names = FALSE)

