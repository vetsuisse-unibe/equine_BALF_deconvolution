############################################################
# Comparison with a second deconvolution method (CIBERSORTx)
# Deconvolve the same bulk mixtures with CIBERSORTx and compare against
# cell-count and mRNA-weighted scRNA-seq proportions, to test whether the
# mRNA-content bias is BayesPrism-specific or a general property.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

data_dir    <- "data"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

## ===== Load inputs =====
cs  <- read.csv(file.path(data_dir, "CIBERSORTx_Job3_Results.csv"), check.names = FALSE)
val <- read.csv(file.path(data_dir, "validation_real_bulk_all_methods.csv"), check.names = FALSE)

celltypes <- c("Mo/Ma","T cells","Neutrophils","Mast cells","B cells","Dendritic cells")

# CIBERSORTx -> long; strip "A" prefix to match scRNA sample IDs ("A30" -> "30")
cs_long <- cs %>%
  select(Mixture, all_of(celltypes)) %>%
  pivot_longer(-Mixture, names_to = "CellType", values_to = "CIBERSORTx") %>%
  mutate(Sample = sub("^A", "", Mixture))

# Keep the 11 matched samples (those that also have an scRNA-seq reference)
m <- inner_join(val %>% mutate(Sample = as.character(Sample)),
                cs_long %>% select(Sample, CellType, CIBERSORTx),
                by = c("Sample", "CellType"))

## ===== Agreement metrics =====
score <- function(est, obs) data.frame(
  Pearson_r    = cor(est, obs),
  Spearman_rho = cor(est, obs, method = "spearman"),
  MAE          = mean(abs(est - obs))
)

metrics <- bind_rows(
  cbind(Method = "BayesPrism",           Reference = "cell_counts",   score(m$BayesPrism,           m$Truth)),
  cbind(Method = "BayesPrism_corrected", Reference = "cell_counts",   score(m$BayesPrism_Corrected, m$Truth)),
  cbind(Method = "CIBERSORTx",           Reference = "cell_counts",   score(m$CIBERSORTx,           m$Truth)),
  cbind(Method = "BayesPrism",           Reference = "mRNA_weighted", score(m$BayesPrism,           m$Truth_mRNA)),
  cbind(Method = "CIBERSORTx",           Reference = "mRNA_weighted", score(m$CIBERSORTx,           m$Truth_mRNA))
)
write.csv(metrics, file.path(results_dir, "cibersortx_metrics.csv"), row.names = FALSE)

means <- m %>%
  group_by(CellType) %>%
  summarise(CellCounts = mean(Truth), BayesPrism = mean(BayesPrism),
            CIBERSORTx = mean(CIBERSORTx), .groups = "drop")
write.csv(means, file.path(results_dir, "cibersortx_celltype_means.csv"), row.names = FALSE)

print(metrics)
print(as.data.frame(means))

## ===== Figure =====
# Cell-type palette matched to Figure 3
cell_colors <- c(
  "B cells"         = "#4477AA",  # blue
  "T cells"         = "#EE6677",  # coral
  "Mo/Ma"           = "#228833",  # green
  "Neutrophils"     = "#CCBB44",  # gold
  "Dendritic cells" = "#66CCEE",  # cyan
  "Mast cells"      = "#AA3377"   # magenta
)
method_colors <- c("Cell counts" = "gray70", "BayesPrism" = "#228833", "CIBERSORTx" = "#E64B35")

theme_publication <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    axis.title       = element_text(face = "bold", size = 11),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position  = "none"
  )

# r / Spearman rho / MAE label for the scatter panels
metric_label <- function(est, obs) {
  sprintf("r = %.2f\nρ = %.2f\nMAE = %.3f",
          cor(est, obs), cor(est, obs, method = "spearman"), mean(abs(est - obs)))
}

# ---- Panels A and B: CIBERSORTx vs the two scRNA-seq references ----
pA <- ggplot(m, aes(x = Truth, y = CIBERSORTx, color = CellType)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = cell_colors) +
  xlim(0, 1) + ylim(0, 1) +
  annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0, size = 3.4,
           label = metric_label(m$CIBERSORTx, m$Truth)) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 6, fontface = "bold", label = "A") +
  labs(x = "Experimental scRNA-seq (cell counts)", y = "CIBERSORTx estimate") +
  theme_publication

pB <- ggplot(m, aes(x = Truth_mRNA, y = CIBERSORTx, color = CellType)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = cell_colors) +
  xlim(0, 1) + ylim(0, 1) +
  annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0, size = 3.4,
           label = metric_label(m$CIBERSORTx, m$Truth_mRNA)) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 6, fontface = "bold", label = "B") +
  labs(x = "Experimental scRNA-seq (mRNA-weighted)", y = "CIBERSORTx estimate") +
  theme_publication

# ---- Panel C: mean proportion per cell type, by source ----
bar <- means %>%
  pivot_longer(-CellType, names_to = "Source", values_to = "Prop") %>%
  mutate(
    Source   = factor(Source, levels = c("CellCounts", "BayesPrism", "CIBERSORTx"),
                       labels = names(method_colors)),
    CellType = factor(CellType, levels = celltypes)
  )

pC <- ggplot(bar, aes(x = CellType, y = Prop, fill = Source)) +
  geom_col(position = position_dodge(0.8), width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = method_colors) +
  annotate("text", x = 0.6, y = 0.96, hjust = 0, vjust = 1, size = 6, fontface = "bold", label = "C") +
  labs(x = NULL, y = "Mean proportion") +
  theme_publication +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

# Shared cell-type legend along the bottom of panels A and B (matched to Figure 3)
cell_legend <- get_legend(
  pA + labs(color = "Cell type") +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 4))) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 10),
          legend.text  = element_text(size = 9))
)

scatter_row <- plot_grid(pA, pB, ncol = 2, align = "h")
top <- plot_grid(scatter_row, cell_legend, ncol = 1, rel_heights = c(1, 0.1))
fig <- plot_grid(top, pC, ncol = 1, rel_heights = c(1, 1.05))

ggsave(file.path(results_dir, "Figure_CIBERSORTx_comparison.pdf"), fig,
       width = 10, height = 9, device = cairo_pdf)
ggsave(file.path(results_dir, "Figure_CIBERSORTx_comparison.png"), fig,
       width = 10, height = 9, dpi = 300)

message("Wrote metrics, cell-type means, and Figure_CIBERSORTx_comparison to ", results_dir)
