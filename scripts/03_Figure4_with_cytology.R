############################################################
# Figure 4: Cell Type Composition with Cytology Data
# Replaces mRNA-weighted truth with cytology counts
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggplot2)
})

## ===== Setup =====
data_dir <- "data"

## ===== Load Cytology Data =====
message("Loading cytology data")

cyto_df <- read_excel(file.path(data_dir, "DCC_cyto.xlsx"))

# Filter to cytopre (pre-challenge cytology)
cyto_pre <- cyto_df %>%
  filter(type == "cytopre") %>%
  select(horse, status, macrophages, lymphocytes, neutrophils, mastocytes, eosinophils) %>%
  mutate(across(c(macrophages, lymphocytes, neutrophils, mastocytes, eosinophils), as.numeric))

# Convert to proportions (they appear to be percentages)
cyto_pre <- cyto_pre %>%
  mutate(
    total = macrophages + lymphocytes + neutrophils + mastocytes + eosinophils,
    macrophages = macrophages / total,
    lymphocytes = lymphocytes / total,
    neutrophils = neutrophils / total,
    mastocytes = mastocytes / total,
    eosinophils = eosinophils / total
  ) %>%
  select(-total)

cyto_pre$Sample <- sub("^A", "", cyto_pre$horse)

message("  Loaded cytopre data for ", nrow(cyto_pre), " samples")

# Filter to cytopost (post-challenge cytology)
cyto_post <- cyto_df %>%
  filter(type == "cytopost") %>%
  select(horse, status, macrophages, lymphocytes, neutrophils, mastocytes, eosinophils) %>%
  mutate(across(c(macrophages, lymphocytes, neutrophils, mastocytes, eosinophils), as.numeric))

# Convert to proportions
cyto_post <- cyto_post %>%
  mutate(
    total = macrophages + lymphocytes + neutrophils + mastocytes + eosinophils,
    macrophages = macrophages / total,
    lymphocytes = lymphocytes / total,
    neutrophils = neutrophils / total,
    mastocytes = mastocytes / total,
    eosinophils = eosinophils / total
  ) %>%
  select(-total)

cyto_post$Sample <- sub("^A", "", cyto_post$horse)

message("  Loaded cytopost data for ", nrow(cyto_post), " samples")

## ===== Load BayesPrism Results =====
message("Loading BayesPrism results")

load("bp.res.real_bulk_matched.RData")

# Get cell type fractions
library(BayesPrism)
theta <- get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")

# Convert to dataframe
bp_df <- as.data.frame(theta)
bp_df$Sample <- rownames(bp_df)

## ===== Load scRNA Ground Truth =====
message("Loading scRNA ground truth")

library(Seurat)
scRNA_ref <- readRDS(file.path(data_dir, "sc22_all_seed.rds"))

# Calculate cell count proportions per sample
sample_id <- as.character(scRNA_ref@meta.data$sample_id)
cell_type <- as.character(scRNA_ref@meta.data$MajCellType)

ground_truth <- table(sample_id, cell_type)
ground_truth_prop <- prop.table(ground_truth, margin = 1)
truth_df <- as.data.frame.matrix(ground_truth_prop)
truth_df$Sample <- rownames(truth_df)

## ===== Calculate mRNA-corrected BayesPrism =====
message("Calculating mRNA-corrected estimates")

# Get mRNA content per cell type
sc.dat <- t(as.matrix(GetAssayData(scRNA_ref, assay = "RNA", layer = "counts")))
cell_mRNA <- rowSums(sc.dat)
mRNA_per_celltype <- tapply(cell_mRNA, cell_type, mean)

# Correct BayesPrism estimates
theta_corrected <- theta
for (ct in colnames(theta)) {
  if (ct %in% names(mRNA_per_celltype)) {
    theta_corrected[, ct] <- theta[, ct] / mRNA_per_celltype[ct]
  }
}
theta_corrected <- theta_corrected / rowSums(theta_corrected)

corrected_df <- as.data.frame(theta_corrected)
corrected_df$Sample <- rownames(corrected_df)

## ===== Find Matched Samples =====
message("Finding matched samples")

# Find samples with all data types
matched_samples <- Reduce(intersect, list(
  truth_df$Sample,
  bp_df$Sample,
  cyto_pre$Sample,
  cyto_post$Sample
))

print(paste0("  Matched samples: ", paste(matched_samples, collapse = ", ")))

## ===== Prepare Data for Plotting =====
message("Preparing data for plotting")

# Cell type color palette (updated - merged categories)
cell_colors <- c(
  "Lymphocytes"     = "#EE6677",  # Red/coral (T cells dominant)
  "Mo/Ma/DC"        = "#228833",  # Green (Mo/Ma dominant)
  "Neutrophils"     = "#CCBB44",  # Yellow/gold
  "Mast cells"      = "#AA3377"   # Purple/magenta
)

# Function to combine T cells + B cells into Lymphocytes AND Mo/Ma + DC into Mo/Ma/DC
combine_cell_types <- function(df, sample_col = "Sample") {
  df %>%
    mutate(
      Lymphocytes = `T cells` + `B cells`,
      `Mo/Ma/DC` = `Mo/Ma` + `Dendritic cells`
    ) %>%
    select(-`T cells`, -`B cells`, -`Mo/Ma`, -`Dendritic cells`)
}

# Prepare scRNA ground truth - combine T+B into Lymphocytes and Mo/Ma+DC into Mo/Ma/DC
truth_combined <- truth_df %>%
  filter(Sample %in% matched_samples) %>%
  combine_cell_types()

truth_long <- truth_combined %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "Cell Count (scRNA)")

# Prepare BayesPrism original - combine cell types
bp_combined <- bp_df %>%
  filter(Sample %in% matched_samples) %>%
  combine_cell_types()

bp_long <- bp_combined %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "BayesPrism Original")

# Prepare BayesPrism corrected - combine cell types
corrected_combined <- corrected_df %>%
  filter(Sample %in% matched_samples) %>%
  combine_cell_types()

corrected_long <- corrected_combined %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "BayesPrism Corrected")

# Prepare cytology PRE data (long format)
# Note: Cytology macrophages = Mo/Ma/DC (they don't distinguish DC)
cyto_pre_long <- cyto_pre %>%
  filter(Sample %in% matched_samples) %>%
  select(Sample, macrophages, lymphocytes, neutrophils, mastocytes) %>%
  rename(
    "Mo/Ma/DC" = macrophages,
    "Lymphocytes" = lymphocytes,
    "Neutrophils" = neutrophils,
    "Mast cells" = mastocytes
  ) %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "Cytology Pre")

# Prepare cytology POST data (long format)
cyto_post_long <- cyto_post %>%
  filter(Sample %in% matched_samples) %>%
  select(Sample, macrophages, lymphocytes, neutrophils, mastocytes) %>%
  rename(
    "Mo/Ma/DC" = macrophages,
    "Lymphocytes" = lymphocytes,
    "Neutrophils" = neutrophils,
    "Mast cells" = mastocytes
  ) %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "Cytology Post")

# Combine all data
all_data <- bind_rows(truth_long, cyto_pre_long, cyto_post_long, bp_long, corrected_long)

# Set method order
all_data$Method <- factor(all_data$Method,
                          levels = c("Cell Count (scRNA)", "Cytology Pre", "Cytology Post",
                                    "BayesPrism Original", "BayesPrism Corrected"))

# Set cell type order (merged categories)
all_data$CellType <- factor(all_data$CellType,
                            levels = c("Mo/Ma/DC", "Lymphocytes", "Neutrophils",
                                      "Mast cells"))

## ===== Add Condition Information =====
message("Adding condition information")

# Get condition for each sample from scRNA-seq metadata
sample_condition <- scRNA_ref@meta.data %>%
  select(sample_id, disease_state) %>%
  distinct() %>%
  mutate(Sample = as.character(sample_id),
         Condition = disease_state) %>%
  select(Sample, Condition)

# Add condition to all_data
all_data <- all_data %>%
  left_join(sample_condition, by = "Sample")

print("  Samples per condition:")
print(table(unique(all_data[, c("Sample", "Condition")])$Condition))

## ===== Calculate Mean Composition =====
message("Calculating mean composition")

# Overall composition (for supplementary figure)
composition_summary <- all_data %>%
  group_by(Method, CellType) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

print(composition_summary %>% pivot_wider(names_from = CellType, values_from = Proportion))

# Composition by condition (for main Figure 4)
composition_by_condition <- all_data %>%
  group_by(Condition, Method, CellType) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

message("\nComposition by condition:")
print(composition_by_condition %>%
        pivot_wider(names_from = CellType, values_from = Proportion) %>%
        arrange(Condition, Method))

## ===== Create Figures =====
message("Creating figures")

theme_publication <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# ---- Main Figure 4: Composition by Condition (CTL vs SEA) ----
fig4_by_condition <- ggplot(composition_by_condition, aes(x = Method, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  facet_wrap(~Condition, ncol = 2) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_publication +
  labs(
    x = "Method",
    y = "Mean Proportion"
  )

dir.create("manuscript", showWarnings = FALSE)
ggsave("manuscript/Figure4_composition_by_condition.pdf", fig4_by_condition, width = 12, height = 6, device = cairo_pdf)
ggsave("manuscript/Figure4_composition_by_condition.png", fig4_by_condition, width = 12, height = 6, dpi = 300)
message("  Saved: Figure4_composition_by_condition.pdf/png (Main Figure)")

# ---- Supplementary Figure: Overall Composition (Figure 5 in manuscript) ----
figS_overall <- ggplot(composition_summary, aes(x = Method, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_publication +
  labs(
    x = "Method",
    y = "Mean Proportion"
  )

ggsave("manuscript/Figure4_composition_barplot_with_cytology.pdf", figS_overall, width = 8, height = 6, device = cairo_pdf)
ggsave("manuscript/Figure4_composition_barplot_with_cytology.png", figS_overall, width = 8, height = 6, dpi = 300)
message("  Saved: Figure4_composition_barplot_with_cytology.pdf/png")

## ===== Print Summary Statistics =====
print("\n=== Summary Statistics ===")

summary_wide <- composition_summary %>%
  pivot_wider(names_from = Method, values_from = Proportion)

print(summary_wide)

# Calculate ranges for BayesPrism Original
bp_ranges <- all_data %>%
  filter(Method == "BayesPrism Original") %>%
  group_by(CellType) %>%
  summarise(
    Min = min(Proportion) * 100,
    Max = max(Proportion) * 100,
    Mean = mean(Proportion) * 100,
    .groups = "drop"
  )

print("\nBayesPrism Original ranges (%):")
print(bp_ranges)

