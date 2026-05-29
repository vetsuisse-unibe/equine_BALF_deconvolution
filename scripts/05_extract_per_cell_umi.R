############################################################
# Per-cell UMI totals + cell type / sample / condition
# from the scRNA-seq reference (input for the mRNA sensitivity analyses)
############################################################

suppressPackageStartupMessages({
  library(Seurat)
})

data_dir <- "data"

# The Seurat object (~5 GB) is not in the repo; see the data accession in the
# manuscript. data/per_cell_umi.csv is provided so scripts 06-07 run without it.
scRNA_ref <- readRDS(file.path(data_dir, "sc22_all_seed.rds"))
DefaultAssay(scRNA_ref) <- "RNA"

meta <- scRNA_ref@meta.data
umi_per_cell <- Matrix::colSums(GetAssayData(scRNA_ref, assay = "RNA", layer = "counts"))

# Disease state is in the object metadata; fall back to merged_metadata.csv
disease_col <- intersect(c("DiseaseState", "Condition", "Group"), colnames(meta))[1]
if (!is.na(disease_col)) {
  condition <- as.character(meta[[disease_col]])
} else {
  mm <- read.csv(file.path(data_dir, "merged_metadata.csv"))
  condition <- mm$DiseaseState[match(rownames(meta), mm$CellID)]
}

per_cell <- data.frame(
  CellID    = rownames(meta),
  CellType  = as.character(meta$MajCellType),
  SampleID  = as.character(meta$sample_id),
  Condition = condition,
  UMI       = umi_per_cell[rownames(meta)]
)

write.csv(per_cell, file.path(data_dir, "per_cell_umi.csv"), row.names = FALSE)
message(nrow(per_cell), " cells written to per_cell_umi.csv")
