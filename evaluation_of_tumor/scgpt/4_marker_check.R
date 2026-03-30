# 02_dotplot_with_knn_csv.R
# Function: Read CSV files, build a Seurat object, and use the 'knn_prediction' column to draw a DotPlot.
# Extra feature: Export the DotPlot drawing data to a TSV file.

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix) 

# --- User-configurable parameters ---
# Expression matrix file path
EXPRESSION_CSV_PATH <- "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/3_expression_matrix.csv"
# Metadata file path
METADATA_CSV_PATH <- "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/3_cell_metadata_merged.csv"
# Annotation column name (consistent with the Python script)
cell_type_col_name <- "knn_prediction" 
# Output filenames
PLOT_FILENAME <- "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/4_DotPlot_KnnPrediction_CSV.pdf"
DATA_FILENAME <- "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/4_DotPlot_KnnPrediction_CSV_Data.tsv"
RDS_FILENAME  <- "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/4_scRNA_processed.rds"
# --------------------------

# 1. Read data
print("Reading expression matrix and metadata CSV files...")

# 1.1 Read expression matrix
expression_matrix <- read.csv(EXPRESSION_CSV_PATH, row.names = 1, check.names = FALSE) 
counts <- as(as.matrix(expression_matrix), "dgCMatrix")

# 1.2 Read cell metadata
cell_metadata <- read.csv(METADATA_CSV_PATH, row.names = 1, check.names = FALSE)

# 2. Build Seurat object
print("Building Seurat object...")
scRNA <- CreateSeuratObject(counts = counts, project = "KnnDotPlot", assay = "RNA")

# 3. Add metadata to Seurat object
print("Merging metadata...")
if (all(rownames(cell_metadata) %in% colnames(scRNA))) {
    cols_to_merge <- setdiff(colnames(cell_metadata), colnames(scRNA@meta.data))
    if (!cell_type_col_name %in% cols_to_merge) {
        cols_to_merge <- c(cols_to_merge, cell_type_col_name)
    }
    # Select columns to merge, drop = FALSE ensures format
    metadata_to_add <- cell_metadata[colnames(scRNA), cols_to_merge, drop = FALSE]
    scRNA@meta.data <- cbind(scRNA@meta.data, metadata_to_add)
} else {
    stop("Error: Cell barcodes in metadata do not match those in the expression matrix!")
}

# 4. Preprocessing
print("Normalizing and scaling data...")
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA) 

# 5. Check annotation column
if (!cell_type_col_name %in% colnames(scRNA@meta.data)) {
    stop(paste0("Error: Annotation column '", cell_type_col_name, "' not found in the Seurat object."))
}
print(paste0("Found annotation column: '", cell_type_col_name, "', preparing to plot."))

# 6. Define Marker list 
markers <- list(
     "T_Cells" = c("CD3E","CD3D"),
     "B_Cells" = c("CD79A","CD79B","MS4A1"),
     "Epithelial_Cells" = c("EPCAM","KRT7"),
     "NK_Cells" = c("NKG7","GNLY","GZMA"),
     "Netrophils" = c("S100A8","S100A9","CSF3R"),
     "Dendritic_Cells" = c("CD1C","CD1A","ITGAX"),
     "Macrophages" = c("C1QC","C1QB","C1QA"),
     "Fibroblasts" = c("COL1A2", "COL1A1", "DCN","LUM","COL6A1"),
     "Endothelial_Cells" = c("PECAM1","CLDN5","VWF"),
     "Monocytes" = c("CD14","CCR2","CXCR4","FCN1"),
     "Mast_Cells" = c("KIT","TPSAB1","TPSB2"),
     "Plasma_Cells" = c("MZB1","IGKC","JCHAIN","IGHG1","IGHA1")
)

# 7. Plot and save
print("Generating DotPlot...")

# Assign the DotPlot object to variable p
p <- DotPlot(scRNA, features = markers, group.by = cell_type_col_name) + 
      RotatedAxis() + 
      ggtitle(paste("Marker Genes by", cell_type_col_name, "Annotation")) +
      # Note: If markers is a named list, DotPlot automatically adds feature.groups
      facet_grid(~feature.groups, scales = "free_x", space = "free_x") +
      theme(
        axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size = 10, angle = 45, hjust = 0, vjust = 0) 
      )

# Save PDF
pdf(PLOT_FILENAME, width = 25, height = 10) 
print(p)
dev.off()
print(paste0("Image saved to: ", PLOT_FILENAME))

# ---------------------------------------------------------
# [New Feature] Extract and save data
# ---------------------------------------------------------
print("Exporting underlying DotPlot data...")

# Extract data from the ggplot object
dotplot_data <- p$data

# Simply rename columns to improve readability (Seurat default column names are clear, usually don't need renaming, but good to confirm)
# avg.exp: Average expression (color intensity)
# pct.exp: Percent expressed (dot size)
# features.plot: Gene name
# id: Grouping (i.e., knn_prediction)
# feature.groups: Gene grouping (e.g., T_Cells, B_Cells, etc.)

# Sort by id (cell type) and features.plot (gene) to make the table neat
dotplot_data <- dotplot_data %>% 
  arrange(id, features.plot) %>%
  select(features.plot, id, avg.exp, pct.exp, avg.exp.scaled, feature.groups)

# Save as TSV
write.table(dotplot_data, file = DATA_FILENAME, sep = "\t", quote = FALSE, row.names = FALSE)

print(paste0("Data saved to: ", DATA_FILENAME))
print("Script execution completed!")

# 9. Save Seurat object (RDS) - [New step]
print("9. Saving Seurat object as RDS file...")
# This step will include counts, data (normalized), scale.data, and complete meta.data
saveRDS(scRNA, file = RDS_FILENAME)

print(paste0("   - RDS object saved to: ", RDS_FILENAME))
print("=== All tasks completed ===")