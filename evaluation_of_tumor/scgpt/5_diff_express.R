library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel) # Used to prevent overlap of volcano plot labels. If not installed, run install.packages("ggrepel")

# ==========================================
# 1. Preparation
# ==========================================
# Assume the scRNA object already exists (from your previous script)
# If this is a new script, load the saved object via readRDS first, or rerun the build steps

# Set parameters
target_col <- "knn_prediction"  # Grouping column
group_1 <- "Malignant"          # Target group (positive values mean high expression here)
group_2 <- "Epithelial"         # Control group (negative values mean high expression here)

# Ensure the Seurat object uses the correct column as identity
print("Loading Seurat object...")
scRNA <- readRDS("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/4_scRNA_processed.rds")

Idents(scRNA) <- target_col

print(paste0("Current grouping column: ", target_col))
print(paste0("Comparing: ", group_1, " vs ", group_2))

# Check if these two groups exist in the data
available_groups <- levels(Idents(scRNA))
if (!(group_1 %in% available_groups && group_2 %in% available_groups)) {
  stop(paste("Error: Cannot find", group_1, "or", group_2, "in the data. Existing groups include:", paste(available_groups, collapse = ", ")))
}

# ==========================================
# 2. Perform Differential Expression Analysis (FindMarkers)
# ==========================================
print("Calculating differentially expressed genes (Wilcoxon test)...")

# logfc.threshold = 0.25 is Seurat's default recommendation, filtering genes with at least 0.25 log2 difference
# min.pct = 0.1 means the gene is expressed in at least 10% of cells
de_markers <- FindMarkers(scRNA, 
                          ident.1 = group_1, 
                          ident.2 = group_2, 
                          test.use = "wilcox", 
                          logfc.threshold = 0.25,
                          min.pct = 0.1)

# Organize results
de_markers$gene <- rownames(de_markers)
# Add up/down regulation status
de_markers$diff_type <- "NS"
de_markers$diff_type[de_markers$avg_log2FC > 0.5 & de_markers$p_val_adj < 0.05] <- paste0("Up in ", group_1)
de_markers$diff_type[de_markers$avg_log2FC < -0.5 & de_markers$p_val_adj < 0.05] <- paste0("Up in ", group_2)

# ==========================================
# 3. Save Results
# ==========================================
output_filename <- paste0("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/5_DE_", group_1, "_vs_", group_2, ".csv")
write.csv(de_markers, output_filename, row.names = FALSE)
print(paste("Differential analysis results saved as:", output_filename))

# ==========================================
# 4. Visualization: Volcano Plot
# ==========================================
print("Drawing Volcano plot...")

# Extract the most significant genes for labeling (e.g., top 10 highest and lowest FC)
top_genes <- de_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(diff_type) %>%
  top_n(n = 10, wt = abs(avg_log2FC)) %>%
  pull(gene)

# Plotting
p_volcano <- ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  # Draw points
  geom_point(aes(color = diff_type), alpha = 0.6, size = 1.5) +
  # Set colors
  scale_color_manual(values = c("NS" = "grey", 
                                "Up in Malignant" = "#E41A1C", # Red
                                "Up in Epithelial" = "#377EB8")) + # Blue
  # Add threshold lines
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha=0.5) +
  # Label gene names (only significant ones)
  geom_text_repel(data = subset(de_markers, gene %in% top_genes & diff_type != "NS"),
                  aes(label = gene), 
                  max.overlaps = 15,
                  box.padding = 0.5) +
  # Titles and labels
  labs(title = paste0("Volcano Plot: ", group_1, " vs ", group_2),
       x = "Average Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Expression Status") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Save image
ggsave(paste0("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/5_Volcano_", group_1, "_vs_", group_2, ".pdf"), p_volcano, width = 8, height = 6)
print(paste0("Volcano plot saved as: Volcano_", group_1, "_vs_", group_2, ".pdf"))

# ==========================================
# 5. Visualization: Violin Plot (Top Markers)
# ==========================================
# Select the top 6 genes with the largest difference to display
top_6_genes <- head(de_markers[order(de_markers$p_val_adj, -abs(de_markers$avg_log2FC)), "gene"], 6)

if (length(top_6_genes) > 0) {
    pdf(paste0("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/5_VlnPlot_", group_1, "_vs_", group_2, ".pdf"), width = 10, height = 8)
    print(VlnPlot(scRNA, features = top_6_genes, group.by = target_col, pt.size = 0.1, idents = c(group_1, group_2)))
    dev.off()
    print("Top genes Violin plot saved.")
}

# ==========================================
# 6. Extra Task: Exclusive Volcano Plot for Genes of Interest
# ==========================================
print("Drawing exclusive Volcano plot for genes of interest...")

# Define your list of genes of interest
interest_genes <- c(
  "PYROXD1", "PIK3C2G", "IQGAP1", "CRTC3", "BACH1", 
  "USP16",   "ADD2",    "CLEC4F", "FIGLA", "BCAS3", 
  "ST6GALNAC2", "CKAP2", "ELAVL2", "HSD17B7", "PLD5P1", 
  "PCDH9",   "B3GLCT"
)

# Find which of these genes exist in the differential analysis results table
# Note: If gene expression is too low or difference is too small, they might be filtered by FindMarkers, so we need to intersect
genes_to_label <- intersect(interest_genes, de_markers$gene)

if (length(genes_to_label) > 0) {
  
  # Create a new column for labeling, only genes of interest get a Label, others get NA
  de_markers$label_text <- NA
  de_markers$label_text[de_markers$gene %in% genes_to_label] <- de_markers$gene[de_markers$gene %in% genes_to_label]
  
  # To make genes of interest more obvious, we can bold these points
  de_markers$is_interest <- ifelse(de_markers$gene %in% genes_to_label, "Yes", "No")

  p_volcano_interest <- ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    # 1. Draw background points (non-interest genes, grey or light colored)
    geom_point(data = subset(de_markers, is_interest == "No"), 
               aes(color = diff_type), alpha = 0.3, size = 1) +
    
    # 2. Draw points for genes of interest (opaque, slightly larger, black border)
    geom_point(data = subset(de_markers, is_interest == "Yes"), 
               aes(fill = diff_type), shape = 21, color = "black", size = 3, stroke = 1) +
    
    # Color settings (keep consistent with previous)
    scale_color_manual(values = c("NS" = "grey", "Up in Malignant" = "#E41A1C", "Up in Epithelial" = "#377EB8")) +
    scale_fill_manual(values = c("NS" = "grey", "Up in Malignant" = "#E41A1C", "Up in Epithelial" = "#377EB8")) +
    
    # Guidelines
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha=0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha=0.5) +
    
    # 3. Labels: only label genes of interest
    geom_text_repel(aes(label = label_text), 
                    max.overlaps = 50, # Allow more label overlap
                    box.padding = 0.5,
                    min.segment.length = 0, # Always show line segments
                    fontface = "bold",
                    color = "black") +
    
    labs(title = paste0("Volcano Plot (Interest Genes): ", group_1, " vs ", group_2),
         subtitle = "Only selected genes are labeled",
         x = "Average Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  # Save
  ggsave(paste0("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/5_Volcano_InterestGenes_", group_1, "_vs_", group_2, ".pdf"), p_volcano_interest, width = 8, height = 6)
  print(paste0("Volcano plot for genes of interest saved. Labeled ", length(genes_to_label), " genes."))
  
} else {
  print("Note: None of your genes of interest passed the filtering threshold for differential analysis (possibly due to low expression or insignificant difference), cannot be marked on the Volcano plot.")
}

# ==========================================
# 7. Extra Task: Violin Plot for Genes of Interest
# ==========================================
print("Drawing Violin plot for genes of interest...")

# For Violin plots, we need to find genes in the original Seurat object (can be plotted even if not significantly different)
valid_vln_genes <- intersect(interest_genes, rownames(scRNA))

if (length(valid_vln_genes) > 0) {
  
  # To prevent the image from being too large, we page the genes or adjust the image height
  # Here we set 3 per row, automatically calculate height
  n_genes <- length(valid_vln_genes)
  plot_height <- ceiling(n_genes / 3) * 4 # Dynamic height
  
  pdf(paste0("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/5_VlnPlot_InterestGenes_", group_1, "_vs_", group_2, ".pdf"), width = 12, height = plot_height)
  
  # stack = TRUE (Seurat V5) or combine (older versions)
  # Using purrr::map or directly passing the vector usually works, here we use the base method
  print(
    VlnPlot(scRNA, 
            features = valid_vln_genes, 
            group.by = target_col, 
            idents = c(group_1, group_2), # Only display these two groups
            pt.size = 0.1, 
            ncol = 3) # 3 per row
  )
  
  dev.off()
  print(paste0("Violin plot for genes of interest saved, drawing a total of ", n_genes, " genes."))
  
} else {
  print("Error: None of your genes of interest were found in the Seurat object.")
}

print("=== All analysis tasks completed ===")