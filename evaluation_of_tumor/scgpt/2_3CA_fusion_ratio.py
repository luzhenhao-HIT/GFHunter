import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/scgpt/3CA_test_embedded_with_fusion.h5ad")

# 1. Generate crosstab (Rows: Cell Type, Columns: Fusion Status)
# This step will automatically count how many 'Reliable Fusion' and 'Other Cells' are in each cell type
df_stats = pd.crosstab(adata.obs['knn_prediction'], adata.obs['fusion_status'])

# 2. Ensure 'Reliable Fusion' column exists
# (Prevents column missing errors in extreme cases where no cells have fusions)
if 'Reliable Fusion' not in df_stats.columns:
    df_stats['Reliable Fusion'] = 0

# 3. Calculate ratio
# Numerator: Number of Fusions for this type
# Denominator: Total number of cells for this type (i.e., sum of the row)
df_stats['Total_Cells'] = df_stats.sum(axis=1)
df_stats['Fusion_Ratio'] = df_stats['Reliable Fusion'] / df_stats['Total_Cells']

# 4. Organize results: keep only needed columns and sort by ratio in descending order (optional)
output = df_stats[['Reliable Fusion', 'Total_Cells', 'Fusion_Ratio']].sort_values('Fusion_Ratio', ascending=False)
output.columns = ['Fusion_Count', 'Total_Count', 'Ratio'] # Rename columns to make them more readable

# 5. Print preview
print(output)

# 6. Save as TSV
output.to_csv("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/2_fusion_ratio_by_celltype.tsv", sep='\t')
print("\nResults saved to 2_fusion_ratio_by_celltype.tsv")