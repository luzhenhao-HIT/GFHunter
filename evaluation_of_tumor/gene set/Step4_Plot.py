import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scanpy as sc
import os
import numpy as np
from scipy.stats import mannwhitneyu

# ================= Configuration Paths =================
BASE_DIR = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output"
GSEA_SIG_CSV = os.path.join(BASE_DIR, "Step2_GSEA/2_gsea_significant.csv")
GSVA_MATRIX_CSV = os.path.join(BASE_DIR, "Step3_GSVA/3_gsva_matrix.csv")
ADATA_PATH = os.path.join(BASE_DIR, "Step1_DEG/1_malignant_fusion.h5ad")

# Output directories
OUT_DIR = os.path.join(BASE_DIR, "Step4_Plot_Final_Fixed")
SINGLE_FIG_DIR = os.path.join(OUT_DIR, "Single_Violins")
DATA_DIR = os.path.join(OUT_DIR, "Plot_Data")

os.makedirs(SINGLE_FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# ================= 1. Read and Integrate Data =================
print("Reading data...")
adata = sc.read_h5ad(ADATA_PATH)
meta = adata.obs[['fusion_status']].copy()

# Read GSVA matrix
gsva_mat = pd.read_csv(GSVA_MATRIX_CSV, index_col=0)

# Read and sort GSEA results
gsea_df = pd.read_csv(GSEA_SIG_CSV)
nes_col = next(c for c in gsea_df.columns if 'nes' in c.lower())
term_col = next(c for c in gsea_df.columns if 'term' in c.lower())

# Get all significant pathways (descending by NES)
target_pathways = gsea_df.sort_values(nes_col, ascending=False)[term_col].tolist()
print(f"Planning to plot {len(target_pathways)} pathways in total.")

# Merge data
plot_data = gsva_mat.T.join(meta, how='inner')

if plot_data.empty:
    print("Error: Cell IDs do not match, unable to plot.")
    exit()

# ================= 2. Loop Plotting + Statistics + Export Data =================
print(f"Starting to plot, images will be saved to: {SINGLE_FIG_DIR}")
print(f"Source data will be saved to: {DATA_DIR}")

my_palette = {"Fusion-Positive": "#E41A1C", "Non-Fusion": "#377EB8"}

def get_significance_stars(p_value):
    if p_value < 0.0001:
        return "****"
    elif p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"  # Or return "" to indicate no stars

ps = []
for i, pathway in enumerate(target_pathways):
    if pathway not in plot_data.columns:
        continue
        
    # Prepare data for the current pathway
    subset = plot_data[[pathway, 'fusion_status']].copy()
    subset.columns = ['Score', 'Status']
    
    # --- A. Export plot data ---
    clean_name = pathway.replace("HALLMARK_", "")
    safe_name = f"{i+1:02d}_{clean_name}"
    
    csv_filename = os.path.join(DATA_DIR, f"{safe_name}.csv")
    subset.to_csv(csv_filename)
    
    # --- B. Perform statistical tests ---
    group_pos = subset[subset['Status'] == 'Fusion-Positive']['Score']
    group_neg = subset[subset['Status'] == 'Non-Fusion']['Score']
    
    p_text = "P = N/A"
    if len(group_pos) > 0 and len(group_neg) > 0:
        # Use Mann-Whitney U test (Wilcoxon rank-sum test)
        stat, pval = mannwhitneyu(group_pos.values, group_neg.values,
                                  alternative='two-sided')
        # Get significance stars
        stars = get_significance_stars(pval)
        
        # Format output
        if pval < 0.0001:
            p_text = f"P < 0.0001{stars}"
        elif pval < 0.001:
            p_text = f"P < 0.001{stars}"
        elif pval < 0.01:
            p_text = f"P < 0.01{stars}"
        elif pval < 0.05:
            p_text = f"P < 0.05{stars}"
        else:
            p_text = f"P = {pval:.3f}{stars}"
        
        ps.append((safe_name, pval, p_text))

    # --- C. Plotting ---
    plt.figure(figsize=(5, 6))
    
    sns.violinplot(
        data=subset, 
        x='Status', 
        y='Score', 
        hue='Status',        # [Critical Fix] Explicitly specify hue
        legend=False,        # [Critical Fix] Hide redundant legend
        palette=my_palette,
        inner="box",         
        linewidth=1.5,
        density_norm='width'              
    )
    
    # --- D. Beautification ---
    plt.title(f"{clean_name}\n({p_text})", fontsize=12, fontweight='bold')
    plt.ylabel("GSVA Score", fontsize=11)
    plt.xlabel("", fontsize=11)
    
    sns.despine()
    
    # Force x-axis labels to be horizontal (no rotation)
    plt.xticks(rotation=0, fontsize=10)
    
    # Add margins to prevent cropping
    plt.margins(y=0.1)
    
    # Save image
    save_path = os.path.join(SINGLE_FIG_DIR, f"{safe_name}.png")
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.savefig(save_path.replace('.png', '.pdf'), bbox_inches='tight')
    
    plt.close()

with open(os.path.join(OUT_DIR, "Significance_Results.tsv"), 'w') as f:
    f.write("Pathway\tP-value\tAnnotation\n")
    for name, pval, ptext in ps:
        f.write(f"{name}\t{pval:.6e}\t{ptext}\n")
print("=== All plotting and data exporting completed ===")