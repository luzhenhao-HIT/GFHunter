import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ================= Configuration Paths =================
FUSION_PATH = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/SeekGene_CA.fusion_find.tsv"
EMB_H5AD_PATH = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/scgpt/3CA_test_embedded_with_fusion.h5ad"
RAW_H5AD_PATH = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/CA_SeekGene.umap.h5ad"
OUT_DIR = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output/Step1_DEG"
os.makedirs(OUT_DIR, exist_ok=True)

# ================= 1. Load and Preprocess =================
print("[1/5] Loading data...")
adata_emb = sc.read_h5ad(EMB_H5AD_PATH)
adata_raw = sc.read_h5ad(RAW_H5AD_PATH)

# Transfer annotations
common = adata_emb.obs_names.intersection(adata_raw.obs_names)
adata = adata_raw[common].copy()
adata.obs['knn_prediction'] = adata_emb.obs.loc[common, 'knn_prediction']

# Keep only Malignant cells
adata = adata[adata.obs['knn_prediction'] == 'Malignant'].copy()
print(f"      Number of Malignant cells: {adata.n_obs}")

# ================= 2. Robust Fusion Parsing Logic =================
print("[2/5] Parsing Fusion file...")

# Build ID Mapper (remove suffix for matching)
def strip_suffix(bc): return bc.split('-')[0]
id_mapper = {strip_suffix(bc): bc for bc in adata.obs_names} # Clean -> Real

# Also include Real -> Real mapping in case they already match
for bc in adata.obs_names: id_mapper[bc] = bc

adata.obs['fusion_status'] = 'Non-Fusion'
rf_count, sf_count = 0, 0

if os.path.exists(FUSION_PATH):
    with open(FUSION_PATH, 'r') as f:
        f.readline() # Skip Header
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            if len(parts) < 2: 
                line = f.readline(); continue
            
            fusion = parts[0]
            cells = parts[-2].strip().split(';')
            flag = parts[-1].strip()
            
            for cell in cells:
                cell = cell.strip()
                if not cell: continue
                
                # ID matching
                matched_id = id_mapper.get(cell) or id_mapper.get(strip_suffix(cell))
                
                if matched_id:
                    if flag == 'Reliable Fusion':
                        adata.obs.loc[matched_id, 'fusion_status'] = 'Fusion-Positive' # Unified labeling
                        rf_count += 1
                    elif flag == 'Suspected Fusion':
                        # Do not overwrite Reliable
                        if adata.obs.loc[matched_id, 'fusion_status'] != 'Fusion-Positive':
                            adata.obs.loc[matched_id, 'fusion_status'] = 'Fusion-Positive'
                            sf_count += 1
            line = f.readline()
else:
    print("Error: Fusion file not found.")
    exit()

print(f"      Matching stats: Total Fusion+ cells = {sum(adata.obs['fusion_status']=='Fusion-Positive')}")

# ================= 3. Differential Expression Analysis =================
print("[3/5] Performing differential expression analysis (Malignant: Fusion+ vs Non-Fusion)...")

# Normalization (if needed)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
adata.raw = adata

sc.tl.rank_genes_groups(
    adata, 
    groupby='fusion_status', 
    groups=['Fusion-Positive'], 
    reference='Non-Fusion', 
    method='wilcoxon',
    pts=True
)

# ================= 4. Save Results =================
deg_df = sc.get.rank_genes_groups_df(adata, group="Fusion-Positive")
csv_path = f"{OUT_DIR}/1_deg_results.csv"
deg_df.to_csv(csv_path, index=False)
print(f"      DEG results saved to: {csv_path}")

# Save annotated h5ad for downstream steps
h5ad_path = f"{OUT_DIR}/1_malignant_fusion.h5ad"
adata.write(h5ad_path)
print(f"      H5AD saved to: {h5ad_path}")

# ================= 5. Volcano Plot =================
plt.figure(figsize=(8, 6))
df = deg_df.dropna()
colors = ['red' if (p < 0.05 and l > 0.5) else 'blue' if (p < 0.05 and l < -0.5) else 'grey' 
          for p, l in zip(df['pvals_adj'], df['logfoldchanges'])]
plt.scatter(df['logfoldchanges'], -np.log10(df['pvals_adj'] + 1e-300), c=colors, s=5)
plt.title("Volcano: Malignant Fusion+ vs Non-Fusion")
plt.savefig(f"{OUT_DIR}/1_volcano.png")