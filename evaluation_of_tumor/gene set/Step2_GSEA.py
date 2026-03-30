import gseapy as gp
import pandas as pd
import os
import matplotlib.pyplot as plt

# ================= Configuration =================
DEG_CSV = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output/Step1_DEG/1_deg_results.csv"
GMT_PATH = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/code/h.all.v2025.1.Hs.symbols.gmt"
OUT_DIR = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output/Step2_GSEA"
os.makedirs(OUT_DIR, exist_ok=True)

# ================= Run GSEA =================
print("Reading DEG results...")
df = pd.read_csv(DEG_CSV)

# 1. Prepare Ranked List
# Using scores (Z-score) is highly recommended as it combines fold change and P-value
rank_data = df[['names', 'scores']].sort_values('scores', ascending=False)

print("Running GSEA Preranked...")
pre_res = gp.prerank(
    rnk=rank_data,
    gene_sets=GMT_PATH,
    threads=8,
    min_size=5,
    max_size=500,
    permutation_num=1000,
    outdir=OUT_DIR,
    seed=42,
    verbose=True
)

# Save significant results
res_df = pre_res.res2d
# Compatible column name
fdr_col = next(c for c in res_df.columns if 'fdr' in c.lower() or 'q-val' in c.lower())
sig_df = res_df[res_df[fdr_col] < 0.25].sort_values(fdr_col) # FDR < 0.25 (GSEA officially recommended relaxed threshold)

sig_path = f"{OUT_DIR}/2_gsea_significant.csv"
sig_df.to_csv(sig_path)
print(f"GSEA analysis complete. Significant pathways saved to: {sig_path}")