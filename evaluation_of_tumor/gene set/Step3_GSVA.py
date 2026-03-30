import scanpy as sc
import gseapy as gp
import pandas as pd
import numpy as np
import os

# ================= Configuration =================
# Read the H5AD generated in the first step directly to save the trouble of repeated filtering and mapping
INPUT_H5AD = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output/Step1_DEG/1_malignant_fusion.h5ad"
GMT_PATH = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/code/h.all.v2025.1.Hs.symbols.gmt"
OUT_DIR = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output/Step3_GSVA"
os.makedirs(OUT_DIR, exist_ok=True)

# ================= Prepare Data =================
print("Reading Malignant H5AD...")
adata = sc.read_h5ad(INPUT_H5AD)

print("Preparing GSVA input matrix...")
# Deduplicate gene names
if not adata.var_names.is_unique:
    adata.var_names_make_unique()

# Optional for acceleration: use only highly variable genes (e.g., Top 5000)
# sc.pp.highly_variable_genes(adata, n_top_genes=5000)
# genes_use = adata.var_names[adata.var['highly_variable']]
genes_use = adata.var_names # Using all genes for demonstration

gsva_df = pd.DataFrame(
    adata[:, genes_use].X.toarray().T if hasattr(adata.X, "toarray") else adata[:, genes_use].X.T,
    index=genes_use.astype(str), # Must be converted to strings
    columns=adata.obs_names
)

# ================= Run GSVA =================
print("Running GSVA...")
try:
    gsva_res = gp.gsva(
        data=gsva_df,
        gene_sets=GMT_PATH,
        outdir=OUT_DIR,
        method='gsva',
        kcdf='Gaussian', # Use Gaussian for log-transformed data
        threads=32
    )
    
    # Save results
    res_path = f"{OUT_DIR}/3_gsva_matrix.csv"
    # res2d is a long table (Term, Name, ES) or a matrix, depending on the version
    # Here we force saving as a matrix format: Index=Term, Columns=Cells
    if hasattr(gsva_res, 'res2d') and 'Name' in gsva_res.res2d.columns: # If it is a long table
        matrix = gsva_res.res2d.pivot(index='Term', columns='Name', values='ES')
    else:
        matrix = gsva_res.res2d
        
    matrix.to_csv(res_path)
    print(f"GSVA result matrix saved to: {res_path}")

except Exception as e:
    print(f"GSVA error: {e}")