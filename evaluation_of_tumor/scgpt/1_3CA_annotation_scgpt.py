import sys
import os
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mode
import faiss

import scgpt as scg
import torch 

# Ignore warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=UserWarning, module='scanpy')

# --- 1. Global Configuration ---
CONFIG = {
    # Path settings
    "model_dir": Path("/home/user/luzhenhao/cite_handle/scgpt-model/whole-human"),
    "ref_adata_path": Path("/home/user/luzhenhao/time/data/Cancer_altas/Data_Lung/3CA_Lung_Atlas_Integrated_cleaned.h5ad"),
    "test_adata_path": Path("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/CA_SeekGene.umap.h5ad"),
    "fusion_file": Path("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/SeekGene_CA.fusion_find.tsv"),
    
    # Output paths
    "ref_embed_save_path": Path("/home/user/luzhenhao/cite_handle/Cancer_map/data_Map/reference_embedded_cleaned.h5ad"),
    "test_embed_save_path": Path("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/scgpt/3CA_test_embedded_with_fusion.h5ad"),
    "output_pdf_path": Path("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/1_3CA_SeekGene_CA.scGPT_and_Fusion_Analysis.pdf"),

    # Parameter settings
    "cell_type_key": "cell_type",
    "gene_col": "index",
    "batch_size": 64, # Keep small to prevent OOM (Out of Memory)
    "n_neighbors": 15,
    "knn_k": 10,
}

# --- 2. Helper Functions ---
def load_or_embed(adata_path, save_path, model_dir, gene_col, obs_to_save=None):
    recompute = False
    
    if save_path.exists():
        # Check timestamp: if the input file is newer than the cache, recompute
        try:
            input_mtime = adata_path.stat().st_mtime
            cache_mtime = save_path.stat().st_mtime
            if input_mtime > cache_mtime:
                print("⚠️ Input file update detected, preparing to recompute embeddings...")
                recompute = True
            else:
                print(f"Valid embedding cache detected: {save_path}, loading directly...")
                adata = sc.read_h5ad(save_path)
        except Exception:
            recompute = True
    else:
        recompute = True

    if recompute:
        print(f"Loading raw data: {adata_path}")
        adata = sc.read_h5ad(adata_path)
        
        # Clear GPU memory
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        print("Generating embeddings (ScGPT)...")
        adata = scg.tasks.embed_data(
            adata,
            model_dir,
            gene_col=gene_col,
            obs_to_save=obs_to_save,
            batch_size=CONFIG["batch_size"],
            return_new_adata=True,
        )
        print(f"Saving embedding results to: {save_path}")
        save_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(save_path)
    return adata

def sort_adata_for_plotting(adata_obj, col_name, order_map):
    adata_sorted = adata_obj.copy()
    temp_series = adata_sorted.obs[col_name].astype(str)
    adata_sorted.obs['plot_order'] = temp_series.map(order_map).fillna(0)
    return adata_sorted[adata_sorted.obs['plot_order'].argsort()]

def strip_suffix(barcode):
    """Remove barcode suffix"""
    return barcode.split('-')[0].strip()

# --- 3. Data Preparation & Fixes ---
print("=== Step 1/5: Prepare Data ===")
ref_embed_adata = load_or_embed(
    CONFIG["ref_adata_path"], CONFIG["ref_embed_save_path"],
    CONFIG["model_dir"], gene_col=CONFIG["gene_col"], obs_to_save=CONFIG["cell_type_key"]
)

test_embed_adata = load_or_embed(
    CONFIG["test_adata_path"], CONFIG["test_embed_save_path"], 
    CONFIG["model_dir"], gene_col=CONFIG["gene_col"], obs_to_save=None
)

# === [Fix]: Restore original cell barcodes ===
print("\nRestoring cell indices from original file (Fixing Barcodes)...")
original_test_adata = sc.read_h5ad(CONFIG["test_adata_path"], backed='r')
if test_embed_adata.shape[0] != original_test_adata.shape[0]:
    raise ValueError("Critical Error: Inconsistent cell counts!")
test_embed_adata.obs_names = original_test_adata.obs_names
print("✅ Cell barcodes fixed.")


# --- 4. KNN Prediction ---
print("\n=== Step 2/5: KNN Prediction ===")
if 'knn_prediction' not in test_embed_adata.obs.columns:
    X_ref, X_test = ref_embed_adata.X.astype('float32'), test_embed_adata.X.astype('float32')
    
    # Faiss Indexing
    index = faiss.IndexFlatL2(X_ref.shape[1])
    index.add(X_ref)
    distances, indices = index.search(X_test, CONFIG['knn_k'])

    # [Core Fix] Force conversion to string to solve '<' not supported between 'str' and 'float'
    # .astype(str) converts NaN to 'nan', allowing for sorting
    ref_labels = ref_embed_adata.obs[CONFIG["cell_type_key"]].astype(str).to_numpy()
    
    unique_classes, encoded_ref_labels = np.unique(ref_labels, return_inverse=True)
    neighbor_encoded = encoded_ref_labels[indices]
    mode_result = mode(neighbor_encoded, axis=1, keepdims=True)
    preds = unique_classes[mode_result.mode.flatten()]
    
    test_embed_adata.obs['knn_prediction'] = preds
    print("KNN prediction complete.")
else:
    print("Skipping KNN calculation (already exists).")

if 'X_umap' not in test_embed_adata.obsm:
    print("Calculating UMAP...")
    sc.pp.neighbors(test_embed_adata, use_rep='X', n_neighbors=CONFIG["n_neighbors"])
    sc.tl.umap(test_embed_adata)

# --- 5. Parse Fusion Data ---
print("\n=== Step 3/5: Parse Fusion File and Label Cells ===")
fusion_details = {} 
test_embed_adata.obs['fusion_status'] = 'Other Cells'

# Build standard ID mapping dictionary (original and suffix-removed only)
id_mapper = {}
for real_id in test_embed_adata.obs_names:
    id_mapper[real_id] = real_id # Original full
    clean_key = strip_suffix(real_id)
    id_mapper[clean_key] = real_id # Suffix removed
    
rf_count = 0
sf_count = 0

if CONFIG["fusion_file"].exists():
    with open(CONFIG["fusion_file"], 'r') as f:
        f.readline()
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            if len(parts) < 2: 
                line = f.readline(); continue
                
            fusion = parts[0]
            cells = parts[-2].strip().split(';')
            flag = parts[-1].strip()

            if fusion not in fusion_details:
                fusion_details[fusion] = {'RF': [], 'SF': []}

            for cell in cells:
                cell = cell.strip()
                if not cell: continue

                # Matching logic
                matched_real_id = None
                if cell in id_mapper:
                    matched_real_id = id_mapper[cell]
                else:
                    clean_cell = strip_suffix(cell)
                    if clean_cell in id_mapper:
                        matched_real_id = id_mapper[clean_cell]

                if matched_real_id:
                    if flag == 'Reliable Fusion':
                        fusion_details[fusion]['RF'].append(matched_real_id)
                        test_embed_adata.obs.loc[matched_real_id, 'fusion_status'] = 'Reliable Fusion'
                        rf_count += 1
                    elif flag == 'Suspected Fusion':
                        fusion_details[fusion]['SF'].append(matched_real_id)
                        if test_embed_adata.obs.loc[matched_real_id, 'fusion_status'] != 'Reliable Fusion':
                            test_embed_adata.obs.loc[matched_real_id, 'fusion_status'] = 'Suspected Fusion'
                            sf_count += 1
            
            line = f.readline()
            
    print(f"Matching stats: Reliable matches = {rf_count}, Suspected matches = {sf_count}")
    if rf_count == 0 and sf_count == 0:
        print("❌ Critical Warning: Match count is still 0! Please check your data.")
else:
    print("Fusion file not found.")

# --- 6. Plot and Save PDF ---
print(f"\n=== Step 4/5: Generate PDF ({CONFIG['output_pdf_path']}) ===")
CONFIG["output_pdf_path"].parent.mkdir(parents=True, exist_ok=True)
test_embed_adata.obs['knn_prediction'] = test_embed_adata.obs['knn_prediction'].astype('category')

with PdfPages(CONFIG["output_pdf_path"]) as pdf:
    
    # Page 1
    print("Plotting Page 1...")
    fig, ax = plt.subplots(figsize=(10, 8)) 
    sc.pl.umap(
        test_embed_adata,
        color='knn_prediction',
        ax=ax, show=False, title="Page 1: scGPT KNN Predicted Cell Types", frameon=False,
        legend_loc="right margin"
    )
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

    # Page 2
    print("Plotting Page 2...")
    order_map = {'Other Cells': 0, 'Suspected Fusion': 1, 'Reliable Fusion': 2}
    adata_sorted = sort_adata_for_plotting(test_embed_adata, 'fusion_status', order_map)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata_sorted,
        color="fusion_status",
        palette={'Other Cells': 'lightgrey', 'Suspected Fusion': 'orange', 'Reliable Fusion': 'red'},
        title='Page 2: Global Fusion Status',
        show=False, ax=ax, frameon=False,
        legend_loc="right margin" 
    )
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

    # Page 3+
    sorted_fusions = sorted(fusion_details.items(), key=lambda x: len(x[1]['RF']) + len(x[1]['SF']), reverse=True)
    print(f"Plotting individual fusions ({len(sorted_fusions)} total)...")
    
    for i, (fusion_name, cell_dict) in enumerate(sorted_fusions):
        if len(cell_dict['RF']) + len(cell_dict['SF']) == 0: continue
        
        test_embed_adata.obs['current'] = 'Other Cells'
        for c in cell_dict['SF']: test_embed_adata.obs.loc[c, 'current'] = 'Suspected'
        for c in cell_dict['RF']: test_embed_adata.obs.loc[c, 'current'] = 'Reliable'
            
        adata_plot = sort_adata_for_plotting(test_embed_adata, 'current', {'Other Cells':0, 'Suspected':1, 'Reliable':2})
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(
            adata_plot,
            color='current',
            palette={'Other Cells': 'lightgrey', 'Suspected': 'orange', 'Reliable': 'red'},
            title=f'{fusion_name} (R:{len(cell_dict["RF"])}, S:{len(cell_dict["SF"])})',
            show=False, ax=ax, frameon=False,
            legend_loc="right margin"
        )
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()
        del test_embed_adata.obs['current']

# --- 7. Save ---
print("\n=== Step 5/5: Save Results ===")
if 'plot_order' in test_embed_adata.obs: del test_embed_adata.obs['plot_order']
test_embed_adata.write_h5ad(CONFIG['test_embed_save_path'])
print("🎉 Done!")