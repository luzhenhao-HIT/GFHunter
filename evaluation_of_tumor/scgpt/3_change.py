# 01_merge_annotation_and_export_csv.py
# Function: Read the original H5AD and the predicted H5AD, align cells, merge the 'knn_prediction' column, and export the expression matrix and metadata to CSV.

import anndata
import pandas as pd
import numpy as np
import sys
import os

# --- User-configurable parameters (Be sure to modify ORIGINAL_H5AD_PATH) ---
# Original H5AD file path [!!! PLEASE MODIFY THIS PATH !!!]
ORIGINAL_H5AD_PATH = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/CA_SeekGene.umap.h5ad" 
# H5AD file path containing new annotation information
PREDICTION_H5AD_PATH = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/scgpt/3CA_test_embedded_with_fusion.h5ad"

# New annotation column name (extracted from the prediction file)
NEW_CELL_TYPE_COL = "knn_prediction"
# ----------------------------------------------------

try:
    # 1. Read data
    print(f"Reading original H5AD file: {ORIGINAL_H5AD_PATH}")
    original_adata = anndata.read_h5ad(ORIGINAL_H5AD_PATH)
    print(f"Reading predicted H5AD file: {PREDICTION_H5AD_PATH}")
    prediction_adata = anndata.read_h5ad(PREDICTION_H5AD_PATH)
    
    # 2. Check and extract prediction column
    if NEW_CELL_TYPE_COL not in prediction_adata.obs.columns:
        print(f"Error: Required column '{NEW_CELL_TYPE_COL}' not found in prediction file '{PREDICTION_H5AD_PATH}'.")
        sys.exit(1)
        
    print(f"Prediction column '{NEW_CELL_TYPE_COL}' found, preparing to merge.")
    new_predictions = prediction_adata.obs[[NEW_CELL_TYPE_COL]]
    
    # 3. Align cells and merge annotations
    # Use reindex to align barcodes
    aligned_predictions = new_predictions.reindex(original_adata.obs_names)

    # Add new annotations to the obs of the original adata
    original_adata.obs[NEW_CELL_TYPE_COL] = aligned_predictions[NEW_CELL_TYPE_COL]
    
    print(f"New annotation column merged into original data with column name: '{NEW_CELL_TYPE_COL}'")

    # 4. Export expression matrix CSV
    print("Exporting 1. expression_matrix.csv...")
    
    # Convert to dense matrix and transpose (R/Seurat convention: rows are genes, columns are cells)
    if isinstance(original_adata.X, (np.ndarray, pd.DataFrame)):
        expression_data = original_adata.X
    else:
        expression_data = original_adata.X.toarray()
        
    expression_df = pd.DataFrame(
        expression_data.T, 
        index=original_adata.var_names,
        columns=original_adata.obs_names
    )
    expression_df.to_csv("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/3_expression_matrix.csv", index=True, header=True)
    
    # 5. Export metadata CSV (including the new knn_prediction column)
    print("Exporting 2. cell_metadata_merged.csv...")
    original_adata.obs.to_csv("/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/scgpt/output/3_cell_metadata_merged.csv", index=True, header=True)
    
    print("-" * 40)
    print("Step 1 (Data merging and CSV export) complete!")
    print("Output files: expression_matrix.csv and cell_metadata_merged.csv")

except FileNotFoundError as e:
    print(f"File not found error: Please check if the path is correct. Details: {e}")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred during processing: {e}")
    sys.exit(1)