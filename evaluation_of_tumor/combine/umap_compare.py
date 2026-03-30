import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import warnings

warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Liberation Sans', 'sans-serif']
plt.rcParams['pdf.fonttype'] = 42

def get_custom_palette(all_types, highlight_group='Malignant', highlight_color='#D62728'):
    color_map = {}
    base_muted_colors = [
        '#4E79A7', '#A0CBE8', '#59A14F', '#8CD17D', '#B6992D', 
        '#499894', '#86BCB6', '#E15759', '#FF9D9A', '#79706E', 
        '#BAB0AC', '#D37295', '#FABFD2', '#B07AA1', '#D4A6C8',
        '#9D7660', '#D7B5A6', '#728FCE', '#C5D0E6', '#667C26'
    ]
    other_groups = sorted([g for g in all_types if g != highlight_group])
    colors_for_others = []
    if len(other_groups) > 0:
        colors_for_others = base_muted_colors * (len(other_groups) // len(base_muted_colors) + 1)
    
    idx = 0
    for g in all_types:
        if g == highlight_group:
            color_map[g] = highlight_color
        else:
            color_map[g] = colors_for_others[idx]
            idx += 1
    return color_map

def parse_fusion_file(file_path):
    fusion_dict = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                fusion_name = parts[0].strip()
                cell_str = parts[-2].strip() 
                if not cell_str: continue
                cells = [c.strip() for c in cell_str.split(';') if c.strip()]
                fusion_dict[fusion_name] = cells
    except Exception as e:
        print(f"Error parsing file: {e}")
        return {}
    return fusion_dict

def plot_dual_umap_final(h5ad_path1, h5ad_path2, 
                          fusion_file1, fusion_file2,
                          cell_type_key='knn_prediction',
                          title1="ONT", title2="PacBio",
                          figsize=(18, 7), 
                          dot_size=None,
                          save_path=None):
    
    print("Loading data...")
    adata1 = sc.read_h5ad(h5ad_path1)
    adata2 = sc.read_h5ad(h5ad_path2)
    
    for ad in [adata1, adata2]:
        if cell_type_key not in ad.obs.columns:
            raise ValueError(f"Column {cell_type_key} does not exist!")
        ad.obs[cell_type_key] = ad.obs[cell_type_key].astype(str).replace('nan', 'Unknown')

    # --- 2. Determine common fused cells ---
    print("Parsing fusion files...")
    f1_dict = parse_fusion_file(fusion_file1)
    f2_dict = parse_fusion_file(fusion_file2)
    
    common_fusions = set(f1_dict.keys()) & set(f2_dict.keys())
    
    # === Filter specific fusions ===
    filter_target = "GLCCI1:UMAD1"
    if filter_target in common_fusions:
        common_fusions.remove(filter_target)
    
    print(f"Final number of fusion events for plotting: {len(common_fusions)}")
    
    def get_target_cells(fusion_dict, common_keys, all_obs_names):
        target_cells = set()
        dataset_barcodes = set(all_obs_names)
        for fusion in common_keys:
            for cell in fusion_dict[fusion]:
                if cell in dataset_barcodes: target_cells.add(cell)
                elif f"{cell}-1" in dataset_barcodes: target_cells.add(f"{cell}-1")
                elif cell.split('-')[0] in dataset_barcodes: target_cells.add(cell.split('-')[0])
        return target_cells

    fusion_cells_1 = get_target_cells(f1_dict, common_fusions, adata1.obs_names)
    fusion_cells_2 = get_target_cells(f2_dict, common_fusions, adata2.obs_names)

    # --- 3. Generate color map ---
    all_types = sorted(list(set(adata1.obs[cell_type_key].unique()) | 
                            set(adata2.obs[cell_type_key].unique())))
    color_map = get_custom_palette(all_types, highlight_group='Malignant')

    # --- 4. Automatic dot size ---
    if dot_size is None:
        max_cells = max(adata1.shape[0], adata2.shape[0])
        dot_size = 120000 / max_cells
        dot_size = max(2, min(dot_size, 20))

    # --- 5. Plotting ---
    fig = plt.figure(figsize=figsize, dpi=300)
    gs = plt.GridSpec(1, 3, width_ratios=[1, 1, 0.4], wspace=0.1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax_legend = fig.add_subplot(gs[0, 2])

    def plot_layered_umap(ax, adata, fusion_set, title):
        # 1. Distinguish between fusion and non-fusion
        mask_fusion = adata.obs_names.isin(fusion_set)
        adata_normal = adata[~mask_fusion]
        adata_fusion = adata[mask_fusion]
        
        # 2. Draw background: normal cells (circles)
        if adata_normal.shape[0] > 0:
            colors = [color_map[t] for t in adata_normal.obs[cell_type_key]]
            ax.scatter(adata_normal.obsm['X_umap'][:, 0],
                       adata_normal.obsm['X_umap'][:, 1],
                       c=colors,
                       s=dot_size,
                       edgecolor='none',
                       alpha=0.8, # Background slightly transparent to highlight foreground
                       rasterized=True,
                       zorder=1) # Ensure background is at the bottom
            
        '''# 3. Draw foreground: fusion cells (star shape)
        if adata_fusion.shape[0] > 0:
            # Get the corresponding cell type color for fusion cells
            fusion_colors = [color_map[t] for t in adata_fusion.obs[cell_type_key]]
            
            ax.scatter(adata_fusion.obsm['X_umap'][:, 0],
                       adata_fusion.obsm['X_umap'][:, 1],
                       c=fusion_colors,     # Fill with the cell's own color
                       marker='o',          # Change to star (or 'D' for diamond, 'P' for plus)
                       s=dot_size,          # Enlarge size to make it prominent
                       edgecolors='black',  # Black edge
                       linewidths=0.5,      # Edge linewidth
                       alpha=1,             # Opaque
                       label='Fusions',
                       zorder=10,           # Force to top layer
                       rasterized=False)'''

        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.axis('off')

    print("Plotting Dataset 1...")
    plot_layered_umap(ax1, adata1, fusion_cells_1, title1)
    print("Plotting Dataset 2...")
    plot_layered_umap(ax2, adata2, fusion_cells_2, title2)

    # --- 6. Generate legend (modified) ---
    print("Generating legend...")
    ax_legend.axis('off')
    
    legend_handles = []
    
    # 1. Fusion (changed to shape legend)
    # Create a star legend with a black edge
    '''fusion_handle = mlines.Line2D([], [], color='white', 
                                  marker='o', 
                                  markerfacecolor='white', # Inner white color indicates generic shape
                                  markeredgecolor='black', 
                                  markeredgewidth=1.0,
                                  markersize=12, 
                                  label='Fusion Cells')
    legend_handles.append(fusion_handle)'''
    
    # 2. Malignant
    if 'Malignant' in all_types:
        legend_handles.append(mpatches.Patch(facecolor=color_map['Malignant'], edgecolor='none', label='Malignant'))
    
    # 3. Others
    for t in all_types:
        if t != 'Malignant':
            # === Key modification: replace "_" with " " in t ===
            clean_label = t.replace('_', ' ')
            legend_handles.append(mpatches.Patch(facecolor=color_map[t], edgecolor='none', label=clean_label))
            
    ax_legend.legend(handles=legend_handles,
                     loc='center left',
                     #title='Cell Types & Fusions',
                     #title_fontsize=12,
                     fontsize=10,
                     frameon=False,
                     ncol=2 if len(all_types) > 15 else 1)

    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=600)
        print(f"Saved successfully: {save_path}")
    
    plt.show()

# ================= Execution Section (keep your paths unchanged) =================
if __name__ == "__main__":
    # Please ensure the paths here are your actual paths
    h5ad_1 = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/scgpt/3CA_test_embedded_with_fusion.h5ad"
    h5ad_2 = "/home/user/luzhenhao/time/GFHunter/SeekGene_PB_CA/work/result/scgpt/3CA_test_embedded_with_fusion.h5ad"
    
    fusion_file_1 = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/result/SeekGene_CA.fusion_find.tsv"
    fusion_file_2 = "/home/user/luzhenhao/time/GFHunter/SeekGene_PB_CA/work/result/SeekGene_PB_CA.fusion_find.tsv"
    
    out_path = "/home/user/luzhenhao/time/GFHunter/SeekGene_PB_CA/work/code/combine/result/dual_umap.pdf"

    plot_dual_umap_final(
        h5ad_path1=h5ad_1,
        h5ad_path2=h5ad_2,
        fusion_file1=fusion_file_1,
        fusion_file2=fusion_file_2,
        cell_type_key='knn_prediction',
        title1="ONT",
        title2="PacBio",
        save_path=out_path
    )