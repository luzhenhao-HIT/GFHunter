import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import glob
import numpy as np
import matplotlib.patches as mpatches

# ================= Configuration Paths =================
DATA_DIR = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output/Step4_Plot_Final_Fixed/Plot_Data"
OUT_DIR = "/home/user/luzhenhao/time/GFHunter/SeekGene_CA/work/code/last_analysis/output/Step5_Nature_Grid_Bright_Aligned" # Modify output directory name to distinguish

# Classification logic
PATHWAY_GROUPS = {
    "Proliferative_Signaling": [
        "03_MYC_TARGETS_V1", "05_MYC_TARGETS_V2", "06_E2F_TARGETS", 
        "17_PI3K_AKT_MTOR_SIGNALING", "08_MTORC1_SIGNALING"
    ],
    "Cell_Cycle": [
        "09_G2M_CHECKPOINT", "11_MITOTIC_SPINDLE"
    ],
    "Metabolic_Reprogramming": [
        "14_GLYCOLYSIS", "02_OXIDATIVE_PHOSPHORYLATION", "12_FATTY_ACID_METABOLISM"
    ],
    "Stress_Adaptive": [
        "01_UNFOLDED_PROTEIN_RESPONSE", "04_PROTEIN_SECRETION", "07_DNA_REPAIR"
    ],
    "Others": [] 
}

for subfolder in PATHWAY_GROUPS.keys():
    os.makedirs(os.path.join(OUT_DIR, subfolder), exist_ok=True)

# ================= Visual Style Settings =================
label_map = {
    "Fusion-Positive": "F+",
    "Non-Fusion": "F-"
}

my_palette = {
    "F+": "#E64B35",  
    "F-": "#448AFF"   
}

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 12,
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': False,
    'axes.spines.bottom': False,
})

# ================= 1. Pre-calculate global Y-axis range =================
print("Scanning all files to calculate a unified Y-axis range...")
csv_files = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))

if not csv_files:
    print("No CSV files found.")
    exit()

all_scores = []
for csv_file in csv_files:
    df_temp = pd.read_csv(csv_file, index_col=0)
    all_scores.extend(df_temp['Score'].values)

# Calculate global min/max and leave a 5% buffer space
y_min_raw = min(all_scores)
y_max_raw = max(all_scores)
y_range = y_max_raw - y_min_raw
GLOBAL_Y_MIN = y_min_raw - (y_range * 0.05)
GLOBAL_Y_MAX = y_max_raw + (y_range * 0.05)

print(f"Global Y-axis range set: [{GLOBAL_Y_MIN:.2f}, {GLOBAL_Y_MAX:.2f}]")

# ================= 2. Main Plotting Logic =================
for csv_file in csv_files:
    # Prepare data
    filename = os.path.basename(csv_file).replace(".csv", "")
    subset = pd.read_csv(csv_file, index_col=0)
    subset['Status'] = subset['Status'].map(label_map)
    
    # Determine folder
    target_folder = "Others"
    for group, p_list in PATHWAY_GROUPS.items():
        if any(pid in filename for pid in p_list):
            target_folder = group
            break
    save_dir = os.path.join(OUT_DIR, target_folder)

    # Plotting
    fig, ax = plt.subplots(figsize=(3.5, 4.5))

    # Grid
    ax.grid(axis='both', linestyle='-', alpha=0.6, color='#E0E0E0', zorder=0)

    # A. Violin plot
    sns.violinplot(
        data=subset,
        x='Status',
        y='Score',
        hue='Status',
        palette=my_palette,
        inner=None,
        linewidth=0,
        density_norm='width',
        cut=0,
        alpha=0.9,
        zorder=1,
        ax=ax,
        legend=False
    )

    # B. Box plot
    sns.boxplot(
        data=subset,
        x='Status',
        y='Score',
        width=0.15,
        boxprops={'facecolor': '#2C3E50', 'edgecolor': 'black', 'linewidth': 1, 'zorder': 2},
        medianprops={'color': 'white', 'linewidth': 2, 'zorder': 3},
        whiskerprops={'color': 'black', 'linewidth': 1.5, 'zorder': 2},
        capprops={'color': 'black', 'linewidth': 1.5, 'zorder': 2},
        showfliers=False,
        zorder=2,
        ax=ax
    )

    # C. Scatter plot
    np.random.seed(42)
    sns.stripplot(
        data=subset,
        x='Status',
        y='Score',
        color='black',
        jitter=0.2,
        size=2.5,
        alpha=0.5,
        linewidth=0,
        zorder=3,
        ax=ax
    )

    # Embellishments
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    # Set unified Y-axis range [Core modification]
    ax.set_ylim(GLOBAL_Y_MIN, GLOBAL_Y_MAX)

    # Label formatting
    labels = [item.get_text() for item in ax.get_xticklabels()]
    clean_labels = [label.replace(' (', '\n(') for label in labels] # There are actually no parentheses here, but keeping the logic is fine
    # Force F+ / F- to be non-italic
    ax.set_xticklabels(clean_labels, rotation=0, fontsize=11, fontweight='medium')
    
    # Remove tick marks
    ax.tick_params(axis='both', which='both', length=0, pad=5)
    
    # Save
    save_path = os.path.join(save_dir, filename)
    plt.savefig(save_path + ".png", dpi=600, bbox_inches='tight')
    plt.savefig(save_path + ".pdf", bbox_inches='tight')
    
    plt.close(fig)

# ================= Generate Legend =================
print("Generating legend...")
fig_leg, ax_leg = plt.subplots(figsize=(4, 1))

legend_patches = [
    mpatches.Patch(color='#E64B35', label='F+'),
    mpatches.Patch(color='#448AFF', label='F-')
]

ax_leg.legend(
    handles=legend_patches, 
    loc='center', 
    ncol=2, 
    frameon=False, 
    fontsize=12
)
ax_leg.axis('off')

leg_save_path = os.path.join(OUT_DIR, "Legend_Bright")
plt.savefig(leg_save_path + ".png", dpi=600, bbox_inches='tight')
plt.savefig(leg_save_path + ".pdf", bbox_inches='tight')
plt.close(fig_leg)

print("=== Plotting Completed ===")
print(f"All images have aligned Y-axes ({GLOBAL_Y_MIN:.2f} ~ {GLOBAL_Y_MAX:.2f})")
print(f"Images saved to: {OUT_DIR}")