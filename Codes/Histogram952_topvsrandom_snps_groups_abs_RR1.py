import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------
# 1. Load and clean data
# -----------------------------
file_name = "/content/drive/MyDrive/952_topvsrandom_snps_groups_abs.csv"
df = pd.read_csv(file_name)

# Clean column names
df.columns = ["Index", "SNP", "Abs_Effect_Magnitude", "Group"]
df = df.drop(columns=["Index"])

# -----------------------------
# 2. Compute MEAN absolute effect per unique SNP (one value per SNP)
# -----------------------------
mean_per_snp = (
    df.groupby(['SNP', 'Group'])['Abs_Effect_Magnitude']
      .mean()
      .reset_index()
)

# Split into two groups
random_means = mean_per_snp[mean_per_snp['Group'] == 'Random SNPs']['Abs_Effect_Magnitude']
top_means    = mean_per_snp[mean_per_snp['Group'] == 'Target/Top SNPs']['Abs_Effect_Magnitude']

# -----------------------------
# 3. Apply capping: all values >= 20 → 21 (for binning)
# -----------------------------
CAP_THRESHOLD = 20
CAP_VALUE = 21  # Represents the "≥20" bin

random_capped = np.where(random_means > CAP_THRESHOLD, CAP_VALUE, random_means)
top_capped    = np.where(top_means > CAP_THRESHOLD, CAP_VALUE, top_means)

# -----------------------------
# 4. Define bins: 0–1, 1–2, ..., 19–20, 20–21 (where 21 = ≥20)
# -----------------------------
bins = np.arange(0, CAP_THRESHOLD + 2)  # 0,1,2,...,20,21 → creates bins up to [20,21]

# -----------------------------
# 5. Plot step histograms (outlines only)
# -----------------------------
plt.figure(figsize=(11, 6.5))

# Top SNPs - blue
plt.hist(top_capped, bins=bins, histtype='step', color='#d62728',
         label=f'Target/Top SNPs (n=952)', linewidth=2.8, zorder=3)

# Random SNPs - red/orange
plt.hist(random_capped, bins=bins, histtype='step', color='#1f77b4',
         label=f'Random SNPs (n=952)', linewidth=2.8, zorder=3)

# -----------------------------
# 6. Customize x-axis: label last tick as "≥20"
# -----------------------------
tick_positions = np.arange(0, CAP_THRESHOLD + 1, 1)  # 0,2,4,...,20
tick_labels = [str(int(x)) for x in tick_positions]
tick_labels[-1] = '≥20'  # Replace 20 with "≥20"

plt.xticks(tick_positions, tick_labels, fontsize=11)

# Optional: add minor ticks for all integers
plt.gca().set_xticks(np.arange(0, 22, 1), minor=True)
plt.gca().tick_params(axis='x', which='minor', length=0)

# -----------------------------
# 7. Final plot styling
# -----------------------------
plt.xlim(0, CAP_THRESHOLD + 1)
plt.ylim(0, None)

plt.title("Distribution of Mean Absolute Effect per SNP\n"
          "952 Top vs 952 Random SNPs", fontsize=15, pad=20)
plt.xlabel("Mean Absolute Effect Magnitude ", fontsize=13)
plt.ylabel("Number of SNPs", fontsize=13)

plt.legend(fontsize=12, loc='upper right')
plt.grid(axis='y', alpha=0.4, linestyle='--', zorder=0)

# Add count of extreme SNPs in corner
n_top_extreme = (top_means >= 20).sum()
n_rand_extreme = (random_means >= 20).sum()

text = f"SNPs with mean ≥ 20:\nTop: {n_top_extreme}\nRandom: {n_rand_extreme}"
plt.text(0.68, 0.75, text, transform=plt.gca().transAxes,
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.8'),
         fontsize=12, verticalalignment='top')

plt.tight_layout()

# -----------------------------
# 8. Save the plot
# -----------------------------
output_plot_file = "/content/drive/MyDrive/mean_abs_effect_histogram_top_vs_random_cap20.png"
plt.savefig(output_plot_file, dpi=300, bbox_inches='tight')
print(f"\nPlot saved to: {output_plot_file}")

plt.show()