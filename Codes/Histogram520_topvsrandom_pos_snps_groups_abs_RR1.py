import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load data
file_path = "/content/drive/MyDrive/topandrandom_pos_means_df.csv"
df = pd.read_csv(file_path)

# 2. Extract means for each group (already mean values in this dataset)
# Using .abs() to ensure absolute magnitude as requested
top_means = df[df['Group'] == 'Top Positive SNPs']['Mean_Positive_Effect'].abs()
random_means = df[df['Group'] == 'Random Positive SNPs']['Mean_Positive_Effect'].abs()

n_top = len(top_means)
n_random = len(random_means)

# -----------------------------
# 3. Apply capping: all values >= 20 -> 21 (for binning)
# -----------------------------
CAP_THRESHOLD = 20
CAP_VALUE = 21  # Represents the ">=20" bin

random_capped = np.where(random_means > CAP_THRESHOLD, CAP_VALUE, random_means)
top_capped = np.where(top_means > CAP_THRESHOLD, CAP_VALUE, top_means)

# -----------------------------
# 4. Define bins: 0-1, 1-2, ..., 19-20, 20-21 (where 21 = >=20)
# -----------------------------
bins = np.arange(0, CAP_THRESHOLD + 2)  # 0,1,2,...,20,21

# -----------------------------
# 5. Plot step histograms (outlines only)
# -----------------------------
plt.figure(figsize=(11, 6.5))

# Top SNPs - color #d62728
plt.hist(top_capped, bins=bins, histtype='step', color='#d62728',
         label=f'Top Positive SNPs (n={n_top})', linewidth=2.8, zorder=3)

# Random SNPs - color #1f77b4
plt.hist(random_capped, bins=bins, histtype='step', color='#1f77b4',
         label=f'Random Positive SNPs (n={n_random})', linewidth=2.8, zorder=3)

# -----------------------------
# 6. Customize x-axis: label last tick as ">=20"
# -----------------------------
tick_positions = np.arange(0, CAP_THRESHOLD + 1, 1)
tick_labels = [str(int(x)) for x in tick_positions]
tick_labels[-1] = '≥20'

plt.xticks(tick_positions, tick_labels, fontsize=11)
plt.gca().set_xticks(np.arange(0, 22, 1), minor=True)
plt.gca().tick_params(axis='x', which='minor', length=0)

# -----------------------------
# 7. Final plot styling
# -----------------------------
plt.xlim(0, CAP_THRESHOLD + 1)
plt.ylim(0, None)

plt.title(f"Distribution of Mean Absolute Effect per SNP\n"
          f"{n_top} Top vs {n_random} Random Positive SNPs", fontsize=15, pad=20)
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
plt.savefig("/content/drive/MyDrive/pos_snps_absolute_effect_histogram.png")
plt.show()