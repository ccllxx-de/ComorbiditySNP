import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import os
import sys

# --- File Loading ---
# Load the disease-SNP matrix
disease_snp_sparse = scipy.sparse.load_npz("./DiseaseInteractiveNetwork/data_bioinfo/disease_snp.npz")
disease_snp_array = disease_snp_sparse.toarray()

# Load SNP and disease lists
snp_list = pd.read_csv("./DiseaseInteractiveNetwork/data_bioinfo/snp_list_rsid.csv", header=None)[0]
disease_list = pd.read_csv("./DiseaseInteractiveNetwork/data_bioinfo/phecode_list.csv", header=None)[0]

# Convert to DataFrame for easier indexing
disease_snp = pd.DataFrame(disease_snp_array, index=disease_list, columns=snp_list)
N_SNPS = disease_snp.shape[1] # Total number of original SNP features

# Load disease pairs with labels
disease_pairs_df = pd.read_csv("./disease_pairs_with_labels.csv")

# Create a mapping from disease name to index
disease_to_index = {disease: i for i, disease in enumerate(disease_list)}


# -----------------------------------------------------------
# PCA AND VARIANCE REDUCTION CALCULATION
# -----------------------------------------------------------

# 1. Standardize the entire disease_snp matrix
scaler = StandardScaler()
disease_snp_scaled = scaler.fit_transform(disease_snp.values)

# 2. PCA for dimensionality reduction
n_components = 50
pca = PCA(n_components=n_components, random_state=42)
SNP_reduced = pca.fit_transform(disease_snp_scaled)

# 3. Calculate Variance Metrics
total_variance_original = disease_snp_scaled.shape[1] 

# Explained variance by the 50 components
explained_variance_sum = np.sum(pca.explained_variance_)
explained_variance_ratio_sum = np.sum(pca.explained_variance_ratio_)

# Variance Reduction (Information Loss)
variance_loss_absolute = total_variance_original - explained_variance_sum
variance_loss_percentage = (1 - explained_variance_ratio_sum) * 100


# Prepare output directory
output_dir = "./results"
os.makedirs(output_dir, exist_ok=True)


# 4. Save the Variance Info to a file
output_file_path = os.path.join(output_dir, "PCA_Variance_Reduction_RR1SNPvectors.txt")

with open(output_file_path, "w") as f:
    f.write("PCA Dimensionality Reduction Summary RR1 SNP-V \n")
    f.write("====================================\n")
    f.write(f"Original SNP Features: {total_variance_original}\n")
    f.write(f"Reduced PCA Components: {n_components}\n")
    f.write("-" * 30 + "\n")
    f.write(f"Total Explained Variance (Sum): {explained_variance_sum:.4f}\n")
    f.write(f"Explained Variance Ratio: {explained_variance_ratio_sum:.4f} ({explained_variance_ratio_sum*100:.2f}%)\n")
    f.write("-" * 30 + "\n")
    f.write(f"Variance Reduced (Information Loss): {variance_loss_absolute:.4f}\n")
    f.write(f"Percentage of Information Lost: {variance_loss_percentage:.2f}%\n")

print(f"Variance reduction analysis saved to: {output_file_path}")
