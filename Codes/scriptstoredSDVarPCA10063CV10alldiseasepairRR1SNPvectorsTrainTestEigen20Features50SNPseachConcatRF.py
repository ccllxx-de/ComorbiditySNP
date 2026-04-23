import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, precision_recall_curve, auc
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
# PCA AND FEATURE CREATION
# -----------------------------------------------------------

# 1. Standardize the entire disease_snp matrix (SNPs are features)
scaler = StandardScaler()
disease_snp_scaled = scaler.fit_transform(disease_snp.values)

# 2. PCA for dimensionality reduction (reducing SNP features from N_SNPs to 50)
pca = PCA(n_components=50, random_state=42)
SNP_reduced = pca.fit_transform(disease_snp_scaled)
# SNP_reduced dimensions: (num_diseases, 50)

# 3. Initialize the final feature list and label list
X_final, y_final = [], []

# 4. Extract and concatenate features
for _, row in disease_pairs_df.iterrows():
    disease1, disease2 = row["disease pair ID"].split("-")

    if disease1 in disease_to_index and disease2 in disease_to_index:
        disease1_idx = disease_to_index[disease1]
        disease2_idx = disease_to_index[disease2]
    else:
        continue 

    label = row["label"] 

    snp_reduced_disease1 = SNP_reduced[disease1_idx, :]
    snp_reduced_disease2 = SNP_reduced[disease2_idx, :]

    # Concatenate the reduced features (100 components total)
    X_pair = np.concatenate((snp_reduced_disease1, snp_reduced_disease2))
    
    X_final.append(X_pair)
    y_final.append(label)

X_reduced = np.array(X_final)
y = np.array(y_final)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_reduced, y, test_size=0.2, stratify=y, random_state=42)

# Define models
models = {
    "Random Forest": RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1),
    #"Neural Network": MLPClassifier(hidden_layer_sizes=(100,), activation="relu", solver="adam", max_iter=500, random_state=42)
}

# Prepare output directory
output_dir = "./results"
os.makedirs(output_dir, exist_ok=True)

# -----------------------------------------------------------
# 5. FEATURE IMPORTANCE ANALYSIS AND SNP EXTRACTION
# -----------------------------------------------------------

RF_model = models["Random Forest"]
RF_model.fit(X_train, y_train)

# Get feature importance
importances = RF_model.feature_importances_
feature_indices = np.argsort(importances)[::-1] # Indices sorted by importance (descending)
top_n_features = 20
top_20_indices = feature_indices[:top_n_features]
top_20_importances = importances[top_20_indices]

# Map the 100 feature indices back to their PCA component index (0-49) and disease (D1 or D2)
analysis_records = []
MAX_SNPS_TO_EXTRACT = 1000

# Combined SNP importance score array (size = N_SNPS)
snp_loading_scores = np.zeros(N_SNPS) 
pca_component_contributions = {} 

for idx in top_20_indices:
    importance = importances[idx]
    
    # Determine which PC (0-49) and which disease (D1 or D2) the feature corresponds to
    if idx < 50:
        component_idx = idx 
        disease_key = "Disease 1"
    else:
        component_idx = idx - 50 
        disease_key = "Disease 2"

    # Get the corresponding eigenvector (SNP loadings)
    eigenvector = pca.components_[component_idx]
    eigenvalue = pca.explained_variance_[component_idx]

    # Accumulate SNP importance: the contribution of a SNP to the classification is
    # proportional to its loading on the PC, weighted by the PC's importance for the task.
    # We use the absolute loading since both positive and negative correlation indicate influence.
    snp_loading_scores += np.abs(eigenvector) * importance
    
    # Track contribution to PC analysis file
    if component_idx not in pca_component_contributions:
        pca_component_contributions[component_idx] = 0
    pca_component_contributions[component_idx] += importance

    # Rank the vector components (loadings) by absolute value for this specific PC
    loading_indices_sorted = np.argsort(np.abs(eigenvector))[::-1]
    top_50_snp_indices = loading_indices_sorted[:50]
    top_50_snps = snp_list.iloc[top_50_snp_indices].tolist()

    analysis_records.append({
        "Feature_Index": idx,
        "Disease_Key": disease_key,
        "PC_Index (0-49)": component_idx,
        "Importance": importance,
        "Eigenvalue": eigenvalue,
        "Top_50_SNP_Indices": top_50_snp_indices.tolist(),
        "Top_50_SNPs": top_50_snps,
    })

# Identify the top 1000 SNPs based on the combined weighted loading score
top_1000_snp_indices = np.argsort(snp_loading_scores)[::-1][:MAX_SNPS_TO_EXTRACT]
top_1000_snps_list = snp_list.iloc[top_1000_snp_indices].tolist()

# -----------------------------------------------------------
# 6. EXECUTION AND SAVING RESULTS
# -----------------------------------------------------------

# Redirect printed output to file
sys.stdout = open(os.path.join(output_dir, "Concat_model_results_SNP_PCA_Importance.txt"), "w")

# Store metrics and curves
pr_curves = {}
all_metrics = []

# Set up 10-fold Stratified CV
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# Initialize per-model CV storage
cv_results = {name: {
    "accuracy": [], "precision": [], "recall": [],
    "f1": [], "roc_auc": [], "pr_auc": []
} for name in models}

cv_fold_records = []

# Perform manual 10-fold CV
for fold_idx, (train_idx, val_idx) in enumerate(cv.split(X_train, y_train)):
    print(f"\n=== CV Fold {fold_idx + 1} ===")
    X_train_cv, X_val_cv = X_train[train_idx], X_train[val_idx]
    y_train_cv, y_val_cv = y_train[train_idx], y_train[val_idx]

    for name, model in models.items():
        model.fit(X_train_cv, y_train_cv)
        y_pred = model.predict(X_val_cv)

        try:
            y_prob = model.predict_proba(X_val_cv)[:, 1]
        except AttributeError:
            y_prob = model.decision_function(X_val_cv)
            y_prob = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min()) 

        # Calculate CV metrics
        acc = accuracy_score(y_val_cv, y_pred)
        prec = precision_score(y_val_cv, y_pred, zero_division=0)
        rec = recall_score(y_val_cv, y_pred, zero_division=0)
        f1 = f1_score(y_val_cv, y_pred, zero_division=0)
        roc = roc_auc_score(y_val_cv, y_prob)
        precisions, recalls, _ = precision_recall_curve(y_val_cv, y_prob)
        pr_auc = auc(recalls, precisions)

        # Store CV metrics
        cv_results[name]["accuracy"].append(acc)
        cv_results[name]["precision"].append(prec)
        cv_results[name]["recall"].append(rec)
        cv_results[name]["f1"].append(f1)
        cv_results[name]["roc_auc"].append(roc)
        cv_results[name]["pr_auc"].append(pr_auc)

        cv_fold_records.append({
            "Model": name, "Fold": fold_idx + 1, "Accuracy": acc,
            "Precision": prec, "Recall": rec, "F1": f1, "ROC_AUC": roc, "PR_AUC": pr_auc
        })

# Save per-fold CV results
cv_fold_df = pd.DataFrame(cv_fold_records)
cv_fold_df.to_csv(os.path.join(output_dir, "Concat_CV10Fold_Metrics_SNP_PCA_Detailed_RF_RR1.csv"), index=False)

# Evaluate models on test set
final_metrics_summary = []
for name, model in models.items():
    model.fit(X_train, y_train) # Retrain on full training set
    y_pred = model.predict(X_test)

    try:
        y_prob = model.predict_proba(X_test)[:, 1]
    except AttributeError:
        y_prob = model.decision_function(X_test)
        y_prob = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min()) 

    # Calculate Test metrics
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, zero_division=0)
    recall = recall_score(y_test, y_pred, zero_division=0)
    f1 = f1_score(y_test, y_pred, zero_division=0)
    roc_auc = roc_auc_score(y_test, y_prob)
    precisions, recalls, _ = precision_recall_curve(y_test, y_prob)
    pr_auc = auc(recalls, precisions)

    pr_curves[name] = (recalls, precisions, pr_auc)

    print(f"\n=== {name} Model Results (SNP PCA 50, Concatenated) ===")
    print(f"CV Accuracy: {np.mean(cv_results[name]['accuracy']):.4f} ± {np.std(cv_results[name]['accuracy']):.4f}")
    print(f"CV Precision: {np.mean(cv_results[name]['precision']):.4f} ± {np.std(cv_results[name]['precision']):.4f}")
    print(f"CV Recall: {np.mean(cv_results[name]['recall']):.4f} ± {np.std(cv_results[name]['recall']):.4f}")
    print(f"CV F1 Score: {np.mean(cv_results[name]['f1']):.4f} ± {np.std(cv_results[name]['f1']):.4f}")
    print(f"CV ROC AUC: {np.mean(cv_results[name]['roc_auc']):.4f} ± {np.std(cv_results[name]['roc_auc']):.4f}")
    print(f"CV PR AUC: {np.mean(cv_results[name]['pr_auc']):.4f} ± {np.std(cv_results[name]['pr_auc']):.4f}")
    print(f"Test Accuracy: {accuracy:.4f}")
    print(f"Test Precision: {precision:.4f}")
    print(f"Test Recall: {recall:.4f}")
    print(f"Test F1 Score: {f1:.4f}")
    print(f"Test ROC AUC: {roc_auc:.4f}")
    print(f"Test PR AUC: {pr_auc:.4f}\n")

    # Store metrics for final CSV
    final_metrics_summary.append({
        "Model": name, 
        "CV_Accuracy_Mean": np.mean(cv_results[name]["accuracy"]), "CV_Accuracy_Std": np.std(cv_results[name]["accuracy"]),
        "Test_Accuracy": accuracy,
        "Test_ROC_AUC": roc_auc,
        "Test_PR_AUC": pr_auc
        # Other metrics omitted for brevity in summary, but calculated above
    })


# Save the detailed PC analysis results (Top 20 RF Features to 50 most loaded SNPs)
pd.DataFrame(analysis_records).to_csv(os.path.join(output_dir, "RF_Top20_PC_to_SNP_Analysis.csv"), index=False)

# Save the top 1000 SNPs list
pd.Series(top_1000_snps_list).to_csv(os.path.join(output_dir, "Top_1000_SNPs_Ranked.csv"), index=False, header=['SNP_ID'])

# Add the top 1000 SNPs to the final metrics file
metrics_df = pd.DataFrame(final_metrics_summary)
metrics_df['Top_1000_SNPs'] = ', '.join(top_1000_snps_list) 
metrics_df.to_csv(os.path.join(output_dir, "Concat_all_model_metrics_SNP_PCA_Final.csv"), index=False)

# Print the final summary of the feature importance and SNP extraction
print("\n" + "="*80)
print("=== FEATURE IMPORTANCE AND SNP RANKING SUMMARY ===")
print("="*80)
print(f"Total Original SNPs: {N_SNPS}")
print(f"Total Principal Components Used: 50")
print("Top 20 Most Important Features (PCA Components) from Random Forest:")
for i in range(top_n_features):
    idx = top_20_indices[i]
    if idx < 50:
        d_key, pc_idx = "Disease 1", idx
    else:
        d_key, pc_idx = "Disease 2", idx - 50
    print(f"Rank {i+1}: Feature Index {idx} (PC_{pc_idx} for {d_key}) | Importance: {importances[idx]:.4f}")

print(f"\nTop {MAX_SNPS_TO_EXTRACT} Most Relevant SNPs (Ranked by Weighted PC Loading):")
print(top_1000_snps_list[:5] + ["...", f"({MAX_SNPS_TO_EXTRACT} total)"])


# Save PR curve plot
plt.figure(figsize=(8, 6))
for name, (recalls, precisions, pr_auc) in pr_curves.items():
    plt.plot(recalls, precisions, marker='.', label=f'{name} (PR AUC = {pr_auc:.4f})')

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve SNP PCA 50 (Concatenated Features)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pr_curve_SNP_PCA_Concat.png"))

# Close stdout redirection
sys.stdout.close()
