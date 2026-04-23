import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.svm import SVC
#from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, precision_recall_curve, auc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import os
import sys
#sys.stdout.reconfigure(line_buffering=True)

# Load the disease-SNP matrix
disease_snp = scipy.sparse.load_npz("./DiseaseInteractiveNetwork/data_bioinfo/disease_snp.npz").toarray()

# Load SNP and disease lists
snp_list = pd.read_csv("./DiseaseInteractiveNetwork/data_bioinfo/snp_list_rsid.csv", header=None)[0]
disease_list = pd.read_csv("./DiseaseInteractiveNetwork/data_bioinfo/phecode_list.csv", header=None)[0]

# Convert to DataFrame for easier indexing
disease_snp = pd.DataFrame(disease_snp, index=disease_list, columns=snp_list)

# Load disease pairs with labels
disease_pairs_df = pd.read_csv("./disease_pairs_with_labels.csv")

# Create a mapping from disease name to index for fast lookup in the PCA result matrix
disease_to_index = {disease: i for i, disease in enumerate(disease_list)}

# all disease pairs

# -----------------------------------------------------------
# REVISED FEATURE SCALING AND PCA BLOCK
# PCA is done on the disease-SNP matrix to reduce SNP features to 50 components.
# -----------------------------------------------------------

# 1. Standardize the entire disease_snp matrix
# Diseases are samples (rows), SNPs are features (columns).
scaler = StandardScaler()
# Use the underlying numpy array for standardization
disease_snp_scaled = scaler.fit_transform(disease_snp.values)

# 2. PCA for dimensionality reduction (reducing SNP features from N_SNPs to 50)
pca = PCA(n_components=50, random_state=42)
SNP_reduced = pca.fit_transform(disease_snp_scaled)
# SNP_reduced dimensions: (num_diseases, 50)

# 3. Initialize the final feature list and label list
X_reduced, y = [], []

# 4. Extract and concatenate features for each of the selected disease pairs
for _, row in disease_pairs_df.iterrows():
    disease1, disease2 = row["disease pair ID"].split("-")

    if disease1 in disease_to_index and disease2 in disease_to_index:
        disease1_idx = disease_to_index[disease1]
        disease2_idx = disease_to_index[disease2]
    else:
        continue  # Skip if disease names are not found

    label = row["label"]  # Get the label

    # Extract reduced SNP vectors (the 50 PCA components) for both diseases
    snp_reduced_disease1 = SNP_reduced[disease1_idx, :]
    snp_reduced_disease2 = SNP_reduced[disease2_idx, :]

    # Concatenate the reduced features (REVISED STEP)
    # The final feature vector has 100 components (50 from each disease)
    X_pair = np.concatenate((snp_reduced_disease1, snp_reduced_disease2))
    
    # Store the final feature vector and label
    X_reduced.append(X_pair)
    y.append(label)

# Convert to numpy arrays
X_reduced = np.array(X_reduced)
y = np.array(y)

# -----------------------------------------------------------
# NOTE: The variable 'X_reduced' is the final 100-component feature matrix,
# and 'y' is the corresponding label vector.
# -----------------------------------------------------------

'''
## Quantify Variance Reduction 
#variance_retained = np.sum(pca.explained_variance_ratio_)
#variance_lost = 1 - variance_retained

explained_var = pca.explained_variance_ratio_
cumulative_var = np.cumsum(explained_var)
#variance_retained = cumulative_var[-1] * 100
#variance_loss = 100 - variance_retained
variance_retained = cumulative_var[-1] 
variance_loss = 1 - variance_retained

pca_variance_df = pd.DataFrame({
    "Component": np.arange(1, len(explained_var) + 1),
    "ExplainedVarianceRatio": explained_var,
    "CumulativeExplainedVariance": cumulative_var
})
'''

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_reduced, y, test_size=0.2, stratify=y, random_state=42)

# Define models
models = {
    "TrapezoidNN": MLPClassifier(hidden_layer_sizes=(100, 50, 25, 5), activation="relu",
                      solver="adam", max_iter=500, random_state=42)
}

# Prepare output directory
output_dir = "./results"
os.makedirs(output_dir, exist_ok=True)

# Redirect printed output to file
sys.stdout = open(os.path.join(output_dir, "model_resultsRR1SNPvectorsNN.txt"), "w")		

'''
pca_variance_df.to_csv(os.path.join(output_dir, "PCA_Variance_Reduction_SummaryRR1SNPvectors.csv"), index=False)
'''

# Store metrics and curves
pr_curves = {}
roc_auc_scores = {}
all_metrics = []

pr_curves = {}
roc_auc_scores = {}
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
    print(f"\n=== Fold {fold_idx + 1} ===")
    X_train_cv, X_val_cv = X_train[train_idx], X_train[val_idx]
    y_train_cv, y_val_cv = y_train[train_idx], y_train[val_idx]

    for name, model in models.items():
        model.fit(X_train_cv, y_train_cv)
        y_pred = model.predict(X_val_cv)

        try:
            y_prob = model.predict_proba(X_val_cv)[:, 1]
        except AttributeError:
            y_prob = model.decision_function(X_val_cv)
            y_prob = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min())  # normalize

        acc = accuracy_score(y_val_cv, y_pred)
        prec = precision_score(y_val_cv, y_pred, zero_division=0)
        rec = recall_score(y_val_cv, y_pred, zero_division=0)
        f1 = f1_score(y_val_cv, y_pred, zero_division=0)
        roc = roc_auc_score(y_val_cv, y_prob)
        precisions, recalls, _ = precision_recall_curve(y_val_cv, y_prob)
        pr_auc = auc(recalls, precisions)

        cv_results[name]["accuracy"].append(acc)
        cv_results[name]["precision"].append(prec)
        cv_results[name]["recall"].append(rec)
        cv_results[name]["f1"].append(f1)
        cv_results[name]["roc_auc"].append(roc)
        cv_results[name]["pr_auc"].append(pr_auc)


        cv_fold_records.append({
            "Model": name,
            "Fold": fold_idx + 1,
            "Accuracy": acc,
            "Precision": prec,
            "Recall": rec,
            "F1": f1,
            "ROC_AUC": roc,
            "PR_AUC": pr_auc
        })

# ==============================
# Save per-fold CV results
# ==============================
cv_fold_df = pd.DataFrame(cv_fold_records)
cv_fold_df.to_csv(os.path.join(output_dir, "Metrics_RR1SNPvectorsNN.csv"), index=False)


# Evaluate models on test set
for name, model in models.items():
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    try:
        y_prob = model.predict_proba(X_test)[:, 1]
    except AttributeError:
        y_prob = model.decision_function(X_test)
        y_prob = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min())  # normalize

    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, zero_division=0)
    recall = recall_score(y_test, y_pred, zero_division=0)
    f1 = f1_score(y_test, y_pred, zero_division=0)
    roc_auc = roc_auc_score(y_test, y_prob)
    precisions, recalls, _ = precision_recall_curve(y_test, y_prob)
    pr_auc = auc(recalls, precisions)

    pr_curves[name] = (recalls, precisions, pr_auc)
    roc_auc_scores[name] = roc_auc

    print(f"\n=== {name} Model Results RR1 SNP vectors NNtrapezoid ===")
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

all_metrics.append({
        "Model": name,
        "CV_Accuracy_Mean": np.mean(cv_results[name]["accuracy"]),
        "CV_Accuracy_Std": np.std(cv_results[name]["accuracy"]),
        "CV_Precision_Mean": np.mean(cv_results[name]["precision"]),
"CV_precision_Std": np.std(cv_results[name]["precision"]),
        "CV_Recall_Mean": np.mean(cv_results[name]["recall"]),
"CV_recall_Std": np.std(cv_results[name]["recall"]),
"CV_F1_Mean": np.mean(cv_results[name]["f1"]),
"CV_F1_Std": np.std(cv_results[name]["f1"]),
        "CV_ROC_AUC_Mean": np.mean(cv_results[name]["roc_auc"]),
        "CV_ROC_AUC_Std": np.std(cv_results[name]["roc_auc"]),
"CV_PR_AUC_Mean": np.mean(cv_results[name]["pr_auc"]),
"CV_PR_AUC_Std": np.std(cv_results[name]["pr_auc"]),
        "Test_Accuracy": accuracy,
        "Test_Precision": precision,
        "Test_Recall": recall,
        "Test_F1": f1,
        "Test_ROC_AUC": roc_auc,
        "Test_PR_AUC": pr_auc
    })


# Save metrics to CSV
pd.DataFrame(all_metrics).to_csv(os.path.join(output_dir, "model_metricsRR1SNPvectorsNN.csv"), index=False)

