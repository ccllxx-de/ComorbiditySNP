#W rows of disease vectors, RR1, Concat
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, precision_recall_curve, auc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import os
import sys

# Load the disease-SNP matrix (Assuming file is accessible at this relative path)
disease_snp = scipy.sparse.load_npz("./DiseaseInteractiveNetwork/data_bioinfo/disease_snp.npz").toarray()

# Load SNP and disease lists (Assuming files are accessible)
snp_list = pd.read_csv("./DiseaseInteractiveNetwork/data_bioinfo/snp_list_rsid.csv", header=None)[0]
disease_list = pd.read_csv("./DiseaseInteractiveNetwork/data_bioinfo/phecode_list.csv", header=None)[0]

# Convert to DataFrame for easier indexing
disease_snp = pd.DataFrame(disease_snp, index=disease_list, columns=snp_list)

# Compute Cosine Similarity Matrix (W rows of disease vectors, RR1)
# W has dimensions (num_diseases, num_diseases)
W = metrics.pairwise.cosine_similarity(disease_snp)

# Ensure W is symmetric
assert np.allclose(W, W.T), "Error: W matrix is not symmetric!"

# Load disease pairs with labels (Assuming file is accessible)
disease_pairs_df = pd.read_csv("./disease_pairs_with_labels.csv")

# Convert selected disease pairs to a set for fast lookup
valid_pairs = set(disease_pairs_df["disease pair ID"])

# Create a mapping from disease name to index
disease_to_index = {disease: i for i, disease in enumerate(disease_list)}

# -----------------------------------------------------------
# REVISED FEATURE SCALING AND PCA BLOCK
# PCA is performed on the full W matrix rows (disease feature vectors) once.
# -----------------------------------------------------------

# 1. Standardize the entire W matrix (features for all diseases)
# Each row of W is the feature vector for a disease (similarity to all other diseases).
scaler = StandardScaler()
W_scaled = scaler.fit_transform(W)

# 2. PCA for dimensionality reduction (50 components for each disease vector)
pca = PCA(n_components=50, random_state=42)
W_reduced = pca.fit_transform(W_scaled)
# W_reduced dimensions: (num_diseases, 50)

# 3. Initialize the final feature list and label list
X_final, y_final = [], []

# 4. Extract features for each of the selected disease pairs from the reduced W
for _, row in disease_pairs_df.iterrows():
    disease1, disease2 = row["disease pair ID"].split("-")

    if disease1 in disease_to_index and disease2 in disease_to_index:
        disease1_idx = disease_to_index[disease1]
        disease2_idx = disease_to_index[disease2]
    else:
        continue  # Skip if disease names are not found in the list

    label = row["label"]  # Get the label

    # Extract W_reduced row vectors (the 50 PCA components) for both diseases
    w_reduced_disease1 = W_reduced[disease1_idx, :]
    w_reduced_disease2 = W_reduced[disease2_idx, :]

    # Concatenate the reduced features (REVISED STEP)
    # The resulting vector has 100 components (50 from each disease)
    X_pair = np.concatenate((w_reduced_disease1, w_reduced_disease2))
    
    # Store the final feature vector and label
    X_final.append(X_pair)
    y_final.append(label)

# Convert to numpy arrays
X_reduced = np.array(X_final)
y = np.array(y_final)

# -----------------------------------------------------------
# NOTE: The variable 'X_reduced' is the final 100-component feature matrix,
# replacing the need for the original 'X_scaled' variable in subsequent steps.
# 'y' is the corresponding label vector.
# -----------------------------------------------------------



# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_reduced, y, test_size=0.2, stratify=y, random_state=42)

# Define models
models = {
    "Random Forest": RandomForestClassifier(n_estimators=100, random_state=42),
    #"SVM": SVC(kernel="rbf", C=1.0, probability=True, random_state=42),
    #"Logistic Regression": LogisticRegression(C=1.0, n_jobs=-1, solver="saga",max_iter=200, random_state=42),
    #"Neural Network": MLPClassifier(hidden_layer_sizes=(100,), activation="relu", solver="adam", max_iter=500, random_state=42)
}

# Prepare output directory
output_dir = "./results"
os.makedirs(output_dir, exist_ok=True)

# Redirect printed output to file
sys.stdout = open(os.path.join(output_dir, "Concat_model_resultsCV10alldiseasepairRR1diseaserowvectorsPCA100TrainTest.txt"), "w")

'''
pca_variance_df.to_csv(os.path.join(output_dir, "PCA_Variance_Reduction_SummaryRR1diseaserowvectorsCV10alldiseasepairSDVarPCA100TrainTest.csv"), index=False)
'''

# Store metrics and curves
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

# NEW: to store all fold results for CSV export
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

# NEW: record per-fold data for CSV export
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
cv_fold_df.to_csv(os.path.join(output_dir, "Concat_CV10Fold_Metrics_PerModelSDVarCV10RR1diseaserowvectorsalldiseasepairPCA100TrainTest.csv"), index=False)



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


    print(f"\n=== {name} Model Results RR1 diseaserowvectors CV10 PCA100 alldiseasepair ===")
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
pd.DataFrame(all_metrics).to_csv(os.path.join(output_dir, "Concat_all_model_metricsCV10RR1diseaserowvectorsalldiseasepairPCA100TrainTest.csv"), index=False)

# Save PR curve plot
plt.figure(figsize=(8, 6))
for name, (recalls, precisions, pr_auc) in pr_curves.items():
    plt.plot(recalls, precisions, marker='.', label=f'{name} (PR AUC = {pr_auc:.4f})')

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve PCA100 RR1 diseaserowvectors CV10 alldiseasepair")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pr_curveRR1diseaserowvectorsCV10alldiseasepairTrainTestConcat.png"))

# Close stdout redirection
#sys.stdout.close()
