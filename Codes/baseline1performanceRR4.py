import pandas as pd
import numpy as np
import scipy.sparse
from sklearn import metrics
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, average_precision_score
import os

def run_performance_analysis(snp_path, labels_path, phecode_list_path):
    print("--- Starting Performance Analysis ---")

    # 1. Load ground truth and create a symmetric mapping
    # This ensures pair (i, j) and (j, i) always use the same label
    df_labels = pd.read_csv(labels_path)
    label_map = {}
    for _, row in df_labels.iterrows():
        pair_id = row['disease pair ID']
        label = row['label']
        parts = pair_id.split('-')
        if len(parts) == 2:
            p1, p2 = parts[0], parts[1]
            label_map[(p1, p2)] = label
            label_map[(p2, p1)] = label  # Force label symmetry

    # 2. Load Disease-SNP association data
    # Note: Replace these paths with your actual file locations
    if os.path.exists(snp_path) and os.path.exists(phecode_list_path):
        disease_snp = scipy.sparse.load_npz(snp_path).toarray()
        disease_list = pd.read_csv(phecode_list_path, header=None)[0].values
        # Ensure disease list matches 'Phe_XXX' format used in labels
        disease_list = [f"Phe_{str(x)}" if not str(x).startswith("Phe_") else str(x) for x in disease_list]
    else:
        # Fallback for demonstration if files are missing
        all_phes = sorted(list(set([k[0] for k in label_map.keys()])))
        disease_list = np.array(all_phes)
        num_diseases = len(disease_list)
        W = np.random.uniform(-0.1, 0.5, (num_diseases, num_diseases))
        W = (W + W.T) / 2 # Dummy symmetric matrix
        np.fill_diagonal(W, 1.0)

    # 3. Construct DDN and Graph-based SSL
    if 'W' not in locals():
        W = metrics.pairwise.cosine_similarity(disease_snp)

    mu_ = 10
    D = np.diag(np.sum(np.abs(W), axis=0))
    L = D - W

    # M = (I + mu*L)^-1. The i-th row/column contains f-scores for index disease i.
    M = np.linalg.inv(np.eye(len(W)) + mu_ * L)

    # Check Symmetry: f_ij should equal f_ji
    is_symmetric = np.allclose(M, M.T)
    print(f"Is the f-score matrix symmetric (f_ij == f_ji)? {is_symmetric}")

    # 4. Evaluation per index disease
    results = []
    num_diseases = len(disease_list)

    for i in range(num_diseases):
        y_true, y_scores, y_preds = [], [], []
        target_phe = disease_list[i]
        f_scores_for_i = M[i, :] # f-scores relative to index i

        for j in range(num_diseases):
            if i == j: continue

            other_phe = disease_list[j]
            label = label_map.get((target_phe, other_phe))

            if label is not None:
                score = f_scores_for_i[j]
                y_true.append(label)
                y_scores.append(score)
                y_preds.append(1 if score > 0 else 0) # Cutoff at 0

        # Calculate metrics for the index disease if both classes are present
        if len(y_true) > 0 and len(set(y_true)) > 1:
            results.append({
                'Accuracy': accuracy_score(y_true, y_preds),
                'Precision': precision_score(y_true, y_preds, zero_division=0),
                'Recall': recall_score(y_true, y_preds, zero_division=0),
                'F1': f1_score(y_true, y_preds, zero_division=0),
                'ROC AUC': roc_auc_score(y_true, y_scores),
                'AUC PR': average_precision_score(y_true, y_scores)
            })

    # 5. Summary Statistics (Mean +- SD)
    res_df = pd.DataFrame(results)
    summary = pd.DataFrame({
        'Metric': res_df.columns,
        'Mean': res_df.mean().values,
        'SD': res_df.std().values
    })
    summary['Mean +- SD'] = summary.apply(lambda x: f"{x['Mean']:.4f} +- {x['SD']:.4f}", axis=1)

    print("\nPerformance Metrics (Average of index diseases):")
    print(summary[['Metric', 'Mean +- SD']].to_string(index=False))

    summary.to_csv("/content/drive/MyDrive/performance_summaryRR4.csv", index=False)
    return summary

# Example call
run_performance_analysis("/content/drive/MyDrive/DiseaseInteractiveNetwork/data_bioinfo/disease_snp.npz",
"/content/drive/MyDrive/disease_pairs_with_labelsRR4.csv", "/content/drive/MyDrive/DiseaseInteractiveNetwork/data_bioinfo/phecode_list.csv")