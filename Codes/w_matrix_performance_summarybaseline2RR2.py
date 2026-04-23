import pandas as pd
import numpy as np
import scipy.sparse
from sklearn import metrics
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, average_precision_score
import os

def run_w_matrix_performance(snp_path, labels_path, phecode_list_path):
    print("--- Starting W-Matrix Performance Analysis ---")

    # 1. Load ground truth labels and ensure symmetry
    # This creates a map for both (A, B) and (B, A) to the same label.
    df_labels = pd.read_csv(labels_path)
    label_map = {}
    for _, row in df_labels.iterrows():
        p1, p2 = row['disease pair ID'].split('-')
        label_map[(p1, p2)] = row['label']
        label_map[(p2, p1)] = row['label']

    # 2. Load SNP data and PheCode list
    # Assuming the files are named as in your previous context
    disease_snp = scipy.sparse.load_npz(snp_path).toarray()
    disease_list = pd.read_csv(phecode_list_path, header=None)[0].values
    # Ensure IDs are strings and prefixed correctly for label matching
    disease_list = [f"Phe_{str(x)}" if not str(x).startswith("Phe_") else str(x) for x in disease_list]
    num_diseases = len(disease_list)

    # 3. Compute Cosine Similarity Matrix W
    print("Computing Cosine Similarity Matrix W...")
    W = metrics.pairwise.cosine_similarity(disease_snp)

    # 4. Symmetry Check for W (Wij should equal Wji)
    is_symmetric = np.allclose(W, W.T, atol=1e-8)
    if not is_symmetric:
        max_diff = np.abs(W - W.T).max()
        print(f"ALERT: W matrix is NOT symmetric. Max absolute difference: {max_diff}")
    else:
        print("W matrix is confirmed symmetric.")

    # 5. Iterative Evaluation for 427 Rows
    results = []

    for i in range(num_diseases):
        y_true, y_scores, y_preds = [], [], []
        target_phe = disease_list[i]

        # Row i contains similarities of disease i with all others j
        for j in range(num_diseases):
            if i == j: continue  # Skip self-pairing

            other_phe = disease_list[j]
            label = label_map.get((target_phe, other_phe))

            if label is not None:
                wij = W[i, j]
                y_true.append(label)
                y_scores.append(wij)
                # Prediction rule: Wij > 0 -> 1, else 0
                y_preds.append(1 if wij > 0 else 0)

        # Calculate metrics for row i if it contains both positive and negative samples
        if len(y_true) > 0 and len(set(y_true)) > 1:
            results.append({
                'Accuracy': accuracy_score(y_true, y_preds),
                'Precision': precision_score(y_true, y_preds, zero_division=0),
                'Recall': recall_score(y_true, y_preds, zero_division=0),
                'F1': f1_score(y_true, y_preds, zero_division=0),
                'ROC AUC': roc_auc_score(y_true, y_scores),
                'AUC PR': average_precision_score(y_true, y_scores)
            })

    # 6. Aggregate results (Mean +- SD)
    res_df = pd.DataFrame(results)
    summary = pd.DataFrame({
        'Metric': res_df.columns,
        'Mean': res_df.mean().values,
        'SD': res_df.std().values
    })
    summary['Mean +- SD'] = summary.apply(lambda x: f"{x['Mean']:.4f} +- {x['SD']:.4f}", axis=1)

    print("\nSummary Results (W-Matrix Evaluation):")
    print(summary[['Metric', 'Mean +- SD']].to_string(index=False))

    # Save files
    summary.to_csv("/content/drive/MyDrive/w_matrix_performance_summaryRR2.csv", index=False)
    return summary

if __name__ == "__main__":
    run_w_matrix_performance("/content/drive/MyDrive/DiseaseInteractiveNetwork/data_bioinfo/disease_snp.npz",
                             "/content/drive/MyDrive/disease_pairs_with_labelsRR2.csv", "/content/drive/MyDrive/DiseaseInteractiveNetwork/data_bioinfo/phecode_list.csv")