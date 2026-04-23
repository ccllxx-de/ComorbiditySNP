import pandas as pd
import numpy as np
import os
from itertools import combinations

# Load disease_snp matrix (427 diseases x 39382 SNPs)
#disease_snp = np.random.rand(427, 39382)  # Replace with actual data loading
#disease_snp

# Load disease list (427 diseases)
#disease_list = pd.read_csv("/content/drive/MyDrive/DiseaseInteractiveNetwork/data_bioinfo/phecode_list.csv", header=None)[0].tolist()
#disease_list

# Directory containing ground truth files
ground_truth_dir = "/content/drive/MyDrive/DiseaseInteractiveNetwork/data_bioinfo/EHR_comorbidity_UKBB/"

# Initialize a list to store results
results = []

# Iterate through all pairs of diseases
for disease1, disease2 in combinations(disease_list, 2):
    # Construct the ground truth file path for disease1
    ground_truth_file = os.path.join(ground_truth_dir, f"{disease1}.pkl")

    '''
    # Check if the ground truth file exists
    if not os.path.exists(ground_truth_file):
        print(f"Ground truth file not found for disease: {disease1}")
        continue
    '''

    # Load the ground truth file
    ground_truth = pd.read_pickle(ground_truth_file)

    # Find the row corresponding to disease2
    row = ground_truth[ground_truth["PheCode"] == disease2]

    # If disease2 is not found in the ground truth file, skip this pair
    if row.empty:
        print(f"Disease {disease2} not found in ground truth file for {disease1}")
        continue

    # Get the RR value for disease2
    rr_value = row["RR"].values[0]

    # Determine the label based on RR value
    label = 1 if rr_value > 4 else 0

    # Create the disease pair ID
    disease_pair_id = f"{disease1}-{disease2}"

    # Append the result to the list
    results.append({"disease pair ID": disease_pair_id, "label": label})

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)


# Save the results to a CSV file
output_file = "/content/drive/MyDrive/disease_pairs_with_labelsRR4.csv"
results_df.to_csv(output_file, index=False)
print(f"Results saved to {output_file}")