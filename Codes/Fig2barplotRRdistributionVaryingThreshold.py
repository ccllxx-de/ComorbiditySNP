import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Loop through all diseases (idx_ from 0 to 426)
for idx_ in range(len(disease_list)):
    # Predicting direct/inverse comorbidity
    f = np.zeros(len(W))
    y = np.zeros(len(W))
    y[idx_] = 1  # Set seed label for the current disease
    f = inversed_L @ y  # Predict outcomes


    target_disease = disease_list[idx_]
    Labeled_id = np.where(disease_info == target_disease)[0]

    if Labeled_id.size == 0 or Labeled_id[0] != idx_:
        print(f"Index {idx_} and phenotypes do not match for disease {target_disease}")
        continue  # Skip to the next iteration if there's a mismatch

    # Load ground truth data
    RR_phi_location = f"/content/drive/MyDrive/DiseaseInteractiveNetwork/data_bioinfo/EHR_comorbidity_UKBB/{disease_list[Labeled_id][0]}.pkl"
    try:
        Ground_truth = pd.read_pickle(RR_phi_location)
    except FileNotFoundError:
        print(f"Ground truth file not found for disease {target_disease} at {RR_phi_location}")
        continue  # Skip to the next iteration if the file is missing

    # Create summary table
    f_score = pd.DataFrame(f, columns=['f_score'])
    Summary_table = pd.concat([
        Ground_truth["PheCode"],
        disease_info["description"],
        disease_info["category"],
        f_score,
        Ground_truth["RR"],
        Ground_truth["phi"],
        Ground_truth["pval"]
    ], axis=1)


    # Save results
    #out_filename = f"/content/drive/MyDrive/alldiseases/{disease_info['code2'][idx_]}.csv"
    out_filename = f"/content/drive/MyDrive/alldiseasesPairedOther/{disease_info['code2'][idx_]}.csv"
    Summary_table.to_csv(out_filename, index=False)
    print(f"Results saved for {target_disease} at {out_filename}")

print("Processing complete for all diseases.")


# Directory containing the 427 CSV files
input_directory = "/content/drive/MyDrive/alldiseasesPairedOther"  # Replace with your folder path

# Initialize an empty list to store all combined rows
all_rows = []

# Iterate through all CSV files in the directory
for filename in os.listdir(input_directory):
    if filename.endswith(".csv"):
        # Extract the last 4 digits of the filename (before .csv)
        file_last_4_digits = filename[-8:-4]  # Assumes filenames are in the format "XXXX.csv"

        # Load the CSV file
        file_path = os.path.join(input_directory, filename)
        df = pd.read_csv(file_path)

        # Remove the row where the last 4 digits of PheCode match the filename's last 4 digits
        df = df[~df["PheCode"].astype(str).str[-4:].eq(file_last_4_digits)]

        # Append the remaining rows to the list
        all_rows.append(df)

# Combine all remaining rows into a single DataFrame
combined_df = pd.concat(all_rows, ignore_index=True)

# Extract RR values
rr_values = combined_df["RR"].values

# Define thresholds
thresholds = [1, 2, 3, 4, 5, 6,7,8,9,10, 15]

# Calculate the percentage of RR values greater than each threshold
percentages = {}
for threshold in thresholds:
    percentage = (np.sum(rr_values > threshold) / len(rr_values)) * 100
    percentages[threshold] = percentage
    print(f"Percentage of RR values greater than {threshold}: {percentage:.2f}%")

# Create a bar plot for the percentages
plt.figure(figsize=(10, 6))
plt.bar(percentages.keys(), percentages.values(), color='blue')
plt.xlabel("Threshold (RR Value)")
plt.ylabel("Percentage of RR Values > Threshold")
plt.title("Percentage of RR Values Greater Than Different Thresholds")
plt.xticks(list(percentages.keys()))
plt.grid(axis='y')
plt.show()