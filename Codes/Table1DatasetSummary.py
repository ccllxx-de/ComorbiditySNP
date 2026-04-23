import pandas as pd
import os

# --- Google Drive Setup ---
# NOTE: This part is REQUIRED if you are running this in Google Colab.
# You must run this and follow the prompts to link your Drive.
# from google.colab import drive
# drive.mount('/content/drive')

# Define the file paths and their corresponding dataset labels located in your Drive
files = {
    "RR1": "/content/drive/MyDrive/disease_pairs_with_labels.csv",
    "RR2": "/content/drive/MyDrive/disease_pairs_with_labelsRR2.csv",
    "RR4": "/content/drive/MyDrive/disease_pairs_with_labelsRR4.csv"
}

results = []

for dataset_name, file_path in files.items():
    # Read the CSV file using the full Drive path
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}. Ensure your drive is mounted and the path is correct.")
        # Skip to the next file if not found
        continue
    except Exception as e:
        print(f"An error occurred while reading {file_path}: {e}")
        continue

    # Calculate value counts for the 'label' column
    counts = df['label'].value_counts()
    total_count = df.shape[0]

    # Get counts, defaulting to 0 if a label is missing
    positive_count = counts.get(1, 0)
    negative_count = counts.get(0, 0)

    # Calculate percentages
    positive_percentage = (positive_count / total_count) * 100
    negative_percentage = (negative_count / total_count) * 100

    # Store results
    results.append({
        'Dataset': dataset_name,
        'Positive_Count': positive_count,
        'Negative_Count': negative_count,
        'Total_Count': total_count,
        'Positive_Percentage': positive_percentage,
        'Negative_Percentage': negative_percentage
    })

# Create the final summary DataFrame
summary_df = pd.DataFrame(results)

# Format percentages for display
summary_df['Positive_Percentage'] = summary_df['Positive_Percentage'].map('{:.2f}%'.format)
summary_df['Negative_Percentage'] = summary_df['Negative_Percentage'].map('{:.2f}%'.format)

# Print the table
print(summary_df)

output_file = "/content/drive/MyDrive/RRsetsTable1.csv"
summary_df.to_csv(output_file, index=False)
print(f"Results saved to {output_file}")