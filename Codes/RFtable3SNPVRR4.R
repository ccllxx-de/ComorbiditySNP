library(tidyverse)
library(readr)

# Load the dataset
df <- read.csv("Concat_CV10Fold_Metrics_SNP_PCA_Detailed_RF_RR4.csv")

# Filter the data for the Random Forest model
rf_df <- subset(df, Model == "Random Forest")

# List the metrics to be calculated
metrics <- c("Accuracy", "Precision", "Recall", "F1", "ROC_AUC", "PR_AUC")

# Calculate Mean and SD for each metric
means <- sapply(rf_df[metrics], mean)
sds <- sapply(rf_df[metrics], sd)

# Combine into a clean table format
results <- data.frame(
  Metric = metrics,
  Mean = means,
  Std_Dev = sds
)

print("Random Forest Performance Metrics:")
print(results)
