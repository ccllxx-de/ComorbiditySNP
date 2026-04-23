# Load necessary libraries
library(tidyverse)

# --- 1. Define File Names ---
input_file <- "952_top_snps_original_values.csv"
output_positive_data <- "positive_snp_effects_data.csv"
output_negative_data <- "negative_snp_effects_data.csv"
output_summary <- "positive_negative_summary_statistics.csv"

# --- 2. Load and Prepare Data ---
df_original <- read.csv(input_file, check.names = FALSE)
colnames(df_original)[1] <- "Phenotype"
snp_cols <- colnames(df_original)[-1]
num_original_snps <- length(snp_cols)

# Separate the non-SNP column (Phenotype) from the SNP data
df_snp_data <- df_original %>% select(all_of(snp_cols))


# --- 3. Process Positive SNP Effects ---

# Create a dataset where negative values are set to 0
# The logic: if value is > 0, keep it; otherwise (0 or negative), set to 0.
df_positive <- df_snp_data %>%
  mutate(across(everything(), ~ ifelse(. > 0, ., 0)))

# Calculate the mean of raw positive values for each SNP (excluding 0s)
# Note: Since we set negative values to 0, we now exclude 0s to only average the positive effects.
positive_means <- df_positive %>%
  summarise(across(everything(), ~ mean(.[. != 0], na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "SNP_ID", values_to = "Mean_Positive_Effect")

# Filter out SNPs where all values were 0 or negative (resulting in NA mean)
positive_means_nonzero <- positive_means %>%
  filter(!is.na(Mean_Positive_Effect))

# positive_means_nonzero


# Report number of SNPs with at least one positive effect
num_positive_snps <- nrow(positive_means_nonzero)
cat(paste("Number of SNPs with at least one positive effect:", num_positive_snps, "out of", num_original_snps, "\n"))


# Calculate overall summary statistics for the non-zero means
positive_summary <- positive_means_nonzero %>%
  summarise(
    Set = "Positive Effects",
    N_SNPs_with_Effect = n(),
    Mean_of_Nonzero_Means = mean(Mean_Positive_Effect),
    SD_of_Nonzero_Means = sd(Mean_Positive_Effect)
  )
# positive_summary

# Save the positive effects dataset (with negative effects zeroed out)
# Combine with Phenotype column for completeness
# df_positive_final <- bind_cols(df_original %>% select(Phenotype), df_positive)
# df_positive_final
#write.csv(df_positive_final, output_positive_data, row.names = FALSE)


# --- 4. Process Negative SNP Effects ---

# Create a dataset where positive values are set to 0
# The logic: if value is < 0, keep it; otherwise (0 or positive), set to 0.
df_negative <- df_snp_data %>%
  mutate(across(everything(), ~ ifelse(. < 0, ., 0)))

# Calculate the mean of raw negative values for each SNP (excluding 0s)
# Note: Since we set positive values to 0, we now exclude 0s to only average the negative effects.
negative_means <- df_negative %>%
  summarise(across(everything(), ~ mean(.[. != 0], na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "SNP_ID", values_to = "Mean_Negative_Effect")

# Filter out SNPs where all values were 0 or positive (resulting in NA mean)
negative_means_nonzero <- negative_means %>%
  filter(!is.na(Mean_Negative_Effect))

# Report number of SNPs with at least one negative effect
num_negative_snps <- nrow(negative_means_nonzero)
cat(paste("Number of SNPs with at least one negative effect:", num_negative_snps, "out of", num_original_snps, "\n"))

# Calculate overall summary statistics for the non-zero means
negative_summary <- negative_means_nonzero %>%
  summarise(
    Set = "Negative Effects",
    N_SNPs_with_Effect = n(),
    Mean_of_Nonzero_Means = mean(Mean_Negative_Effect),
    SD_of_Nonzero_Means = sd(Mean_Negative_Effect)
  )

# negative_summary

# Save the negative effects dataset (with positive effects zeroed out)
# df_negative_final <- bind_cols(df_original %>% select(Phenotype), df_negative)
#write.csv(df_negative_final, output_negative_data, row.names = FALSE)


# --- 5. Combine and Save Final Summary Statistics ---
final_summary_stats <- bind_rows(positive_summary, negative_summary)
# final_summary_stats
write.csv(final_summary_stats, output_summary, row.names = FALSE)

cat("\n--- Final Summary Statistics ---\n")
print(final_summary_stats)
# cat(paste("\nData for positive effects saved as:", output_positive_data, "\n"))
# cat(paste("Data for negative effects saved as:", output_negative_data, "\n"))
cat(paste("Overall summary statistics saved as:", output_summary, "\n"))


output_overlap_snps <- "overlapping_sign_snps_data.csv"
output_overlap_list <- "overlapping_sign_snps_list.txt"


# List of SNPs that contain at least one positive effect
positive_snp_ids <- positive_means %>%
  filter(!is.na(Mean_Positive_Effect)) %>%
  pull(SNP_ID)

# List of SNPs that contain at least one negative effect
negative_snp_ids <- negative_means %>%
  filter(!is.na(Mean_Negative_Effect)) %>%
  pull(SNP_ID)


# --- 5. Find the Intersection (Overlap) ---

overlap_snp_ids <- intersect(positive_snp_ids, negative_snp_ids)
overlap_snp_ids 

num_overlap <- length(overlap_snp_ids)
num_overlap


cat("--- Overlap Analysis ---\n")
cat(paste("Total SNPs analyzed:", num_original_snps, "\n"))
cat(paste("SNPs with only positive effects or both:", length(positive_snp_ids), "\n"))
cat(paste("SNPs with only negative effects or both:", length(negative_snp_ids), "\n"))
cat(paste("SNPs found in BOTH positive and negative sets (Overlap):", num_overlap, "\n"))

# Save the list of overlapping SNP IDs
writeLines(overlap_snp_ids, output_overlap_list)
cat(paste("The list of", num_overlap, "overlapping SNPs is saved to:", output_overlap_list, "\n"))


# ==============================================================================
# 4. Generate and Save df_positive_final (Filtered Columns, Negative -> 0)
# ==============================================================================

# Select only the Phenotype column and the SNPs identified as having positive values
df_positive_final <- df_original %>%
  select(Phenotype, all_of(positive_snp_ids))

# Apply the transformation: set all negative values in the SNP columns to 0
# The 'Phenotype' column is excluded from this transformation
# df_positive_final <- df_positive_final %>%
#   mutate(across(all_of(positive_snp_ids), ~ ifelse(. < 0, 0, .)))

nrow(df_positive_final)

ncol(df_positive_final)


write.csv(df_positive_final, output_positive_data, row.names = FALSE)
cat(paste("\nFiltered positive dataset saved as:", output_positive_data, "\n"))


# ==============================================================================
# 5. Generate and Save df_negative_final (Filtered Columns, Positive -> 0)
# ==============================================================================

# Select only the Phenotype column and the SNPs identified as having negative values
df_negative_final <- df_original %>%
  select(Phenotype, all_of(negative_snp_ids))

# Apply the transformation: set all positive values in the SNP columns to 0
# The 'Phenotype' column is excluded from this transformation
# df_negative_final <- df_negative_final %>%
  # mutate(across(all_of(negative_snp_ids), ~ ifelse(. > 0, 0, .)))

nrow(df_negative_final)

ncol(df_negative_final)


write.csv(df_negative_final, output_negative_data, row.names = FALSE)
cat(paste("Filtered negative dataset saved as:", output_negative_data, "\n"))
