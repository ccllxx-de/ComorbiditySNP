# Load necessary libraries
library(tidyverse)
library(stats)


# --- 1. Define File Names ---
top_negative_file <- "negative_snp_effects_data.csv"
random_negative_file <- "434_random_non_overlapping_snps_matrix.csv"

output_stats_file <- "negative_snp_comparison_summary.csv"


# --- 2. Helper Function ---

# Function to calculate the mean of non-zero negative values for a single SNP column.
# Since negative values are already zeroed out in the input files, we just exclude 0s.
calculate_non_zero_mean <- function(x) {
  # Select values that are smaller than zero (i.e., the original negative effects)
  x_negative_non_zero <- x[x < 0 & !is.na(x)]
  
  # Return NA if no negative non-zero values are present (although file construction should prevent this)
  if (length(x_negative_non_zero) == 0) {
    return(NA)
  }
  return(mean(x_negative_non_zero))
}

# ==============================================================================
# 3. Calculate Mean negative Effects for Each Dataset
# ==============================================================================

# A. Top negative SNPs
df_top <- read.csv(top_negative_file, check.names = FALSE)
top_snp_cols <- colnames(df_top)[-1]
# top_snp_cols

df_top_snp_data <- df_top %>% select(all_of(top_snp_cols))
# df_top_snp_data

df_top_negative <- df_top_snp_data %>%
  mutate(across(everything(), ~ ifelse(. < 0, ., 0)))
# df_top_negative

top_snps_means <- apply(df_top_negative %>% select(all_of(top_snp_cols)), 2, calculate_non_zero_mean)
# top_snps_means

top_negative_means_df <- tibble(
  Group = "Top negative SNPs",
  Mean_negative_Effect = top_snps_means[!is.na(top_snps_means)]
)
# top_negative_means_df
top_means_vector <- top_snps_means[!is.na(top_snps_means)]
# top_means_vector 

# B. Random negative SNPs
df_random <- read.csv(random_negative_file, check.names = FALSE)
# ncol(df_random)

random_snp_cols <- colnames(df_random)[-1] 
# random_snp_cols

# Separate the non-SNP column (Phenotype) from the SNP data
df_random_snp_data <- df_random %>% select(all_of(random_snp_cols))
ncol(df_random_snp_data)
nrow(df_random_snp_data)

all(colSums(df_random_snp_data < 0) > 0)

df_random_negative <- df_random_snp_data %>%
  mutate(across(everything(), ~ ifelse(. < 0, ., 0)))
# df_random_negative


random_snps_means <- apply(df_random_negative %>% select(all_of(random_snp_cols)), 2, calculate_non_zero_mean)
# random_snps_means

random_negative_means_df <- tibble(
  Group = "Random negative SNPs",
  Mean_negative_Effect = random_snps_means[!is.na(random_snps_means)]
)
# random_negative_means_df
random_means_vector <- random_snps_means[!is.na(random_snps_means)]
# random_means_vector


# Combine data for plotting
comparison_data <- bind_rows(top_negative_means_df, random_negative_means_df)
# comparison_data

write.csv(comparison_data, "topandrandom_neg_means_df.csv", row.names = FALSE)

# ==============================================================================
# 4. Descriptive Statistics Summary
# ==============================================================================

overall_stats <- comparison_data %>%
  group_by(Group) %>%
  summarise(
    N_SNPs_Analyzed = n(),
    Mean_of_Nonzero_Means = mean(Mean_negative_Effect, na.rm = TRUE),
    SD_of_Nonzero_Means = sd(Mean_negative_Effect, na.rm = TRUE)
  )

cat("\n--- Descriptive Statistics Summary ---\n")
print(overall_stats)
#already saved
write.csv(overall_stats, output_stats_file, row.names = FALSE)

# ==============================================================================
# 5. Statistical Comparison (K-S Test)
# ==============================================================================

ks_test_result <- ks.test(top_means_vector, random_means_vector)

cat("\n--- Kolmogorov-Smirnov Test Results ---\n")
print(ks_test_result)




