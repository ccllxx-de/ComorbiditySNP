# Load necessary libraries
library(tidyverse)
library(dplyr)
library(tidyr)

data_file <- "disease_snp_matrixFull.csv"
exclusion_file <- "952_top_snps_original_values.csv"
#output_file <- "952_random_snps_non_overlapping.csv"
output_file_pos <- "520_random_snps_positive_non_overlapping.csv"
output_file_neg <- "434_random_snps_negative_non_overlapping.csv"
num_snps_to_select <- 952
num_pos_snps_to_select <- 520
num_neg_snps_to_select <- 434

# Set the seed for reproducible random sampling
set.seed(42)

# --- 2. Load and Prepare Data ---

# Load the full SNP matrix (Total Pool)
df_original2 <- read.csv(data_file, check.names = FALSE)

colnames(df_original2)[1] <- "Phenotype"
total_snp_pool <- colnames(df_original2)[-1]
cat(paste("Total SNPs in original matrix:", length(total_snp_pool), "\n"))

# Load the 952 Top SNPs (Exclusion Set)
df_exclusion <- read.csv(exclusion_file, check.names = FALSE)
exclusion_snps <- colnames(df_exclusion)[-1]
cat(paste("SNPs to be excluded:", length(exclusion_snps), "\n"))


# --- 3. Create the Non-Overlapping Pool ---

# Use setdiff() to subtract the exclusion set from the total pool
non_overlapping_pool <- setdiff(total_snp_pool, exclusion_snps)
cat(paste("Remaining SNPs available for random selection:", length(non_overlapping_pool), "\n"))

selected_snps <- sample(non_overlapping_pool, num_snps_to_select, replace = FALSE)
selected_snps

write.csv(non_overlapping_pool, "non_overlapping_pool.csv", row.names = TRUE)

typeof(non_overlapping_pool)

df_pool <- read.csv("non_overlapping_pool.csv", check.names = FALSE)

snp_names <- df_pool$x

if ("x" %in% names(df_pool)) {
  pool_ids <- df_pool %>% pull(x)
}

pool_ids <- as.character(pool_ids)
cat(paste("Total SNPs in non-overlapping pool:", length(pool_ids), "\n"))

df_matrix_full <- read.csv(data_file, check.names = FALSE)
matrix_cols <- colnames(df_matrix_full)
selected_cols <- intersect(pool_ids, matrix_cols)
selected_cols

df_working_pool <- df_matrix_full %>% 
  select(all_of(selected_cols))

cat(paste("Number of SNP columns successfully extracted from the matrix:", ncol(df_working_pool), "\n"))

# ==============================================================================
# 3. Identify Potential Candidates
# ==============================================================================

# Helper functions to check for positive/negative values
has_positive <- function(x) { any(x > 0, na.rm = TRUE) }
has_negative <- function(x) { any(x < 0, na.rm = TRUE) }

# Check every column for positive/negative values
snp_checks <- df_working_pool %>%
  summarise(
    across(everything(), .fns = list(
      HasPositive = has_positive,
      HasNegative = has_negative
    ), .names = "{.col}_{.fn}")
  ) %>%
  pivot_longer(everything(), names_to = c("SNP_ID", ".value"), names_sep = "_")

# List of SNPs that have at least one positive value
positive_candidates <- snp_checks %>%
  filter(HasPositive == TRUE) %>%
  pull(SNP_ID)

# List of SNPs that have at least one negative value
negative_candidates <- snp_checks %>%
  filter(HasNegative == TRUE) %>%
  pull(SNP_ID)

cat(paste("SNPs with >= 1 positive effect:", length(positive_candidates), "\n"))
cat(paste("SNPs with >= 1 negative effect:", length(negative_candidates), "\n"))


# ==============================================================================
# 4. Random Sampling
# ==============================================================================
# Set a seed for reproducible random sampling
set.seed(42) 
target_positive_n <- 520
target_negative_n <- 434

# Check if we have enough candidates for sampling
if (length(positive_candidates) < target_positive_n) {
  stop(paste("Insufficient positive candidates (", length(positive_candidates), ") to sample", target_positive_n, "SNPs."))
}
if (length(negative_candidates) < target_negative_n) {
  stop(paste("Insufficient negative candidates (", length(negative_candidates), ") to sample", target_negative_n, "SNPs."))
}


# -------------------------------------------------------------
# NEW: Subset the full matrix to only the non-overlapping SNPs
# -------------------------------------------------------------
non_overlap_matrix <- df_original2[, c("Phenotype", non_overlapping_pool), drop = FALSE]
# Now: rows = diseases, columns = Phenotype + remaining SNPs

# Extract just the numeric SNP columns as a matrix for fast operations
snp_data <- as.matrix(non_overlap_matrix[, non_overlapping_pool])

# Give proper column names (important!)
colnames(snp_data) <- non_overlapping_pool

# -------------------------------------------------------------
# Identify SNPs with at least one positive / negative value
# -------------------------------------------------------------
has_positive <- colSums(snp_data > 0, na.rm = TRUE) > 0
has_negative <- colSums(snp_data < 0, na.rm = TRUE) > 0

snps_with_positive <- colnames(snp_data)[has_positive]
snps_with_negative <- colnames(snp_data)[has_negative]

cat("SNPs with >=1 positive value :", length(snps_with_positive), "\n")
cat("SNPs with >=1 negative value :", length(snps_with_negative), "\n")

# Safety check
if (length(snps_with_positive) < num_pos_snps_to_select)
  stop("Not enough SNPs with positive values!")
if (length(snps_with_negative) < num_neg_snps_to_select)
  stop("Not enough SNPs with negative values!")

# -------------------------------------------------------------
# Random sampling (with replacement = FALSE → no duplicates within group)
# -------------------------------------------------------------
selected_positive <- sample(snps_with_positive, size = num_pos_snps_to_select, replace = FALSE)
selected_negative <- sample(snps_with_negative, size = num_neg_snps_to_select, replace = FALSE)

# -------------------------------------------------------------
# Save results
# -------------------------------------------------------------
# write.csv(data.frame(SNP = selected_positive),
#           output_file_pos, row.names = FALSE, quote = FALSE)
# write.csv(data.frame(SNP = selected_negative),
#           output_file_neg, row.names = FALSE, quote = FALSE)

pos_matrix <- non_overlap_matrix[, c("Phenotype", selected_positive), drop = FALSE]
neg_matrix <- non_overlap_matrix[, c("Phenotype", selected_negative), drop = FALSE]

write.csv(pos_matrix, "520_random_non_overlapping_snps_matrix.csv",
          row.names = FALSE, quote = FALSE)
write.csv(neg_matrix, "434_random_non_overlapping_snps_matrix.csv",
          row.names = FALSE, quote = FALSE)

selected_all_snps <- unique(c(selected_positive, selected_negative))

# -------------------------------------------------------------
# Summary output
# -------------------------------------------------------------
cat("\n=== FINAL SUMMARY ===\n")
cat("Selected 520 SNPs with at least one positive value → saved to:", output_file_pos, "\n")
cat("Selected 434 SNPs with at least one negative value → saved to:", output_file_neg, "\n")
cat("Total unique SNPs selected:", length(selected_all_snps), "\n")
cat("Overlap between positive and negative sets:", 
    length(intersect(selected_positive, selected_negative)), "SNPs\n") 


# --- 5. Create and Save the New Dataset ---

# Select the Phenotype column and the newly chosen random SNP columns
df_selected <- df_original2 %>%
  select(Phenotype, all_of(selected_snps))
df_selected
ncol(df_selected)
dim(df_selected)

# Save the resulting data frame to a CSV file
write.csv(df_selected, output_file, row.names = TRUE)

cat(paste("\nSuccessfully selected", ncol(df_selected), "non-overlapping random SNPs.\n"))
cat(paste("The resulting dataset is saved as:", output_file, "\n"))

# Load necessary libraries
library(tidyverse)
library(stats)

# --- 1. Define File Names ---
top_snps_file <- "952_top_snps_original_values.csv"
random_snps_file <- "952_random_snps_non_overlapping.csv"
output_stats_file <- "952_snp_comparison_summary.csv"
output_data_file <- "952_snp_comparison_means_data.csv"

# --- 2. Helper Function ---

# Function to calculate the mean of non-zero absolute values for a single SNP column
calculate_non_zero_abs_mean <- function(x) {
  x_abs_non_zero <- abs(x[x != 0 & !is.na(x)])
  
  # Return NA if no non-zero values are present
  if (length(x_abs_non_zero) == 0) {
    return(NA)
  }
  return(mean(x_abs_non_zero))
}


# ==============================================================================
# 3. Calculate Mean Absolute Effects for Each Dataset
# ==============================================================================

# A. Top SNPs (File 1)
df_top2 <- read.csv(top_snps_file, check.names = FALSE)
top_snp_cols <- colnames(df_top2)[-1]

top_snps_means <- apply(df_top2 %>% select(all_of(top_snp_cols)), 2, calculate_non_zero_abs_mean)
# 2 for columns
top_snps_means

top_means_df <- tibble(
  Group = "Top 952 SNPs",
  Mean_Abs_Effect = top_snps_means[!is.na(top_snps_means)]
)

top_means_df

# Vector for statistical test
top_means_vector <- top_snps_means[!is.na(top_snps_means)]
top_means_vector

# B. Random SNPs (File 2)
df_random <- read.csv(random_snps_file, check.names = FALSE)
random_snp_cols <- colnames(df_random)[-1]
random_snp_cols

random_snps_means <- apply(df_random %>% select(all_of(random_snp_cols)), 2, calculate_non_zero_abs_mean)
random_snps_means

random_means_df <- tibble(
  Group = "Random 952 SNPs",
  Mean_Abs_Effect = random_snps_means[!is.na(random_snps_means)]
)
random_means_df 

# Vector for statistical test
random_means_vector <- random_snps_means[!is.na(random_snps_means)]
random_means_vector

# Combine data for plotting and saving
comparison_data <- bind_rows(top_means_df, random_means_df)
write.csv(comparison_data, output_data_file, row.names = FALSE)


# ==============================================================================
# 4. Descriptive Statistics Summary
# ==============================================================================

# Calculate overall mean and SD of the mean absolute effects for each group
overall_stats <- comparison_data %>%
  group_by(Group) %>%
  summarise(
    N_SNPs_Analyzed = n(),
    Mean_of_Abs_Means = mean(Mean_Abs_Effect, na.rm = TRUE),
    SD_of_Abs_Means = sd(Mean_Abs_Effect, na.rm = TRUE)
  )

cat("--- Descriptive Statistics Summary ---\n")
print(overall_stats)
write.csv(overall_stats, output_stats_file, row.names = FALSE)


# ==============================================================================
# 5. Statistical Comparison (Kolmogorov-Smirnov Test)
# ==============================================================================

# Compares the entire distribution of the mean absolute effects of the two groups
ks_test_result <- ks.test(top_means_vector, random_means_vector)

cat("\n--- Statistical Comparison (Kolmogorov-Smirnov Test) ---\n")
print(ks_test_result)

cat("\nInterpretation: A small p-value (< 0.05) indicates that the two distributions are statistically different (i.e., the Top SNPs have significantly higher/lower effect magnitudes than the Random SNPs).\n")
