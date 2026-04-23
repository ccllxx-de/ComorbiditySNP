library(tidyverse)


random_snps_file <- "952_random_snps_non_overlapping.csv"

df_random <- read.csv(random_snps_file, check.names = FALSE)
df_random2 <- df_random
colnames(df_random2)[1] <- "Phenotype"
random_snp_cols <- colnames(df_random)[-1]
random_snp_cols

random_snp_cols2 <- names(df_random)[2:ncol(df_random)]
random_snp_cols2
num_snps_analyzed <- length(random_snp_cols)
num_snps_analyzed
cat(paste("Number of target SNPs successfully found and analyzed:", num_snps_analyzed, "\n"))
#952


# Ensure SNP data is numeric
target_snp_data <- df_random %>% select(all_of(random_snp_cols))
target_snp_data <- data.frame(lapply(target_snp_data, function(x) as.numeric(as.character(x))))
# target_snp_data 

# ==============================================================================
# 2. Helper Function for Statistics (Excluding Zeros)
# ==============================================================================

# Custom function to calculate Mean and SD, excluding 0 and NA values
calculate_non_zero_stats <- function(x) {
  x_non_zero <- x[x != 0 & !is.na(x)]
  
  # Return Mean and SD, returning NA if not enough non-zero data
  stats <- c(
    Mean = if (length(x_non_zero) > 0) mean(x_non_zero) else NA,
    SD = if (length(x_non_zero) > 1) sd(x_non_zero) else NA
  )
  return(stats)
}

calculate_non_zero_summary <- function(x) {
  x_non_zero <- x[x != 0 & !is.na(x)]
  
  stats <- c(
    Min. = min(x_non_zero),
    #`1st Qu.` = quantile(x_non_zero, 0.25),
    Median = median(x_non_zero),
    Mean = mean(x_non_zero),
    #`3rd Qu.` = quantile(x_non_zero, 0.75),
    Max. = max(x_non_zero),
    SD = sd(x_non_zero)
  )
  return(stats)
}

# ==============================================================================
# 3. Analysis on ORIGINAL Values (Excluding Zeros)
# ==============================================================================

# Apply the custom summary function
orig_summary_list <- lapply(target_snp_data, calculate_non_zero_summary)
orig_summary_list
orig_summary_df <- do.call(rbind, orig_summary_list) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SNP")
orig_summary_df

# Calculate overall statistics
mean_of_means_orig <- mean(orig_summary_df$Mean, na.rm = TRUE)
mean_of_means_orig
sd_of_means_orig <- sd(orig_summary_df$Mean, na.rm = TRUE)
sd_of_means_orig 

# Save summary, already saved.  952
write.csv(orig_summary_df, "952_random_snps_original_values_summary.csv", row.names = FALSE)


# ==============================================================================
# 4. Analysis on ABSOLUTE Values (Excluding Zeros)
# ==============================================================================

# Convert all non-zero values to their absolute magnitude
target_snp_data_abs <- abs(target_snp_data)
target_snp_data_abs[target_snp_data == 0] <- 0 # Keep zeros as they were

# Apply the custom summary function to absolute values
abs_summary_list <- lapply(target_snp_data_abs, calculate_non_zero_summary)
abs_summary_df <- do.call(rbind, abs_summary_list) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SNP")
# abs_summary_df

# Calculate overall statistics
mean_of_means_abs <- mean(abs_summary_df$Mean, na.rm = TRUE)
# mean_of_means_abs
sd_of_means_abs <- sd(abs_summary_df$Mean, na.rm = TRUE)
# sd_of_means_abs

# Save summary, already saved
write.csv(abs_summary_df, "952_random_snps_absolute_values_summary.csv", row.names = FALSE)


# ==============================================================================
# 5. Summary of Overall Statistics
# ==============================================================================

overall_stats_df <- data.frame(
  Metric = c("Mean of Means (Original)", "SD of Means (Original)",
             "Mean of Means (Absolute)", "SD of Means (Absolute)"),
  Value = c(mean_of_means_orig, sd_of_means_orig,
            mean_of_means_abs, sd_of_means_abs),
  Note = "Calculated from non-zero SNP effect magnitudes only"
)
# overall_stats_df

# Save summary, already saved
write.csv(overall_stats_df, "952_random_snps_overall_statistics.csv", row.names = FALSE)


# ==============================================================================
# 6. Phenotype (PheCode) Information
# ==============================================================================

# Filter the original data for rows where at least one target SNP has a non-zero value
df_random_with_pheno <- df_random2 %>%
  # Temporarily replace NA with 0 for filtering purposes
  mutate(across(all_of(random_snp_cols), ~replace_na(., 0))) %>%
  rowwise() %>%
  # Filter rows where the absolute sum of the target SNP values is greater than 0
  filter(sum(abs(c_across(all_of(random_snp_cols))), na.rm = TRUE) != 0) %>%
  ungroup()

# df_random_with_pheno
# 952 SNPs

# Extract unique associated phenotypes
associated_phenotypes <- unique(df_random_with_pheno$Phenotype)
# associated_phenotypes
df_phenotypes <- data.frame(Phenotype = associated_phenotypes)
# df_phenotypes

# Save summary, already saved
write.csv(df_phenotypes, "952_random_snps_associated_phenotypes.csv", row.names = FALSE)
cat(paste("Extracted", nrow(df_phenotypes), "unique phenotypes associated with these SNPs.\n"))
#360

# ==============================================================================
# 7. Plotting Absolute Value Distribution
# ==============================================================================

# Convert absolute SNP data to long format and filter out zeros/NAs
df_abs_long <- target_snp_data_abs %>%
  # Convert to long format
  pivot_longer(cols = everything(), names_to = "SNP", values_to = "Abs_Effect_Magnitude") %>%
  # Filter out zeros and NAs for the distribution plot
  filter(Abs_Effect_Magnitude != 0 & !is.na(Abs_Effect_Magnitude))

write.csv(df_abs_long, "df_abs_long_952_random_snps_group_abs.csv", row.names = TRUE)

