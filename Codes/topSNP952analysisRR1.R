
# Load necessary libraries
library(tidyverse)
library(patchwork) # For combining plots (optional, but good practice)

library(readr)

setwd("/Users/lixianc/Desktop/LC/UD/Research Liao/")

# --- 1. Define File Names ---
input_file <- "rf_top_1000snp_loadings_RR1SNPvectoralldiseasepairPCA100SDVar.csv"
output_file <- "top_snp_rsid_list.txt"

# --- 2. Load the Data ---
# Note: Using read.csv for standard CSV files
df_loadings <- read.csv(input_file, header = TRUE)
#contain a SNP name that's not valid

# --- 3. Extract the SNP_RSID Column ---
snp_rsids <- df_loadings$SNP_RSID

# --- 4. Save to a Text File ---
# writeLines is used to write a character vector, one element per line
#already saved
writeLines(as.character(snp_rsids), output_file)

length(snp_rsids)
cat(paste("Successfully extracted", length(snp_rsids), "SNP_RSIDs and saved to", output_file, "\n"))

# --- 1. Load SNP List and Main Data ---

snps <- read.csv("rf_top_1000snpIDs_loadings_RR1SNPvectoralldiseasepairPCA100SDVar.csv",
                 header = TRUE, stringsAsFactors = FALSE)
#contain a SNP name that's not valid

# The first column is named "SNP_RSID" – extract it
snp_ids <- snps$SNP_RSID

top_snps_list <- readLines("top_snp_rsid_list.txt")

# Load the main disease SNP matrix
df_original <- read.csv("disease_snp_matrixFull.csv",header = TRUE)
colnames(df_original)[1] <- "Phenotype"

# Select and prepare data
df_target <- df_original %>%
  select(Phenotype, any_of(top_snps_list))
# df_target
ncol(df_target)

# already saved, 
write.csv(df_target, "952_top_snps_original_values.csv", row.names = FALSE)

target_snp_cols <- names(df_target)[2:ncol(df_target)]
# target_snp_cols
num_snps_analyzed <- length(target_snp_cols)
# num_snps_analyzed
cat(paste("Number of target SNPs successfully found and analyzed:", num_snps_analyzed, "\n"))
#952

# Ensure SNP data is numeric
target_snp_data <- df_target %>% select(all_of(target_snp_cols))
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
# orig_summary_list
orig_summary_df <- do.call(rbind, orig_summary_list) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SNP")
# orig_summary_df

# Calculate overall statistics
mean_of_means_orig <- mean(orig_summary_df$Mean, na.rm = TRUE)
# mean_of_means_orig
sd_of_means_orig <- sd(orig_summary_df$Mean, na.rm = TRUE)
# sd_of_means_orig 

# Save summary, already saved.  952
write.csv(orig_summary_df, "952_top_snps_original_values_summary.csv", row.names = FALSE)


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
write.csv(abs_summary_df, "952_top_snps_absolute_values_summary.csv", row.names = FALSE)


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
write.csv(overall_stats_df, "952_top_snps_overall_statistics.csv", row.names = FALSE)


# ==============================================================================
# 6. Phenotype (PheCode) Information
# ==============================================================================

# Filter the original data for rows where at least one target SNP has a non-zero value
df_target_with_pheno <- df_target %>%
  # Temporarily replace NA with 0 for filtering purposes
  mutate(across(all_of(target_snp_cols), ~replace_na(., 0))) %>%
  rowwise() %>%
  # Filter rows where the absolute sum of the target SNP values is greater than 0
  filter(sum(abs(c_across(all_of(target_snp_cols))), na.rm = TRUE) != 0) %>%
  ungroup()

#df_target_with_pheno
#952 SNPs

# Extract unique associated phenotypes
associated_phenotypes <- unique(df_target_with_pheno$Phenotype)
# associated_phenotypes
df_phenotypes <- data.frame(Phenotype = associated_phenotypes)
# df_phenotypes

# Save summary, already saved
write.csv(df_phenotypes, "952_top_snps_associated_phenotypes.csv", row.names = FALSE)
cat(paste("Extracted", nrow(df_phenotypes), "unique phenotypes associated with these SNPs.\n"))


# ==============================================================================
# 7. Plotting Absolute Value Distribution
# ==============================================================================

# Convert absolute SNP data to long format and filter out zeros/NAs
df_abs_long2 <- target_snp_data_abs %>%
  # Convert to long format
  pivot_longer(cols = everything(), names_to = "SNP", values_to = "Abs_Effect_Magnitude") %>%
  # Filter out zeros and NAs for the distribution plot
  filter(Abs_Effect_Magnitude != 0 & !is.na(Abs_Effect_Magnitude))

#already saved
write.csv(df_abs_long2, "df_abs_long2_952_top_snps_group_abs.csv", row.names = TRUE)

# run the other R script that generates df_abs_long before moving on 

library(ggplot2)
library(dplyr)

# Combine data
df_abs_long$Group  <- "Random SNPs"
df_abs_long2$Group <- "Target/Top SNPs"
length(df_abs_long$SNP)
length(df_abs_long2$SNP)
length(unique(df_abs_long$SNP))
length(unique(df_abs_long2$SNP))
combined <- rbind(df_abs_long, df_abs_long2)
combined

write.csv(combined, "952_topvsrandom_snps_groups_abs.csv", row.names = TRUE)