# Load necessary libraries
library(tidyverse)

input_file <- "disease_snp_matrixFull.csv"
output_file <- "952_randomly_selected_snps.csv"
num_snps_to_select <- 952

# Set the seed for reproducibility
set.seed(42)

# --- 2. Load and Prepare Data ---

df_original <- read.csv(input_file, row.names = 1, header = TRUE)
# csv file has no phenotype column
#row.names = 1 the first column should be used as row names

all_cols <- colnames(df_original)

# Identify the columns representing the total pool of SNPs 
snp_cols_pool <- all_cols

# --- 3. Randomly Sample SNP Columns ---

# Check if the desired number of SNPs can be sampled
if (length(snp_cols_pool) < num_snps_to_select) {
  stop("Error: The data file contains fewer SNPs than the requested 952.")
}

# Randomly sample 952 SNP IDs without replacement
selected_snps <- sample(snp_cols_pool, num_snps_to_select, replace = FALSE)

# --- 4. Create and Save the New Dataset ---

# Select the randomly chosen SNP columns
df_selected <- df_original %>%
  select(all_of(selected_snps))
# df_selected

# Save the resulting data frame to a CSV file
write.csv(df_selected, output_file, row.names = TRUE)

cat(paste("Successfully selected", ncol(df_selected), "SNPs.\n"))
cat(paste("The randomly selected dataset is saved as:", output_file, "\n"))


# --- 1. Define File Names ---
top_snps_file <- "952_top_snps_original_values.csv"
random_snps_file <- "952_randomly_selected_snps.csv"
output_file <- "overlapping_snp_ids.txt"

# --- 2. Extract SNP IDs from Top SNPs Data ---
# Load the data, using check.names=FALSE to preserve rsID format
df_top <- read.csv(top_snps_file, check.names = FALSE)

# Extract column names, assuming the first column is the non-SNP column ("Phenotype")
top_snp_ids <- colnames(df_top)[-1]
count_top <- length(top_snp_ids)
count_top

# --- 3. Extract SNP IDs from Random SNPs Data ---
df_random <- read.csv(random_snps_file, check.names = FALSE)

# Extract column names, assuming the first column is the non-SNP column (Phenotype/Index)
random_snp_ids <- colnames(df_random)[-1]
count_random <- length(random_snp_ids)
count_random

# --- 4. Find the Overlap (Intersection) ---
# The intersect() function finds elements common to both vectors
overlap_snps <- intersect(top_snp_ids, random_snp_ids)
count_overlap <- length(overlap_snps)
count_overlap

# --- 5. Print Results and Save Overlap List ---

cat("--- SNP Overlap Comparison ---\n")
cat(paste("Top SNPs Count:", count_top, "\n"))
cat(paste("Random SNPs Count:", count_random, "\n"))
cat(paste("Overlapping (Common) SNPs Count:", count_overlap, "\n"))

if (count_overlap > 0) {
  # Save the list of overlapping SNPs
  writeLines(overlap_snps, output_file)
  cat(paste("\nSuccessfully found", count_overlap, "overlapping SNPs.\n"))
  cat(paste("The list of overlapping SNPs is saved to:", output_file, "\n"))
} else {
  cat("\nNo overlapping SNPs were found between the two datasets.\n")
}

