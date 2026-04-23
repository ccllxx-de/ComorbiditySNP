# ComorbiditySNP

Repository for Detecting disease comorbidity based on SNP association on PheWAS scale

Authors: Lixian Chen and Li Liao

## Key Steps:

### Data Acquisition & Pruning: 

High-dimensional disease-SNP association data were obtained from the UK Biobank PheWAS summary statistics. The dataset was pruned to include 427 diseases and 39,382 SNPs.

Download/obtain all files from the ./data_bioinfo/ folder of https://github.com/dokyoonkimlab/DiseaseInteractiveNetwork;

UK Biobank PheWAS summary statistics were obtained from https://www.leelabsg.org/resources;

Save the following key data files for downstream analysis:

./data_bioinfo/disease_snp.npz  # 427 diseases by 39382 SNPs

./data_bioinfo/full_phecode_list.csv # description of 427 diseases with PheCode

./data_bioinfo/snp_list_rsid.csv # list of used SNPs in this analysis

Save the previously generated results (Predicted scores from graph-based SSL, Relative risk & Phi-correlations obtained from EHR data) for reference: 

./data_bioinfo/predicted_results_for_427_diseases.zip 

Save the previously generated Relative risk & Phi-correlations for 427 diseases: 

./data_bioinfo/EHR_comorbidity_UKBB/ **.pkl 

### Programming Specs: 

Language 1: Python 3 (Python 3.12.13, tested with Python 3.8+) 

Required libraries (install if missing):

pip install pandas numpy matplotlib scipy sklearn os

Recommended environment: Google Colab

Mount your Google Drive at the start of the notebook.

You can also run the scripts locally (Jupyter Notebook, etc.) after changing the file paths.

Language 2: R (RStudio Version 1.4.1717)

Required libraries (install if missing): 

tidyverse dplyr tidyr stats readr 

### Ground Truth Labeling: 

Comorbidity labels were established based on Relative Risk (RR) calculated from EHR-based hospital episode statistics. Three datasets were created with varying thresholds: RR1 (RR > 1.0), RR2 (RR > 2.0), and RR4 (RR > 4.0).

**disease_pairs_with_labels.csv, disease_pairs_with_labelsRR2.csv, disease_pairs_with_labelsRR4.csv**

**diseasepairgroundtruthlabelledRR1.py, diseasepairgroundtruthlabelledRR2.py, and diseasepairgroundtruthlabelledRR4.py**

These three Python scripts are designed to create labeled datasets of disease pairs based on comorbidity relative risk (RR) from UK Biobank EHR data. 
They are nearly identical except for the RR threshold (RR > 1, RR > 2, and RR > 4) used to assign the binary label (0 or 1).

### Feature Engineering (Dimension Reduction): 

Due to the extreme high dimensionality of SNP data (39,382 features), Principal Component Analysis (PCA) was used to reduce the data into lower-dimensional vectors (e.g., 50 components) while modulating information loss.

**RR1SNPvectorsPCAVarianceReduction.py**

This Python script performs PCA dimensionality reduction on the full disease-SNP matrix and calculates how much variance is retained when reducing the original 39,382-dimension SNP features down to 50 principal components.

### The Methods/Modeling Approaches:

SNP Vector (SNP-V): Directly used concatenated reduced SNP vectors for a disease pair as input to a classifier.

Disease-Row Vector (DR-V): Used the entire row from a disease-disease similarity matrix (DDN) as a feature vector to capture broader network-level similarities.

Model Training & Evaluation: Supervised learning was conducted using Neural Networks and Random Forest classifiers, evaluated via 10-fold cross-validation.

SNP Traceability: Using the Random Forest classifier in the SNP-V approach, we traced important features back to specific SNPs via eigenvector projection to identify top SNPs for comorbidity prediction.

**RR1SNPvectorsTopFeaturesTopSNPsRandomForest.py**

This Python script implements a Random Forest classification pipeline using PCA-reduced SNP vectors as features for disease-pair prediction with the RR > 1 labeling threshold. It first standardizes the full disease-SNP matrix, reduces the 39,382 SNP features to 50 principal components per disease via PCA, concatenates the 50 PCs from both diseases in each pair to create 100-dimensional feature vectors, trains a Random Forest classifier, runs 10-fold stratified cross-validation plus a held-out test set evaluation, and finally performs feature importance analysis to rank and extract the top 1,000 most relevant SNPs based on ranked eigenvectors by absolute values.

To use the same pipeline for RR2 or RR4 thresholds of ground truth labeling, replace the RR1 labelled file "./disease_pairs_with_labels.csv" with the corresponding CSV label file, i.e. "./disease_pairs_with_labelsRR2.csv" for RR2 labelled file, "./disease_pairs_with_labelsRR4.csv" for RR4 labelled file. 

### Key Tables:

**Table 1 (Datasets)**: Shows the distribution of comorbid vs. non-comorbid pairs. For example, RR1 is highly imbalanced (88% positive), while RR4 is skewed in the opposite direction (18% positive).

**Table1DatasetSummary.py**

This Python script generates a summary table showing the number and percentage of positive and negative disease pairs across three different relative risk (RR) thresholds.
It produces the summary table **Table 1** using different RR thresholds (RR > 1, RR > 2, and RR > 4) for labeling.

**Table 2 (Performance Evaluation of Neural Network)**: Compares the DR-V and SNP-V approaches against state-of-the-art semi-supervised methods (Baseline-1) and a naïve similarity-based classifier (Baseline-2), of all three cases of RR thresholds.

**baseline1performanceRR1.py, baseline1performanceRR2.py, baseline1performanceRR4.py**

These three Python scripts generates the Baseline-1 performance metrics in **Table 2**. 

The following three Python scripts generates the Baseline-2 performance metrics in **Table 2**.

**w_matrix_performance_summarybaseline2RR1.py**

It evaluates the W-matrix similarity performance against the RR1 labeled dataset. It checks if the similarity matrix is symmetric and calculates predictive performance per disease row. 

**w_matrix_performance_summarybaseline2RR2.py**

It evaluates the W-matrix similarity performance against the RR2 labeled dataset. It checks if the similarity matrix is symmetric and calculates predictive performance per disease row. 

**w_matrix_performance_summarybaseline2RR4.py**

It evaluates the W-matrix similarity performance against the RR4 labeled dataset. It checks if the similarity matrix is symmetric and calculates predictive performance per disease row. 

The following Python script generates the DR-V performance metrics in **Table 2**.

**RR1diseaserowvectorsNeuralNetwork.py**

This Python script trains a Neural Network on disease-row vector features for binary classification of disease pairs using the RR > 1 labeling threshold. It first builds a full cosine similarity matrix W from the disease-SNP matrix (treating each disease as a vector of similarities to all other diseases), standardizes W, reduces each disease vector to 50 principal components via PCA, concatenates the 50 PCs from both diseases in a pair (creating 100-dimensional feature vectors for each pair of diseases), and then performs 10-fold stratified cross-validation plus a held-out test set evaluation.

To use the same pipeline for RR2 or RR4 thresholds of ground truth labeling, replace the RR1 labelled file "./disease_pairs_with_labels.csv" with the corresponding CSV label file, i.e. "./disease_pairs_with_labelsRR2.csv" for RR2 labelled file, "./disease_pairs_with_labelsRR4.csv" for RR4 labelled file. 

The following Python script generates the SNP-V performance metrics in **Table 2**.

**RR1SNPvectorsNeuralNetwork.py**

This Python script trains the Neural Network but uses a different feature representation: direct PCA on the original disease-SNP matrix (SNPs as features). It standardizes the full SNP matrix, reduces the SNP vector from 39,382 in dimension to 50 principal components per disease, concatenates the 50 PCs from both diseases in each pair (again 100 features total), and runs identical 10-fold stratified cross-validation plus test-set evaluation for RR > 1 disease-pair classification.

To use the same pipeline for RR2 or RR4 thresholds of ground truth labeling, replace the RR1 labelled file "./disease_pairs_with_labels.csv" with the corresponding CSV label file, i.e. "./disease_pairs_with_labelsRR2.csv" for RR2 labelled file, "./disease_pairs_with_labelsRR4.csv" for RR4 labelled file. 

**Table 3 (Random Forest Performance Evaluation)**: Compares the performance of DR-V and SNP-V approaches by the Random Forest classifier, of all three cases of RR thresholds. 

For **Table 3** DR-V Random Forest results: 

**ConcatscriptSD53CV10alldiseasepairRR1diseaserowvectorsTrainTestRF.py**

This Python script performs binary classification of disease pairs using the **RR > 1** labeling threshold. Using the disease row vector in a cosine similarity matrix W derived from the disease-SNP data, it applies PCA to reduce each disease row vector to 50 principal components, concatenates the two diseases’ reduced vectors (creating 100-dimensional features per pair), trains a Random Forest classifier, runs 10-fold stratified cross-validation, evaluates on a held-out test set, and saves performance metrics. 

**ConcatscriptSD51CV10alldiseasepairRR2diseaserowvectorsTrainTestRF.py**

This Python script performs binary classification of disease pairs using the **RR > 2** labeling threshold. Using the disease row vector in a cosine similarity matrix W derived from the disease-SNP data, it applies PCA to reduce each disease row vector to 50 principal components, concatenates the two diseases’ reduced vectors (creating 100-dimensional features per pair), trains a Random Forest classifier, runs 10-fold stratified cross-validation, evaluates on a held-out test set, and saves performance metrics.

**ConcatscriptSD52CV10alldiseasepairRR4diseaserowvectorsTrainTestRF.py**

This Python script performs binary classification of disease pairs using the **RR > 4** labeling threshold. Using the disease row vector in a cosine similarity matrix W derived from the disease-SNP data, it applies PCA to reduce each disease row vector to 50 principal components, concatenates the two diseases’ reduced vectors (creating 100-dimensional features per pair), trains a Random Forest classifier, runs 10-fold stratified cross-validation, evaluates on a held-out test set, and saves performance metrics.


For **Table 3** SNP-V Random Forest results: 

**scriptstoredSDVarPCA10063CV10alldiseasepairRR1SNPvectorsTrainTestEigen20Features50SNPseachConcatRF.py**

This Python script gives a Random Forest classification pipeline using PCA-reduced SNP vectors for disease-pair prediction with the **RR > 1** labeling threshold. It standardizes the full disease-SNP matrix, reduces the original SNP features to 50 principal components per disease via PCA, concatenates the two diseases’ 50-PC vectors to create 100-dimensional features per pair, trains a Random Forest classifier, runs 10-fold stratified cross-validation plus a held-out test set evaluation, and conducts feature importance analysis to identify and obtain the top 1,000 most relevant SNPs based on ranked eigenvectors by absolute values.

The results file "RF_Top20_PC_to_SNP_Analysis.csv" gives the 20 most important PCA features, including which disease and PC they correspond to, plus the top 50 SNPs loaded on each feature. 

The results file "Top_1000_SNPs_Ranked.csv" gives the ranked list of the 1,000 top SNPs.

Use the results file "Concat_CV10Fold_Metrics_SNP_PCA_Detailed_RF_RR1.csv" for further analysis of the Random Forest model performance by **RFtable3SNPVRR1.R**. 

**RFtable3SNPVRR1.R**

This R script computes the mean and standard deviation of six classification performance metrics (Accuracy, Precision, Recall, F1, ROC AUC, and PR AUC) from the 10-fold cross-validation results of a Random Forest model that was trained on PCA-reduced SNP vectors using the RR > 1 labeling threshold.

**RFtable3SNPVRR2.R**

This R script computes the mean and standard deviation of six classification performance metrics (Accuracy, Precision, Recall, F1, ROC AUC, and PR AUC) from the 10-fold cross-validation results of a Random Forest model that was trained on PCA-reduced SNP vectors using the RR > 2 labeling threshold.

**RFtable3SNPVRR4.R**

This R script computes the mean and standard deviation of six classification performance metrics (Accuracy, Precision, Recall, F1, ROC AUC, and PR AUC) from the 10-fold cross-validation results of a Random Forest model that was trained on PCA-reduced SNP vectors using the RR > 4 labeling threshold.

**Table 4 a)b)c) (Kolmogorov-Smirnov Test Results)**: Presents the KS-statistics and p-values from the two-sample KS test used to compare the distributions of association strengths between top and random SNPs, of all three cases of RR thresholds. 

To prepare for and generate **Table 4(a)** for the RR1 case, use the following R scripts: 

**952nonOverlapRandomvsTopSNPsRR1.R**

This R script creates non-overlapping random SNP sets for comparison with the 952 top SNPs for the RR1 threshold case. It loads the full disease-SNP matrix, removes the 952 top SNPs to create a clean non-overlapping pool, then randomly samples 952 SNPs overall plus two separate subsets (520 SNPs, each of which has at least one positive effect among 427 phenotypes and 434 SNPs, each of which has at least one negative effect among 427 phenotypes). It also performs statistical comparison (Kolmogorov-Smirnov test) between the top and random groups and prepares data for the absolute effect histograms. 

**topSNP952analysisRR1.R**

This R script analyzes the 952 top-ranked SNPs (based on the Random Forest feature importance) for the RR1 threshold case. It calculates per-SNP statistics on both original effect values and absolute effect magnitudes (excluding zeros), extracts associated phenotypes, and prepares the data for plotting the absolute effect distribution.

**952RandomSNPnonoverlapRR1.R**

This R script performs the similar per-SNP summary statistics analysis as the top-SNP script but on the 952 non-overlapping random SNPs. It computes mean, median, min, max, and SD for both original and absolute non-zero effect values and prepares the data needed for the random-SNP histogram.


**PosNegTopSNPsetRR1.R**

This R script processes the 952 top-ranked SNPs identified by the Random Forest model (with RR1 ground truth labels) and derives two separate groups: SNPs that show positive effects and SNPs that show negative effects. For each subgroup, it calculates mean effect magnitudes (excluding zeros), and identifies SNPs that exhibit both positive and negative effects on different diseases.

**520posSNPtopvsrandomRR1.R**

This R script compares the mean positive effects of “top positive SNPs” against a random set of non-overlapping SNPs. It keeps only positive effect values of SNPs, computes per-SNP mean positive effect, gives comparison tables, and runs a Kolmogorov-Smirnov test to compare the two distributions.

**434negSNPtopvsrandomRR1.R**

This R script compares the mean negative effects of “top negative SNPs” against a random set of non-overlapping SNPs. It mirrors the positive-SNP script, but it keeps only negative values, computes per-SNP mean negative effect, gives comparison tables, and runs a Kolmogorov-Smirnov test to compare the two distributions.

Similar sets of R scripts have been developed to prepare for and generate **Table 4(b)** for the RR2 case and **Table 4(c)** for the RR3 case. 


### Key Figures: 

**Figure 1 (Methodology Overview)**: Illustrates the pipeline and schematic description of proposed methods. 

<img width="644" height="362" alt="Figure1" src="https://github.com/user-attachments/assets/13c93142-88ac-4862-9d2c-97918407f2df" />

**Figure 2 (RR Distribution)**: A histogram showing the percentage of disease pairs at various relative risk thresholds.

**Fig2barplotRRdistributionVaryingThreshold.py**

This Python script generates and displays **Figure 2**, i.e. a bar plot titled:
"Percentage of RR Values Greater Than Different Thresholds"

**Figure 3. Pseudo codes for the proposed methods**.

**Figure 4: Neural Network Architecture**. 

**Figure 5-7: Histograms showing the distribution of mean absolute effect of SNPs (top SNPs vs. random SNPs), of all three cases of RR thresholds**. 

To generate **Figure 5**, run each of the following Python scripts (order does not matter)
(Histogram952_topvsrandom_snps_groups_abs_RR1.py, 
Histogram520_topvsrandom_pos_snps_groups_abs_RR1.py, 
Histogram434_topvsrandom_neg_snps_groups_abs_RR1.py), 
each of which will create a histogram comparing “Top” (red) vs “Random” (blue) SNP groupps with the distribution of mean absolute effect magnitudes of SNPs in three separate analyses (all with the **RR > 1** labeling threshold). Then use the Python script (Fig5Combined3HistogramsSNPsRR1.py) to combine and stack the three histograms vertically for **Figure 5**. 

To generate **Figure 6**,run each of the following Python scripts (order does not matter):
Histogram950_topvsrandom_snps_groups_abs_RR2.py → creates the histogram for overall comparison of 950 Top SNPs versus 950 Random SNPs;

Histogram509_topvsrandom_pos_snps_groups_abs_RR2.py → creates the histogram for the comparison of positive-effect SNPs only (509 Top vs 509 Random);

Histogram499_topvsrandom_neg_snps_groups_abs_RR2.py → creates the histogram for the comparison of negative-effect SNPs only (499 Top vs 499 Random),

all with the **RR > 2** labeling threshold. Then combine and stack the three individual histograms vertically for **Figure 6**. 

To generate **Figure 7**,run each of the following Python scripts (order does not matter):
Histogram950_topvsrandom_snps_groups_abs_RR4.py → creates the histogram for overall comparison of 950 Top SNPs versus 950 Random SNPs;

Histogram505_topvsrandom_pos_snps_groups_abs_RR4.py → creates the histogram for the comparison of positive-effect SNPs only (505 Top vs 505 Random);

Histogram500_topvsrandom_neg_snps_groups_abs_RR4.py → creates the histogram for the comparison of negative-effect SNPs only (500 Top vs 500 Random),

all with the **RR > 3** labeling threshold. Then combine and stack the three individual histograms vertically for **Figure 7**. 

## Key Findings:

The DR-V approach generally provided the best overall performance, suggesting that network-level data can compensate for information lost during SNP dimension reduction.

The SNP-V approach allows for unique biological interpretability by tracing back to the top SNPs that contribute most to comorbidity prediction.

Top-ranked SNPs identified by the model showed statistically significant antagonistic and synergistic patterns compared to random SNPs, potentially uncovering root causes for disease co-occurrence.
