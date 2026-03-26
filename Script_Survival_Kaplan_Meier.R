#------------------------------------------------------------------------------#
# Script Analysis Survival TCGA-GBMLGG Data - RSEM values Histological Type
#------------------------------------------------------------------------------#
# log2(x+1) transformed RSEM normalized count.

# ============================================
#  Load Required Libraries
# ============================================
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)


# Set Working Directory
setwd()

#------------------------------- ALL PATIENTS ---------------------------------#
# ==============================================================================
#  Part 1. Load Count Matrix and Metadata.
# ==============================================================================

# Load Count Matrix.
count_matrix <- read.delim(file.choose(), header = TRUE, row.names = 1, sep = "\t")
colnames(count_matrix)[1:5]
colnames(count_matrix) <- gsub("\\.", "-", colnames(count_matrix)) # Replace dots with hyphens
colnames(count_matrix)[1:5]

# Load Phenotype Data.
clinical_data <- read.delim(file.choose(), header = TRUE, row.names = 1)
colnames(clinical_data) <- gsub("^X_", "", colnames(clinical_data))

survival_data <- read.delim(file.choose(), header = TRUE, row.names = 1)
colnames(survival_data) <- gsub("^X_", "", colnames(survival_data))
survival_data <- survival_data[, -ncol(survival_data)]

# Merge the clinical_data and survival_data
merged_data <- merge(clinical_data, survival_data, by = "row.names", all = TRUE)

# remove(clinical_data, survival_data)
# gc()

# Set rownames back
rownames(merged_data) <- merged_data$Row.names
merged_data$Row.names <- NULL  # Remove the extra row names column

# Identify Common Patients
common_patients <- intersect(colnames(count_matrix), rownames(merged_data))
length(common_patients)

# Filter clinical_data_bis to include only Primary Tumor samples
merged_data <- merged_data %>%
  filter(sample_type == "Primary Tumor") # definition = "Primary solid Tumor"

# Filter clinical_data to include only common patients
merged_data <- merged_data %>%
  filter(rownames(merged_data) %in% common_patients)
nrow(merged_data)

# Identify Common Patients II
common_patients <- intersect(colnames(count_matrix), rownames(merged_data))
length(common_patients)

# Filter count_matrix to include only common patients
count_matrix <- count_matrix[, common_patients]

# Check the resulting data
print(nrow(merged_data))  # Number of rows in clinical data
print(ncol(count_matrix))  # Number of columns in count matrix (patients)



# ==============================================================================
# Part 2. Gene Expression Analysis.
# ==============================================================================

# i) Gene Expression Analysis - Single Analysis.
#-------------------------------------------------------------------------------
# Extract gene expression values from count_matrix
label <- "IL4I1"
gene_sig_exp <- count_matrix[label, ]

# Ensure gene_exp is a numeric vector and sample IDs match between merged_data and count_matrix
gene_sig_exp <- as.numeric(gene_sig_exp)

# Match and correctly order the AHR expression values with merged_data's sample IDs
merged_data$gene_sig_exp <- gene_sig_exp[match(rownames(merged_data), colnames(count_matrix))]

# Verify that gene_exp is correctly added
print(summary(merged_data$gene_sig_exp))
head(merged_data)



# ==============================================================================
# Part 3. Cut-off Median.
# ==============================================================================

# Calculate the median of gene or gene signature values
#-------------------------------------------------------------------------------
median_gene_signature_exp <- median(merged_data$gene_sig_exp, na.rm = TRUE)

# Split patients into low and high expression groups based on the median
merged_data$expression_group <- ifelse(merged_data$gene_sig_exp > median_gene_signature_exp, "High", "Low")

# Check the distribution of groups
table(merged_data$expression_group)

# Verify the first few rows
head(merged_data)



# ==============================================================================
# Part 4. Overall Survival.
# ==============================================================================

# For Overall Survival
merged_data$time  <- merged_data$OS.time
merged_data$event <- merged_data$OS
endpoint_label <- "OS"

# Fit the survival model
surv_object <- Surv(time = merged_data$time, event = merged_data$event)

# Perform Kaplan-Meier survival analysis
fit <- survfit(surv_object ~ expression_group, data = merged_data)
print(fit)

# Plot the Kaplan-Meier survival curves
plot <- ggsurvplot(
  fit,
  data = merged_data,
  pval = TRUE,
  risk.table = TRUE,
  title = bquote(.(endpoint_label) ~ "by" ~ italic(.(label)) ~ "Expression"),
  xlab = "Time (days)",
  ylab = paste(endpoint_label, "probability"),
  legend.labs = c("High Expression", "Low Expression"),
  palette = c("#d1495b", "#2e4057")
)

plot

#------------------------------- GBM PATIENTS ---------------------------------#

# ==============================================================================
# All pipeline only for GBM samples
# ==============================================================================

gbm_data <- merged_data %>%
  filter(primary_disease == "glioblastoma multiforme")

# Calculate the median of gene or gene signature values
median_gene_signature_exp <- median(gbm_data$gene_sig_exp, na.rm = TRUE)

# Split patients into low and high expression groups based on the median
gbm_data$expression_group <- ifelse(gbm_data$gene_sig_exp > median_gene_signature_exp, "High", "Low")

# Check the distribution of groups
table(gbm_data$expression_group)

# Verify the first few rows
head(gbm_data)

surv_object_gbm <- Surv(time = gbm_data$time, event = gbm_data$event)

fit_gbm <- survfit(surv_object_gbm ~ expression_group, data = gbm_data)

plot_gbm <- ggsurvplot(
  fit_gbm,
  data = gbm_data,
  pval = TRUE,
  risk.table = TRUE,
  title = bquote("OS in GBM by" ~ italic(.(label)) ~ "Expression"),
  xlab = "Time (days)",
  ylab = "Survival probability",
  legend.labs = c("High", "Low"),
  palette = c("#d1495b", "#2e4057")
)

plot_gbm


#------------------------------- LGG PATIENTS ---------------------------------#

# ==============================================================================
# All pipeline only for LGG samples
# ==============================================================================

lgg_data <- merged_data %>%
  filter(primary_disease == "brain lower grade glioma")

# Calculate the median of gene or gene signature values
median_gene_signature_exp <- median(lgg_data$gene_sig_exp, na.rm = TRUE)

# Split patients into low and high expression groups based on the median
lgg_data$expression_group <- ifelse(lgg_data$gene_sig_exp > median_gene_signature_exp, "High", "Low")

# Check the distribution of groups
table(lgg_data$expression_group)

# Verify the first few rows
head(lgg_data)

surv_object_lgg <- Surv(time = lgg_data$time, event = lgg_data$event)

fit_lgg <- survfit(surv_object_lgg ~ expression_group, data = lgg_data)

plot_lgg <- ggsurvplot(
  fit_lgg,
  data = lgg_data,
  pval = TRUE,
  risk.table = TRUE,
  title = bquote("OS in LGG by" ~ italic(.(label)) ~ "Expression"),
  xlab = "Time (days)",
  ylab = "Survival probability",
  legend.labs = c("High", "Low"),
  palette = c("#d1495b", "#2e4057")
)

plot_lgg

#------------------------------ GENE SIGNATURE --------------------------------#


# Gene Signature (a set of genes upregulated in some specific cell)

# Gene Signature Analysis
#-------------------------------------------------------------------------------
# Select a file .csv to analyze.
gene_sig_file <- read.csv(file.choose(), header = FALSE, stringsAsFactors = FALSE)
gene_signature <- gene_sig_file[[1]]
label <- "mregDCs Gene Signature"

# Check the length of the String.
length(gene_signature)

# Check which genes in the signature are found in the expression data
found_genes <- gene_signature[gene_signature %in% rownames(count_matrix)]
length(found_genes)

# Compute the mean expression for the found genes in the dataset
gene_signature_exp <- colMeans(count_matrix[found_genes, , drop = FALSE], na.rm = TRUE)

# Match and correctly order the gene expression values with merged_data's sample IDs
merged_data$gene_sig_exp <- gene_signature_exp[match(rownames(merged_data), colnames(count_matrix))]
head(merged_data)


# ==============================================================================
#  Cut-off Median
# ==============================================================================

# Calculate the median of gene or gene signature values
#-------------------------------------------------------------------------------
median_gene_signature_exp <- median(merged_data$gene_sig_exp, na.rm = TRUE)

# Split patients into low and high expression groups based on the median
merged_data$expression_group <- ifelse(merged_data$gene_sig_exp > median_gene_signature_exp, "High", "Low")

# Check the distribution of groups
table(merged_data$expression_group)

# Verify the first few rows
head(merged_data)


# ==============================================================================
# Overall Survival.
# ==============================================================================

# For Overall Survival
merged_data$time  <- merged_data$OS.time
merged_data$event <- merged_data$OS
endpoint_label <- "OS"

# Fit the survival model
surv_object <- Surv(time = merged_data$time, event = merged_data$event)

# Perform Kaplan-Meier survival analysis
fit <- survfit(surv_object ~ expression_group, data = merged_data)
print(fit)

# Plot the Kaplan-Meier survival curves
plot <- ggsurvplot(
  fit,
  data = merged_data,
  pval = TRUE,
  risk.table = TRUE,
  title = bquote(.(endpoint_label) ~ "by" ~ italic(.(label)) ~ "Expression"),
  xlab = "Time (days)",
  ylab = paste(endpoint_label, "probability"),
  legend.labs = c("High Expression", "Low Expression"),
  palette = c("#d1495b", "#2e4057")
)

plot

# End of script.