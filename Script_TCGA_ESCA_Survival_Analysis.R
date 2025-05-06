#-------------------------------------------------------------------------------
# Survival Analysis - TCGA Data (RSEM log2(x+1) normalized values)
#-------------------------------------------------------------------------------

# # Install required libraries
# install.packages("dplyr")
# install.packages("survival")
# install.packages("survminer")

# Load required libraries
library(dplyr)
library(survival)
library(survminer)


# Part 1. Count Matrix and Metadata.
# ------------------------------------------------------------------------------

# Load Count Matrix.
count_matrix <- read.delim(file.choose(), header = TRUE, row.names = 1, sep = "\t")
colnames(count_matrix) <- gsub("\\.", "-", colnames(count_matrix))

# Load Phenotype Data.
clinical_data <- read.delim(file.choose(), header = TRUE, row.names = 1)
survival_data <- read.delim(file.choose(), header = TRUE, row.names = 1)
survival_data <- survival_data[, -ncol(survival_data)]

# Merge clinical and survival data
merged_data <- merge(clinical_data, survival_data, by = "row.names", all = TRUE)
rownames(merged_data) <- merged_data$Row.names
merged_data$Row.names <- NULL

# Filter for Primary Tumor and/or Metastatic samples
merged_data <- merged_data %>% filter(sample_type %in% c("Primary Tumor"))

# Identify common patients
common_patients <- intersect(colnames(count_matrix), rownames(merged_data))
merged_data <- merged_data %>% filter(rownames(merged_data) %in% common_patients)
count_matrix <- count_matrix[, common_patients]


# Part 2. Gene signature Analysis.
# ------------------------------------------------------------------------------

# a) Extract gene expression values from count_matrix
gene_sig_exp <- count_matrix["PDCD1", ] 

# Ensure gene_exp is a numeric vector and sample IDs match between merged_data and count_matrix
gene_sig_exp <- as.numeric(gene_sig_exp)

# Match and correctly order the gene expression values with merged_data's sample IDs
merged_data$gene_sig_exp <- gene_sig_exp[match(rownames(merged_data), colnames(count_matrix))]

# Verify that gene_exp is correctly added
print(summary(merged_data$gene_sig_exp))
head(merged_data)

# ------------------------------------------------------------------------------
# b) Extract more than one gene expression values from count_matrix
# Define your gene signature
gene_sig_exp <- c("PDCD1", "CD274", "CTLA4")

# Check if all genes are present in the count matrix
missing_genes <- gene_sig_exp[!gene_sig_exp %in% rownames(count_matrix)]
if (length(missing_genes) > 0) {
  stop(paste("Missing genes in count matrix:", paste(missing_genes, collapse = ", ")))
}

# Extract and compute average expression of the signature per sample
gene_sig_matrix <- count_matrix[gene_sig_exp, ]
gene_sig_exp <- colMeans(gene_sig_matrix, na.rm = TRUE)

# Add to merged_data
merged_data$gene_sig_exp <- gene_sig_exp[match(rownames(merged_data), names(gene_sig_exp))]

# Verify expression vector
print(summary(merged_data$gene_sig_exp))
head(merged_data)


# Part 3. Cut-off Median
# ------------------------------------------------------------------------------

# Calculate the median expression of the gene signature
median_gene_signature_exp <- median(merged_data$gene_sig_exp, na.rm = TRUE)

# Split patients into low and high expression groups based on the median
merged_data$expression_group <- ifelse(merged_data$gene_sig_exp > 
                                         median_gene_signature_exp, "High", "Low")

# Check the distribution of groups
table(merged_data$expression_group)

# Verify the first few rows
head(merged_data)


# Part 4. Overall Survival Analysis
# ------------------------------------------------------------------------------

# Fit survival model
surv_object <- Surv(time = merged_data$OS.time, event = merged_data$OS)
fit <- survfit(surv_object ~ expression_group, data = merged_data)

# Plot Kaplan-Meier survival curves
plot1 <- ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE, 
                    title = expression(OS ~ "- gene expression"),
                    xlab = "Time (days)", ylab = "Overall Survival Probability",
                    legend.labs = c("High Expression", "Low Expression"),
                    palette = c("#d1495b", "#2e4057"),
                    break.time.by = 365)

# Rotate x-axis labels diagonally
plot1$plot <- plot1$plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot1$table <- plot1$table + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show survival plot
plot1


# Part 5. Disease Specific Survival Analysis
# ------------------------------------------------------------------------------

# Fit survival model
surv_object <- Surv(time = merged_data$DSS.time, event = merged_data$DSS)
fit <- survfit(surv_object ~ expression_group, data = merged_data)

# Plot Kaplan-Meier survival curves
plot2 <- ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE, 
                    title = expression(DSS ~ "- gene expression"),
                    xlab = "Time (days)", ylab = "Disease Specific Survival Probability",
                    legend.labs = c("High Expression", "Low Expression"),
                    palette = c("#d1495b", "#2e4057"),
                    break.time.by = 365)

# Rotate x-axis labels diagonally
plot2$plot <- plot2$plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot2$table <- plot2$table + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show survival plot
plot2

# End of script.
################################################################################