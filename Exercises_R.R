############################################
# R BASICS FOR TRANSCRIPTOMICS
############################################

### SECTION 1 — Calculator Basics

# Exercise 1:
# Perform the following operations:
# 1. Add 5 + 3
# 2. Multiply 7 * 6
# 3. Divide 20 / 4

# Write your code below:


# Exercise 2:
# Compute the following expression:
# (10 + 5) * 2


############################################
### SECTION 2 — Variables

# Exercise 3:
# Create two variables:
# a = 15
# b = 4
# Then compute:
# - sum
# - difference
# - product

# Write your code below:


############################################
### SECTION 3 — Vectors (Gene Expression)

# Exercise 4:
# Create a vector representing gene expression:
gene_expr <- c(10, 50, 30, 80, 25)

# Tasks:
# - Calculate mean expression
# - Find max value
# - Find min value


# Exercise 5:
# Extract values greater than 40 (overexpressed genes)
gene_expr[gene_expr > 40]
genes[gene_expr > 40]


############################################
### SECTION 4 — Data Frames (RNA-seq style)

# Exercise 6:
df <- data.frame(
  gene = c("TP53", "BRCA1", "MYC", "EGFR", "PTEN"),
  control = c(10, 50, 30, 80, 25),
  treated = c(20, 55, 90, 70, 30)
)

# Tasks:
# - View dataset
# - Calculate fold change (treated / control)
# - Add fold_change column

View(df)
df$treated / df$control
df$fold_change <- df$treated / df$control


############################################
### SECTION 5 — Boxplot (Transcriptomics Visualization)

# Exercise 7:
# Create a boxplot comparing control vs treated

boxplot(df$control, df$treated,
        names = c("Control", "Treated"),
        main = "Gene Expression Comparison",
        ylab = "Expression Level")

# Question:
# - Which condition shows higher variability?


############################################
# END OF EXERCISES
############################################