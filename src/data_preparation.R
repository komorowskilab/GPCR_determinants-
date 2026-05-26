# Install necessary packages
if (!requireNamespace("rmcfs", quietly = TRUE)) install.packages("rmcfs")
if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("R.ROSETTA", quietly = TRUE)) {devtools::install_github("komorowskilab/R.ROSETTA")}
if (!requireNamespace("VisuNet", quietly = TRUE)) {devtools::install_github("komorowskilab/VisuNet")}
if (!requireNamespace("rio", quietly = TRUE)) install.packages("rio")

# Load necessary packages
library(rJava)
library(rmcfs)
library(devtools)
library(R.ROSETTA)
library(VisuNet)
library(dplyr)
library(Biostrings)
library(tidyr)
library(readxl)
library(caret)
library(rio)
library(purrr)



### ---------------------------LOAD THE DATA----------------------------------
data <- read.table("aligned.tsv", 
                   header = TRUE, sep = "\t", colClasses = "character", row.names = 1)
#validation <- read.table("validation_aligned.tsv", 
                         #header = TRUE, sep = "\t", colClasses = "character", row.names = 1)
#cluster_data <- read.table("clustering_data_aligned.tsv", 
                           #header = TRUE, sep = "\t", colClasses = "character", row.names = 1)



### ------------------ADD LABELS-----------------------------

labels <- read_excel("labels.xlsx") #Labels for adrenergic and octopamine
data$Receptor <- labels$LABEL[match(rownames(data), labels$SAMPLE)]
table(data$Receptor)


### ------------------SPLIT THE DATA-----------------------------

# Define samples that must be in training set (these are NCBI samples)
# Here they are human adrenergic, skow adrenergic and octopamine, pdumerili adrenergic and octopamine
forced_train_samples <- c("P08913", "XP_002734932", "XP_006823182", "APC23843", "APC23183")

# Separate samples into three groups
normal_samples <- data[!grepl("^M|^CHIMERA", rownames(data)), ]  # NCBI samples
M_samples <- data[grepl("^M", rownames(data)), ]  # Mutant samples
chimera_samples <- data[grepl("^CHIMERA", rownames(data)), ]  # Chimera samples

cat("Sample breakdown:\n")
cat("Normal samples (NCBI):", nrow(normal_samples), "\n")
cat("M samples (MUTANT):", nrow(M_samples), "\n")
cat("CHIMERA samples:", nrow(chimera_samples), "\n\n")

# Set seed for reproducibility
set.seed(38)



# --- HANDLE NCBI SAMPLES ---
# Extract forced samples from normal samples
forced_normal <- normal_samples[rownames(normal_samples) %in% forced_train_samples, ]
remaining_normal <- normal_samples[!rownames(normal_samples) %in% forced_train_samples, ]

cat("Forced to training (from NCBI):", nrow(forced_normal), "\n")
if(nrow(forced_normal) > 0) {
  cat("  Samples:", rownames(forced_normal), "\n")
}

cat("\n=== NCBI SAMPLES CLASS DISTRIBUTION ===\n")
cat("Total NCBI samples:", nrow(normal_samples), "\n")
print(table(normal_samples$Receptor))
print(prop.table(table(normal_samples$Receptor)))

# Split remaining normal samples stratified by Receptor class (80/20)
train_indices_normal <- createDataPartition(remaining_normal$Receptor, 
                                            p = 0.8,
                                            list = FALSE)
train_normal_random <- remaining_normal[train_indices_normal, ]
test_normal <- remaining_normal[-train_indices_normal, ]

# Combine forced + random for normal training set
train_normal <- rbind(forced_normal, train_normal_random)

cat("\nNCBI samples split:\n")
cat("  - Training NCBI:", nrow(train_normal), "(", nrow(forced_normal), "forced )\n")
print(table(train_normal$Receptor))
print(prop.table(table(train_normal$Receptor)))
cat("  - Test NCBI:", nrow(test_normal), "\n")
print(table(test_normal$Receptor))
print(prop.table(table(test_normal$Receptor)))

# Combine into final train and test sets
train_data <- train_normal  # ONLY NCBI samples in training
test_data <- rbind(test_normal, M_samples, chimera_samples)  # NCBI + ALL MUTANT + ALL CHIMERA in test

# Verify forced samples are in training
cat("\n=== VERIFICATION ===\n")
cat("Checking forced NCBI samples are in training set:\n")
for(sample in forced_train_samples) {
  in_train <- sample %in% rownames(train_data)
  in_test <- sample %in% rownames(test_data)
  cat(sprintf("  %s: Training=%s, Test=%s\n", sample, in_train, in_test))
}

# Verify distributions
cat("\n=== TRAINING SET ===\n")
cat("Total samples:", nrow(train_data), "\n")
cat("  - NCBI samples:", nrow(train_normal), "(", nrow(forced_normal), "forced )\n")
cat("Class distribution:\n")
print(table(train_data$Receptor))
print(prop.table(table(train_data$Receptor)))

cat("\n=== TEST SET ===\n")
cat("Total samples:", nrow(test_data), "\n")
cat("  - NCBI samples:", nrow(test_normal), "\n")
cat("  - M samples (MUTANT):", nrow(M_samples), "\n")
cat("  - CHIMERA samples:", nrow(chimera_samples), "\n")
cat("Class distribution:\n")
print(table(test_data$Receptor))
print(prop.table(table(test_data$Receptor)))

# Save the rownames
train_rownames <- data.frame(sample_name = rownames(train_data))
test_rownames <- data.frame(sample_name = rownames(test_data))
write.csv(train_rownames, "train_rownames.csv", row.names = FALSE)
write.csv(test_rownames, "test_rownames.csv", row.names = FALSE)




### ----------------DEFINE AND EXTRACT EXTRACELLULAR REGION------------------------
# This region should be manually defined from the alignment
EC_region <- c(
  paste0("P", 340:358),
  paste0("P", 408:432),
  paste0("P", 438:457),
  paste0("P", 509:516),
  paste0("P", 643:664),
  paste0("P", 1049:1068),
  paste0("P", 1078:1094),
  "Receptor"
)

#full_length_without_loops <- c(
#  paste0("P", 339:386),
#  paste0("P", 392:432),
#  paste0("P", 438:482),
#  paste0("P", 491:516),
#  paste0("P", 643:691),
#  paste0("P", 1028:1068),
#  paste0("P", 1078:1104),
#  "Receptor"
#)


# Keep only columns that exist in the dataframe
existing_cols_train <- intersect(EC_region, colnames(train_data))
existing_cols_test <- intersect(EC_region, colnames(test_data))

# Subset safely
train_data <- train_data[, existing_cols_train]
test_data <- test_data[, existing_cols_test]


###----------------- COLUMN NAMING TO ORIGINAL------------------------------
# First, get the real position names from the necessary excel file 
# Here, we want to change the column names to the original position names that 
# are consensus annotations on GPCRdb.


excel_path <- "position_conversion.xlsx" # should be adjusted if alignment changes

mapping <- excel_sheets(excel_path) %>%
  map_dfr(~ read_excel(excel_path, sheet = .x, col_names = TRUE,
                       col_types = "text"))

name_map <- setNames(mapping$`Position GPCRdb`, mapping$`Position Experimental`)

colnames(train_data) <- ifelse(colnames(train_data) %in% names(name_map), name_map[colnames(train_data)], colnames(train_data))
colnames(test_data) <- ifelse(colnames(test_data) %in% names(name_map), name_map[colnames(test_data)], colnames(test_data))


### ----------- REMOVE COLUMNS WITH ≥90% NA TO GET RID OF UNNECESSARY INFO-------- 

# Replace "-" with NA if needed
train_data[train_data == "-"] <- NA
test_data[test_data == "-"] <- NA

# Calculate NA proportion per column
na_prop_train <- colMeans(is.na(train_data))
na_prop_test  <- colMeans(is.na(test_data))

# Only remove columns that are >=90% NA in BOTH
cols_to_remove <- intersect(
  names(na_prop_train[na_prop_train >= 0.9]),
  names(na_prop_test[na_prop_test >= 0.9])
)

cat("Columns removed (>=90% NA in both):\n")
print(cols_to_remove)

# Remove same columns from both
train_data <- train_data[, !names(train_data) %in% cols_to_remove]
test_data  <- test_data[,  !names(test_data)  %in% cols_to_remove]

cat(length(cols_to_remove), "columns removed\n")

cluster_19 <- cluster_19[, colMeans(is.na(cluster_19)) < 0.9]

# Save the full datasets
write.csv(train_data, "train_set.csv", row.names = TRUE)
write.csv(test_data, "test_set.csv", row.names = TRUE)

