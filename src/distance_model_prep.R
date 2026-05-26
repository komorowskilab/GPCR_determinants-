# Install necessary packages
if (!requireNamespace("rmcfs", quietly = TRUE)) install.packages("rmcfs")
if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("R.ROSETTA", quietly = TRUE)) {devtools::install_github("komorowskilab/R.ROSETTA")}
if (!requireNamespace("VisuNet", quietly = TRUE)) {devtools::install_github("komorowskilab/VisuNet")}

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
library(purrr)




### ---------------------------LOAD THE DATA------------------------------------

all_data <- read.table("aligned_distances.tsv",
                          header = TRUE, sep = "\t", fill = TRUE, row.names = 1)

train_rows <- read.csv("train_rownames.csv", header = TRUE)[,1]
test_rows  <- read.csv("test_rownames.csv", header = TRUE)[,1]

data_train <- all_data[rownames(all_data) %in% train_rows, ]
data_test  <- all_data[rownames(all_data) %in% test_rows, ]

data_validation <- read.table("validation_mutation_distance_aligned.tsv",
                            header = TRUE, sep = "\t", fill = TRUE, row.names = 1)


### ------------------ADD LABELS FOR TRAINING DATA -----------------------------

labels <- read_excel("labels.xlsx")
data_train$Receptor <- labels$LABEL[match(rownames(data_train), labels$SAMPLE)]
data_test$Receptor <- labels$LABEL[match(rownames(data_test), labels$SAMPLE)]


# Check classes
table(data_train$Receptor)
table(data_test$Receptor)



### ------------------ExtraCellular region (EC) ------------------------

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
existing_cols_train <- intersect(EC_region, colnames(data_train))
existing_cols_test <- intersect(EC_region, colnames(data_test))
#existing_cols_validation <- intersect(EC_region, colnames(data_validation))

# Subset safely
data_train <- data_train[, existing_cols_train]
data_test <- data_test[, existing_cols_test]
#data_validation <- data_validation[, existing_cols_validation]




###----------------- COLUMN NAMING TO ORIGINAL------------------------------
#First, get the real position names from the necessary excel file 
#Here, we want to change the column names to the original position names that 
#are consensus annotations on GPCRdb.

excel_path <- "position_conversion.xlsx"

mapping <- excel_sheets(excel_path) %>%
  map_dfr(~ read_excel(excel_path, sheet = .x, col_names = TRUE,
                       col_types = "text"))

name_map <- setNames(mapping$`Position GPCRdb`, mapping$`Position Experimental`)

colnames(data_train) <- ifelse(colnames(data_train) %in% names(name_map), name_map[colnames(data_train)], colnames(data_train))
colnames(data_test) <- ifelse(colnames(data_test) %in% names(name_map), name_map[colnames(data_test)], colnames(data_test))
colnames(data_validation) <- ifelse(colnames(data_validation) %in% names(name_map), name_map[colnames(data_validation)], colnames(data_validation))



### ----------- REMOVE COLUMNS WITH ≥90% NA TO GET RID OF UNNECESSARY INFO-------- 

# Replace "-" with NA if needed
data_train[data_train == "-"] <- NA
data_test[data_test == "-"] <- NA

# Calculate NA proportion per column
na_prop_train <- colMeans(is.na(data_train))
na_prop_test  <- colMeans(is.na(data_test))

# Only remove columns that are >=90% NA in BOTH
cols_to_remove <- intersect(
  names(na_prop_train[na_prop_train >= 0.9]),
  names(na_prop_test[na_prop_test >= 0.9])
)

cat("Columns removed (>=90% NA in both):\n")
print(cols_to_remove)

# Remove same columns from both
data_train <- data_train[, !names(data_train) %in% cols_to_remove]
data_test  <- data_test[,  !names(data_test)  %in% cols_to_remove]

cat(length(cols_to_remove), "columns removed\n")



##----------------------Plot the Binding Site Distance Distribution---------
# Exclude last column
all_distances <- as.numeric(unlist(data_train[, -ncol(data_train)]))
# Remove NAs
all_distances <- all_distances[!is.na(all_distances)]

# Calculate statistics
mean_val <- round(mean(all_distances), 2)
max_val <- round(max(all_distances), 2)
min_val <- round(min(all_distances), 2)
std_val <- round(sd(all_distances), 2)

# Open PNG device
png("binding_site_distance_distribution.png", width = 800, height = 600, res = 100)

# Plot
hist(all_distances, breaks = 50, col = "steelblue", 
     main = "Binding Site Distance Distribution",
     xlab = "Distance (Å)", 
     ylab = "Frequency")

# Add statistics box
legend("topright", 
       legend = c(paste("Mean:", mean_val),
                  paste("Max:", max_val),
                  paste("Min:", min_val),
                  paste("Std:", std_val)),
       bty = "n",
       cex = 0.9)

# Close the device (this actually saves the file)
dev.off()

cat("Plot saved to: binding_site_distance_distribution.png\n")



###-----------------Gaussian noise to improve the signal of the closer regions to ligand----------------
set.seed(123)
noise_for_col <- function(col) {
  col_mode <- round(median(col, na.rm = TRUE), 1)
  
  if (is.na(col_mode) || col_mode <= 9.4) return(col)
  
  if (col_mode <= 14.94) {
    noise_pool <- seq(0, 0.04, length.out = 1e6)
  } else if (col_mode <= 20.48) {
    noise_pool <- seq(0, 0.10, length.out = 1e6)
  } else {
    noise_pool <- seq(0, 0.20, length.out = 1e6)
  }
  
  noise <- sample(noise_pool, length(col), replace = TRUE)
  return(col + (col * noise))
}

# Apply to train
data_train_noisy <- data_train
numeric_cols_train <- sapply(data_train_noisy, is.numeric)
for (colname in names(data_train_noisy)[numeric_cols_train]) {
  data_train_noisy[[colname]] <- noise_for_col(data_train_noisy[[colname]])
}

# Apply to test
data_test_noisy <- data_test
numeric_cols_test <- sapply(data_test_noisy, is.numeric)
for (colname in names(data_test_noisy)[numeric_cols_test]) {
  data_test_noisy[[colname]] <- noise_for_col(data_test_noisy[[colname]])
}

# Apply to validation
#data_validation_noisy <- data_validation
#numeric_cols_test <- sapply(data_validation_noisy, is.numeric)
#for (colname in names(data_validation_noisy)[numeric_cols_test]) {
#  data_validation_noisy[[colname]] <- noise_for_col(data_validation_noisy[[colname]])
#}

# Verify
cat("Are train and noisy identical?", identical(data_train, data_train_noisy), "\n")
cat("Are test and noisy identical?", identical(data_test, data_test_noisy), "\n")
#cat("Are validation and noisy identical?", identical(data_validation, data_validation_noisy), "\n")


##----------------------------------------------------------------------####

# Replace remaining NA with 22
data_train_noisy[is.na(data_train_noisy)] <- 22
data_test_noisy[is.na(data_test_noisy)]   <- 22
#data_validation_noisy[is.na(data_validation_noisy)]   <- 22


###-----------------------SAVE THE OUTPUT-------------------------------------
write.csv(data_test_noisy, "data_test.csv", row.names = TRUE)
write.csv(data_train_noisy, "data_train.csv", row.names = TRUE)
#write.csv(data_validation_noisy, "data_validation.csv", row.names = TRUE)



