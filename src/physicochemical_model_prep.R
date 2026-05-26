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
library(caret)



### ---------------------------LOAD THE DATA----------------------------------
train_data <- read.table("train_set.csv", 
                           header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)
test_data <- read.table("test_set.csv", 
                         header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)
#validation_data <- read.table("validation_set.csv", 
                              #header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)

###--------------------START ADDING PROPERTIES---------------------------

AA <- c("A","R","N","D","C","Q","E","G","H","I",
        "L","K","M","F","P","S","T","W","Y","V")

aa_properties <- data.frame(
  AA       = AA,
  polarity = c("Nonpolar","Polar","Polar","Polar","Polar","Polar","Polar","Nonpolar","Polar","Nonpolar","Nonpolar","Polar","Nonpolar","Nonpolar","Nonpolar","Polar","Polar","Nonpolar","Polar","Nonpolar"),
  charge   = c("Neutral","Positive","Neutral","Negative","Neutral","Neutral","Negative","Neutral","Positive","Neutral","Neutral","Positive","Neutral","Neutral","Neutral","Neutral","Neutral","Neutral","Neutral","Neutral"),
  aromatic = c("Non-aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Non-Aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Aromatic","Non-aromatic","Non-aromatic","Non-aromatic","Aromatic","Aromatic","Non-aromatic"),
  hydrophob= c(1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2),
  ncharge  = c(0,1,0,-1,0,0,-1,0,1,0,0,1,0,0,0,0,0,0,0,0),
  surface  = c(129.0,274.0,195.0,193.0,167.0,225.0,223.0,104.0,224.0,197.0,201.0,236.0,224.0,240.0,159.0,155.0,172.0,285.0,263.0,174.0),
  pka      = c(2.34,2.17,2.02,1.88,1.96,2.17,2.19,2.34,1.82,2.36,2.36,2.18,2.28,1.83,1.99,2.21,2.09,2.83,2.20,2.32),
  pkb      = c(9.69,9.04,8.80,9.60,10.28,9.13,9.67,9.60,9.17,9.60,9.60,8.95,9.21,9.13,10.60,9.15,9.10,9.39,9.11,9.62),
  mweight  = c(89.10,174.20,132.12,133.11,121.16,146.15,147.13,75.07,155.16,131.18,131.18,146.19,149.21,165.19,115.13,105.09,119.12,204.23,181.19,117.15),
  flex     = c(0.783,0.807,0.799,0.822,0.785,0.817,0.826,0.784,0.777,0.776,0.783,0.834,0.806,0.774,0.809,0.811,0.795,0.796,0.788,0.781),
  vdW      = c(13.7,64.9,32.5,30.0,25.0,42.7,40.2,3.5,45.1,44.4,44.4,51.1,45.0,56.1,30.7,18.3,28.5,75.1,61.6,34.1),
  apol     = c(5.34,11.85,7.72,7.25,7.07,8.88,8.52,4.18,10.46,8.83,8.82,9.74,9.68,12.76,7.47,5.83,7.04,16.89,13.64,7.64),
  EN       = c(4.47,4.31,4.87,4.82,4.62,4.61,4.77,4.56,4.19,4.44,4.44,4.33,3.97,4.82,4.60,4.70,4.43,4.31,4.61,4.43),
  Stot     = c(133,221,160,158,150,184,179,112,196,189,188,203,194,212,159,137,159,247,218,168),
  Spol     = c(76,132,121,118,74,126,121,79,98,71,70,96,71,73,60,91,91,85,95,70),
  Snp      = c(58,90,39,41,76,58,59,33,98,118,118,108,125,140,100,47,69,163,123,98),
  HDONR    = c(0,4,2,1,0,2,1,0,1,0,0,2,0,0,0,1,1,1,1,0),
  HACCR    = c(0,3,3,4,0,3,4,0,1,0,0,1,0,0,0,2,2,0,2,0),
  Chpos    = c(0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0),
  Chneg    = c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
  logP     = c(-3.12,-4.79,-3.63,-2.43,-2.35,-3.46,-2.72,-3.21,-3.73,-1.76,-1.67,-3.42,-1.73,-1.56,-2.41,-2.74,-2.43,-1.57,-2.22,-2.29),
  stringsAsFactors = FALSE
)

# ---- FUNCTION TO ADD ALL PROPERTIES ----
add_properties <- function(df, original_cols, aa_properties) {
  prop_names <- names(aa_properties)[names(aa_properties) != "AA"]
  for (pos in original_cols) {
    aa_idx <- match(df[[pos]], aa_properties$AA)
    for (prop in prop_names) {
      df[[paste0(pos, ".", prop)]] <- aa_properties[[prop]][aa_idx]
    }
  }
  return(df)
}

# ---- APPLY ----
original_cols_train <- names(train_data)[names(train_data) != "Receptor"]
original_cols_test  <- names(test_data)[names(test_data) != "Receptor"]

train_data <- add_properties(train_data, original_cols_train, aa_properties)
test_data  <- add_properties(test_data,  original_cols_test,  aa_properties)

###---------------------CLEAN THE TABLE----------------------------------------
# Clean the names of the aminoacids
# First rows are for  the real original excel names
cleaned_train_data <- train_data[ , !grepl("^\\d+\\.\\d+X\\d+$", names(train_data)) ]
cleaned_test_data <- test_data[ , !grepl("^\\d+\\.\\d+X\\d+$", names(test_data)) ]
#cleaned_validation_data <- validation_data[ , !grepl("^\\d+\\.\\d+X\\d+$", names(validation_data)) ]

# Move "receptor" column to the end in cleaned_features_data
cleaned_train_data <- cleaned_train_data[ , c(setdiff(names(cleaned_train_data), "Receptor"), "Receptor") ]
cleaned_test_data <- cleaned_test_data[ , c(setdiff(names(cleaned_test_data), "Receptor"), "Receptor") ]

###-----------------------SAVE THE OUTPUT-------------------------------------
# Write csv for the data that is used in mcfs
write.csv(cleaned_train_data, "train_data_19.csv", row.names = TRUE)
write.csv(cleaned_test_data, "test_data_19.csv", row.names = TRUE)
#write.csv(cleaned_validation_data, "validation_data_19.csv", row.names = TRUE)



