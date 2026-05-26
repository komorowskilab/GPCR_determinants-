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
library(csVisuNet)
library(dplyr)
library(Biostrings)
library(tidyr)
library(readxl)
library(rio)


###------------- LOAD THE DATA-------------------------

train_data <- read.table("/Users/dogya203/Desktop/thesis_folder_final_may/train_set.csv", 
                         header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)
test_data <- read.table("/Users/dogya203/Desktop/thesis_folder_final_may/test_set.csv", 
                        header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)
#validation_data <- read.table("/Users/dogya203/Desktop/thesis_folder_final_may/validation_set.csv", 
                              #header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)
#cluster_19 <- read.table("/Users/dogya203/Desktop/thesis_folder_final_may/clustering/cluster_data.csv", 
                         #header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)


# ----Feature Selection using RMCFS----

mcfs_train <- mcfs(Receptor ~ ., train_data, 
                                 projections = 20000, 
                                 projectionSize = 0.1, 
                                 splits = 5, 
                                 balance = 2,
                                 splitSetSize = 500, 
                                 cutoffPermutations = 30, 
                                 finalCV = TRUE,
                                 seed = 42)             # consistent results 
               

significant_mcfs_train <- mcfs_train$RI[1:mcfs_train$cutoff_value,]  


# ---- Save results ----
# Feature ranking
write.csv(mcfs_train$RI, "mcfs_feature_importance.csv", row.names = FALSE)

# Feature interdependencies
write.csv(mcfs_train$ID, "mcfs_feature_interdependencies.csv", row.names = FALSE)

# Cross-validation results
write.csv(mcfs_train$cv_accuracy, "mcfs_cv_accuracy.csv", row.names = FALSE)

# Significant features
write.csv(significant_mcfs_train, "mcfs_significant_features.csv", row.names = FALSE)


# ---- Plot distances ----
png(file.path("/Desktop", "mcfs_distances.png"), 
    width = 1200, height = 1000, res = 150)
plot(mcfs_train, type = "distances")
dev.off()
cat("✓ Saved distance plot\n")




# ---- Prepare filtered dataset for R.ROSETTA ----
selected_features <- as.character(significant_mcfs_train$attribute) # extract significant feature names from 'significant_result'
selected_features <- c(selected_features, "Receptor")      # add "Receptor" to the list of selected features
mcfs_selected_train <- train_data[, selected_features]
mcfs_selected_train$Receptor <- factor(mcfs_selected_train$Receptor, 
                                                   levels = c("Adrenergic", "Octopamine"))  



# Reverse columns, then move "Receptor" to the end, if needed
#cols <- rev(names(mcfs_selected_train))
#cols_reordered <- c(cols[cols != "Receptor"], "Receptor")
#mcfs_selected_train_reordered <- mcfs_selected_train[, cols_reordered]


# Run R.ROSETTA
out <- rosetta(mcfs_selected_train, 
                 roc =TRUE, 
                 underSample=TRUE,  
                 discrete =TRUE,
                 clroc = "Adrenergic",
                 reducer = "Johnson", 
                 JohnsonParam = list(Modulo=TRUE, BRT=TRUE, BRTprec=0.99, Precompute=TRUE, Approximate=TRUE, Fraction=0.99),)

rules <- out$main
rules_significant <- rules[rules$pValue < 0.05, ]


# Write to a plain text file (tab-separated)
write.table(
  rules,
  file = "rules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Count the number of rules for each decision class (Adrenergic and Octopamine)
rule_count <- table(rules$decision)


# Model evaluation
quality_model <- out$quality


# ROC curve
png("roc_curve_model.png", width = 800, height = 600, res = 150)
plotMeanROC(out)
dev.off()


# Test on external data
# Store decision
test_data_decision <- test_data$Receptor
# Remove decision from table
test_data <- test_data[,-length(test_data)]
# Test
decisions <- predictClass(test_data, rules, discrete=TRUE, normalize = TRUE, normalizeMethod="rss", validate = TRUE, defClass = test_data_decision)
table(Actual = decisions$out$currentClass, Predicted = decisions$out$predictedClass)
table(decisions$out$currentClass == decisions$out$predictedClass)

# To see the support set
recal_data_aa <- recalculateRules(train_data, rules_significant, discrete=TRUE)
View(recal_data_aa[str_detect((str_split(recal_data_aa$supportSetLHS,',')),'JAQ70295'),])


# Reclassification on training data
train_data_decision <- train_data$Receptor
train_data <- train_data[,-length(train_data)]
decisions_train <- predictClass(train_data, rules, discrete=TRUE, normalize = TRUE, normalizeMethod="rss", validate = TRUE, defClass = train_data_decision)
table(Actual = decisions_train$out$currentClass, Predicted = decisions_train$out$predictedClass)
table(decisions_train$out$currentClass == decisions_train$out$predictedClass)


# Test unseen mutant (validaton) data
decisions_validation_aa <- predictClass(validation_data, rules, discrete=TRUE, normalize = TRUE, normalizeMethod="rss", validate = FALSE)



## ----Clustering of models rules only for isoform data----
source("/Users/dogya203/Desktop/Komorowski_Project/src/cluster_rules.R")
cluster_rules(rmcfs_selected_cluster_aa %>% rownames_to_column(.,'SAMPLE') %>% left_join(z,by='SAMPLE') %>% column_to_rownames(.,'SAMPLE')
, recal_significant )




## ----VisuNet----
vis_out <- visunet(rules)


## ----VisuArc---
visuArc(vis_out, decision= 'Adrenergic',feature='6.49X49')




# Testing visualization
library(ggplot2)
library(dplyr)
library(patchwork)

# === CONFUSION MATRIX PLOTTING FUNCTION ===
plot_confusion_matrix <- function(true, pred, title) {
  levels_set <- c("Adrenergic", "Octopamine")
  
  cm <- as.data.frame(table(
    Predicted = factor(pred, levels = levels_set),
    Actual    = factor(true, levels = levels_set)
  ))
  
  cm <- cm %>%
    group_by(Actual) %>%
    mutate(Pct = round(Freq / sum(Freq) * 100, 1)) %>%
    ungroup()
  
  accuracy <- round(sum(true == pred) / length(true) * 100, 1)
  n_total  <- length(true)
  
  ggplot(cm, aes(x = Actual, y = Predicted, fill = Pct)) +
    geom_tile() +
    geom_text(aes(label = paste0(Freq, "\n", Pct, "%")),
              size = 5.5, fontface = "bold",
              color = ifelse(cm$Pct > 55, "white", "black")) +
    scale_fill_gradientn(
      colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      limits = c(0, 100),
      name   = "Recall (%)"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      title    = title,
      subtitle = paste0("Accuracy = ", accuracy, "%  |  n = ", n_total),
      x        = "True Class",
      y        = "Predicted Class"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle    = element_text(hjust = 0.5, face = "plain", size = 11, 
                                      color = "black", margin = margin(b = 8)),
      axis.text        = element_text(face = "bold", size = 12, color = "black"),
      axis.title.x     = element_text(face = "bold", size = 12, margin = margin(b = 6)),
      axis.title.y     = element_text(face = "bold", size = 12, margin = margin(r = 6)),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      legend.title     = element_text(face = "bold", size = 11),
      legend.text      = element_text(size = 10),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# === EXTRACT PREDICTIONS AND LABELS ===
true_test     <- decisions$out$currentClass
pred_test     <- decisions$out$predictedClass
test_rownames <- rownames(test_data)

# === IDENTIFY SAMPLE GROUPS ===
ncbi_idx    <- !grepl("^M|^CHIMERA", test_rownames)  # NCBI/Database samples
mutant_idx  <- grepl("^M", test_rownames) & !grepl("^CHIMERA", test_rownames)  # mutant samples (M only)
chimera_idx <- grepl("^CHIMERA", test_rownames)      # chimera samples

# === VERIFY SAMPLE COUNTS ===
cat("\n=== Sample Counts ===\n")
cat("NCBI samples:    ", sum(ncbi_idx), "\n")
cat("Mutant samples:  ", sum(mutant_idx), "\n")
cat("CHIMERA samples: ", sum(chimera_idx), "\n")
cat("Total test:      ", length(test_rownames), "\n")
cat("Sum check:       ", sum(ncbi_idx) + sum(mutant_idx) + sum(chimera_idx), 
    "(should equal total)\n\n")

# === CREATE CONFUSION MATRICES FOR EACH GROUP ===
p1 <- plot_confusion_matrix(
  true_test[ncbi_idx], 
  pred_test[ncbi_idx], 
  "NCBI Database Samples"
)

p2 <- plot_confusion_matrix(
  true_test[mutant_idx], 
  pred_test[mutant_idx], 
  "Mutant Samples"
)

p3 <- plot_confusion_matrix(
  true_test[chimera_idx], 
  pred_test[chimera_idx], 
  "CHIMERA Samples"
)

# === COMBINE PLOTS ===
combined_plot <- (p1 | p2 | p3) +
  plot_annotation(
    title = "Classification Performance on Test Set by Sample Type",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, 
                                margin = margin(b = 12))
    )
  )

# === DISPLAY ===
print(combined_plot)

# === SAVE HIGH-RESOLUTION OUTPUTS ===
ggsave("confusion_matrices.png", combined_plot, 
       width = 15, height = 5, dpi = 300, bg = "white")

ggsave("confusion_matrices.pdf", combined_plot, 
       width = 15, height = 5, bg = "white")

cat("\n=== Files Saved ===\n")
cat("✓ confusion_matrices.png (300 dpi raster)\n")
cat("✓ confusion_matricesl.pdf (vector format)\n")