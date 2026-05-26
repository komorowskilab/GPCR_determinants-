# Install necessary packages
if (!requireNamespace("rmcfs", quietly = TRUE)) install.packages("rmcfs")
if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("R.ROSETTA", quietly = TRUE)) {devtools::install_github("komorowskilab/R.ROSETTA")}
if (!requireNamespace("VisuNet", quietly = TRUE)) {devtools::install_github("komorowskilab/VisuNet")}
if (!requireNamespace("DT", quietly = TRUE)) install.packages("DT")


library(dplyr)
library(rJava)
library(rmcfs)
library(devtools)
library(R.ROSETTA)
library(csVisuNet)
library(DT)


###------------- LOAD THE DATA-------------------------

# Load the training
data_train <- read.csv("data_train_distance.csv", stringsAsFactors = FALSE, check.names = FALSE)
rownames(data_train) <- data_train[[1]]
data_train[[1]] <- NULL

# Load the testing data
data_test <- read.csv("data_test_distance.csv", stringsAsFactors = FALSE, check.names = FALSE)
rownames(data_test) <- data_test[[1]]
data_test[[1]] <- NULL

# Load the validation data
#data_validation <- read.csv("data_validation.csv", stringsAsFactors = FALSE, check.names = FALSE)
#rownames(data_validation) <- data_validation[[1]]
#data_validation[[1]] <- NULL


names(data_train)[names(data_train) != "Receptor"] <- paste0(names(data_train)[names(data_train) != "Receptor"], ".dist")
names(data_test)[names(data_test) != "Receptor"] <- paste0(names(data_test)[names(data_test) != "Receptor"], ".dist")
#names(data_validation)[names(data_validation) != "Receptor"] <- paste0(names(data_validation)[names(data_validation) != "Receptor"], ".dist")





# ----Feature Selection using RMCFS----
mcfs_train <- mcfs(Receptor ~ ., data_train, 
                                  projections = 20000, 
                                  projectionSize = 0.1, 
                                  splits = 5, 
                                  balance = 2,
                                  splitSetSize = 500, 
                                  cutoffPermutations = 30, 
                                  threadsNumber = 8,
                                  finalCV = TRUE,
                                  seed = 42)               # consistent results 


significant_mcfs <- mcfs_train$RI[1:mcfs_train$cutoff_value,]  



# ---- Save results ----
# Feature ranking
write.csv(mcfs_train$RI, "mcfs_feature_importance.csv", row.names = FALSE)

# Feature interdependencies
write.csv(mcfs_train$ID, "mcfs_feature_interdependencies.csv", row.names = FALSE)

# Cross-validation results
write.csv(mcfs_train$cv_accuracy, "mcfs_cv_accuracy.csv", row.names = FALSE)

# Significant features
write.csv(significant_mcfs, "mcfs_significant_features.csv", row.names = FALSE)


# ---- Plot distances ----
png(file.path("distances.png"), 
    width = 1200, height = 1000, res = 150)
plot(mcfs_train, type = "distances")
dev.off()
cat("✓ Saved distance plot\n")



# ---- Prepare filtered dataset for R.ROSETTA ----
selected_features <- as.character(significant_mcfs$attribute) # extract significant feature names from 'significant_result'
selected_features <- c(selected_features, "Receptor")      # add "Receptor" to the list of selected features
mcfs_selected <- data_train[, selected_features]
mcfs_selected$Receptor <- factor(mcfs_selected$Receptor, 
                                                   levels = c("Adrenergic", "Octopamine"))  



# ---- For full length model, is 22 is not put for NA values, you shuold fill before discretization------#
# Replace NAs with column means for each dataset
# For training data
#for(i in 1:(ncol(mcfs_selected) - 1)) {
#  mcfs_selected[is.na(mcfs_selected[,i]), i] <- 
#    mean(mcfs_selected[,i], na.rm = TRUE)
#}
#
# For test data
#for(i in 1:(ncol(data_test) - 1)) {
#  data_test[is.na(data_test[,i]), i] <- 
#    mean(data_test[,i], na.rm = TRUE)
#}

# For validation data
#for(i in 1:(ncol(data_validation) - 1)) {
#  data_validation[is.na(data_validation[,i]), i] <- 
#    mean(data_validation[,i], na.rm = TRUE)
#}

##------------------------------------------------------------------------###

# Start discretization
# This discretized data frame will be used as input for R.ROSETTA
source("wang.R")
trainros<- cdd_discretize(mcfs_selected, m = 50, alpha = 0.1, lambda = 0.05)
mcfs_selected_discretized <- trainros$discretized



# Reverse columns, then move "Receptor" to the end
#cols <- rev(names(mcfs_selected_discretized))
#cols_reordered <- c(cols[cols != "Receptor"], "Receptor")
#mcfs_selected_discretized_reordered <- mcfs_selected_discretized[, cols_reordered]


# R.ROSETTA
out <- rosetta(mcfs_selected_discretized, 
                 roc =TRUE, 
                 underSample =TRUE,
                 clroc = "Adrenergic",  
                 discrete =TRUE,
                 reducer = "Johnson",
                 JohnsonParam = list(Modulo=TRUE, BRT=TRUE, BRTprec=0.99, Precompute=TRUE, Approximate=TRUE, Fraction=0.99),)

rbm <- out$main
rbm_significant <- rbm[rbm$pValue < 0.05, ]



# Write to a plain text file (tab-separated)
write.table(
  rbm_significant,
  file = "rules_model_significant.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



# Count the number of rules for each decision class (Adrenergic and Octopamine)
rule_count <- table(rbm$decision)

# Model evaluation
quality_model <- out$quality

# ROC curve
png("roc_curve.png", width = 800, height = 600, res = 150)
plotMeanROC(out)
dev.off()


# Apply same discretization cuts to the test set
apply_cdd_cuts <- function(new_data, cutpoints_list) {
  result <- new_data
  
  for (gene in names(cutpoints_list)) {
    
    if (!gene %in% names(new_data)) {
      warning(paste(gene, "not found in new data — skipping."))
      next
    }
    
    cuts <- cutpoints_list[[gene]]
    
    if (is.null(cuts)) {
      # Non-informative gene: was assigned state 1 for all training samples
      result[[gene]] <- factor(rep(1, nrow(new_data)))
      
    } else if (length(cuts) == 1) {
      # 2-state variable
      result[[gene]] <- cut(new_data[[gene]], 
                            breaks = c(-Inf, cuts, Inf), 
                            labels = c(1, 2))
      
    } else if (length(cuts) == 2) {
      # 3-state variable
      result[[gene]] <- cut(new_data[[gene]], 
                            breaks = c(-Inf, cuts, Inf), 
                            labels = c(1, 2, 3))
    }
  }
  
  return(result)
}

# --- Discretization Usage ---
test_discretized <- apply_cdd_cuts(data_test, trainros$cutpoints)
validation_discretized <-apply_cdd_cuts(data_validation, trainros$cutpoints)



# Testing on External Data
# Store decision
data_test_decision <- test_discretized$Receptor
# Remove decision from table
test_discretized <- test_discretized[,-length(test_discretized)]
# Test
decisions <- predictClass(test_discretized, rbm_significant, discrete = TRUE, normalize = TRUE, normalizeMethod="rss", validate = TRUE, defClass = data_test_decision)
table(Actual = decisions$out$currentClass, Predicted = decisions$out$predictedClass)
table(decisions$out$currentClass == decisions$out$predictedClass)


# To see the support set
recal_data <- recalculateRules(mcfs_selected_discretized, rbm_significant, discrete=TRUE)
View(recal_data[str_detect((str_split(recal_data$supportSetLHS,',')),'XP_029642552'),])


# Reclassifying the training data
data_train_decision <- mcfs_selected_discretized$Receptor
mcfs_selected_discretized <- mcfs_selected_discretized[,-length(mcfs_selected_discretized)]
decisions_train <- predictClass(mcfs_selected_discretized, rbm_significant, discrete = TRUE, normalize = TRUE, normalizeMethod="rss", validate = TRUE, defClass = data_train_decision)
table(Actual = decisions_train$out$currentClass, Predicted = decisions_train$out$predictedClass)
table(decisions_train$out$currentClass == decisions_train$out$predictedClass)


# Test on validation (uncharacterized mutants)
decisions_validation <- predictClass(validation_discretized, rbm_significant, discrete=TRUE, normalize = TRUE, normalizeMethod="rss", validate = FALSE)


# VisuNet
vis_out <- visunet(rbm_significant)

# VisuArc
visuArc(vis_out, decision= 'Octopamine', feature='3.33X33.dist=1')



## Plots for testing
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
  
  # Row-wise % (% of actual class = recall perspective)
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
test_rownames <- rownames(data_test)

# === IDENTIFY SAMPLE GROUPS ===
ncbi_idx    <- !grepl("^M|^CHIMERA", test_rownames)  # NCBI/Database samples
mutant_idx  <- grepl("^M", test_rownames) & !grepl("^CHIMERA", test_rownames)  # Mutant samples (M only)
chimera_idx <- grepl("^CHIMERA", test_rownames)      # Chimera samples

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