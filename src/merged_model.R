# We are going to merge the decision tables of physichochemical and distance models
# It is more practical if the model codes were already run and the datasets are 
# in the environment of the R project

library(dplyr)
library(rJava)
library(rmcfs)
library(devtools)
library(R.ROSETTA)
library(VisuNet)
library(DT)
library(ggplot2)
library(dplyr)
library(patchwork)

# Reload, since we deleted the Receptor column in the last testing code
data_test_physicochemical <- read.table("test_data_physichochemical.csv", 
                           header = TRUE, sep = ",", colClasses = "character", row.names = 1,  check.names = FALSE)

# Use row names of the first table as reference order
# To get test_discretized_distance dataset, run the distance model code up to
# discretization step
row_order <- rownames(data_test_physicochemical)
test_discretized_distance <- test_discretized_distance[row_order, ]


merged_table <- cbind(
  data_test_physicochemical[, !names(data_test_physicochemical) %in% "Receptor", drop = FALSE],
  test_discretized_distance[, !names(test_discretized_distance) %in% "Receptor", drop = FALSE],
  Receptor = data_test_physicochemical$Receptor
)

#merged_table_validation <- cbind(
#  validation_data,
#  validation_discretized
#)

rules <- list(rbm_physichochemical_significant, rbm_distance_significant)

#merged_rules <- mergeRBMs(rules, fun = "mean", defClass = merged_table$Receptor)
merged_rules <- dplyr::bind_rows(rbm_physichochemical_significant,rbm_distance_significant)
head(merged_rules)  
merged_rules_significant <- merged_rules[merged_rules$pValue < 0.05, ]

# Test on the external data
# Store decision
merged_table_decision <- merged_table$Receptor
# Remove decision from table
merged_table <- merged_table[,-length(merged_table)]
# Test
decisions_merged <- predictClass(merged_table, merged_rules_significant, discrete = TRUE, normalize = TRUE, normalizeMethod="rss", validate = TRUE, defClass = merged_table_decision)
table(Actual = decisions_merged$out$currentClass, Predicted = decisions_merged$out$predictedClass)
table(decisions_merged$out$currentClass == decisions_merged$out$predictedClass)

# Test on validation (uncharacterized mutant) data
decisions_merged_validation <- predictClass(merged_table_validation, merged_rules_significant, discrete = TRUE, normalize = TRUE, normalizeMethod="rss", validate = FALSE)

## ----VisuNet----
vis_out <- visunet(merged_rules)

## ----VisuArc----
visuArc(vis_out, decision= 'Octopamine',feature='6.61X61.dist=1')




##-------------Check Model Agreement---------------------
get_positions <- function(rules) {
  positions <- sapply(rules$features, function(x) {
    features <- unlist(strsplit(x, ","))
    
    pos <- sapply(features, function(feat) {
      pos_match <- regmatches(feat, regexpr("\\d+\\.\\d+X?\\d+", feat))
      if(length(pos_match) > 0) pos_match else feat
    })
    
    return(pos)
  }, USE.NAMES = FALSE)
  
  unique(unlist(positions))
}

pos_sequence <- get_positions(rules)
pos_physicochemical <- get_positions(rbm_physichochemical_significant)
pos_distance <- get_positions(rbm_distance_significant)

cat("\n=== Positions per model ===\n")
cat("Model sequence positions:\n"); print(head(sort(pos_sequence)))
cat("Model physicochemical positions:\n"); print(head(sort(pos_physicochemical)))
cat("Model distance positions:\n"); print(head(sort(pos_distance)))

# Positions in all three models
in_all_three <- Reduce(intersect, list(pos_sequence, pos_physicochemical, pos_distance))
cat("\n=== Positions in all 3 models ===\n")
print(sort(in_all_three))

# Positions in any two models
in_sequence_and_physicochemical <- intersect(pos_sequence, pos_physicochemical)
in_sequence_and_distance <- intersect(pos_sequence, pos_distance)
in_physicochemical_and_distance <- intersect(pos_physicochemical, pos_distance)

# All combinations
cat("\n=== In sequence and physicochemical only (not distance) ===\n")
print(sort(setdiff(in_sequence_and_physicochemical, pos_distance)))

cat("\n=== In sequence and distance only (not physicochemical) ===\n")
print(sort(setdiff(in_sequence_and_distance, pos_physicochemical)))

cat("\n=== In physicochemical and distance only (not sequence) ===\n")
print(sort(setdiff(in_physicochemical_and_distance, pos_sequence)))

cat("\n=== Summary ===\n")
cat("Total unique positions - sequence:", length(pos_sequence), " | physicochemical:", length(pos_physicochemical), " | distance:", length(pos_distance), "\n")
cat("Overlap sequence & physicochemical:", length(in_sequence_and_physicochemical), "\n")
cat("Overlap sequence & distance:", length(in_sequence_and_distance), "\n")
cat("Overlap physicochemical & distance:", length(in_physicochemical_and_distance), "\n")
cat("All three:", length(in_all_three), "\n")



##----------------Plot---------------------------------------
# Testing visualization

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
true_test     <- decisions_merged$out$currentClass
pred_test     <- decisions_merged$out$predictedClass
test_rownames <- rownames(merged_table)

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