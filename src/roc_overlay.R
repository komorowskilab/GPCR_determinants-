# For overlaying the ROC curves

plotMeanROC_combined <- function(out_list, labels, colors) {
  
  par(mar = c(5, 5, 4, 2), bg = "white", family = "serif")
  
  # set up empty plot
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity (False Positive Rate)",
       ylab = "Sensitivity (True Positive Rate)",
       axes = FALSE, cex.lab = 1.3, font.lab = 2,
       main = "Receiver Operating Characteristic Curves",
       cex.main = 1.4, font.main = 2)
  
  # clean axes
  axis(side = 1, at = seq(0, 1, 0.1), las = 1, cex.axis = 1.1, font = 1)
  axis(side = 2, at = seq(0, 1, 0.1), las = 2, cex.axis = 1.1, font = 1)
  box(lwd = 1.5)
  
  # reference diagonal
  abline(0, 1, lty = 2, col = "gray50", lwd = 1.5)
  
  # grid lines
  grid(nx = 10, ny = 10, col = "gray90", lty = 1, lwd = 0.5)
  
  # plot each model
  auc_values <- c()
  for (i in seq_along(out_list)) {
    out <- out_list[[i]]
    ROCstats <- out$ROCstats
    if (is.null(ROCstats)) ROCstats <- out$ROC.stats
    
    OMSpec <- rowMeans(unstack(ROCstats, form = OneMinusSpecificity ~ CVNumber))
    Sens   <- rowMeans(unstack(ROCstats, form = Sensitivity ~ CVNumber))
    auc_val <- round(out$quality$ROC.AUC.MEAN, digits = 3)
    auc_values[i] <- auc_val
    
    lines(OMSpec, Sens, lwd = 2.5, col = colors[i])
  }
  
  # legend with AUC values
  legend("bottomright",
         legend = paste0(labels, "  (AUC = ", auc_values, ")"),
         col = colors,
         lwd = 2.5,
         bty = "n",           # no legend box
         cex = 1.1,
         seg.len = 1.5,
         y.intersp = 1.3)
}

# Adjust the rosetta outputs
plotMeanROC_combined(
  out_list = list(out1, out2, out3),
  labels   = c("Sequence", "Physicochemical ", "Distance "),
  colors = c("#0072B2",  # blue
             "#E69F00",  # orange/yellow
             "#CC79A7")  # pink/purple
)

tiff("ROC_curves.tiff", width = 3500, height = 3500, res = 300)
plotMeanROC_combined(
  out_list = list(out1, out2, out3),
  labels   = c("Sequence", "Phylo", "Pharma"),
  colors   = c("#0072B2", "#E69F00", "#CC79A7") #color blind friendly
)
dev.off()
