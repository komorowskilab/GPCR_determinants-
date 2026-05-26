install.packages("ape")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("treeio", "ggtree"))

# Load libraries
library(ape)
library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)
library(phangorn)
library(viridis)
library(RColorBrewer)
library(randomcoloR)

# Read the Newick tree
tree <- read.tree("tree.tree")

# Read metadata tables
meta1 <- read.csv("adrenergic_metadata.csv")
meta2 <- read.csv("octopamine_metadata.csv")
meta1$source <- "adrenergic"
meta2$source <- "octopamine"
meta <- bind_rows(meta1, meta2)


# Join metadata with tree, match on 'label'
tree_data <- ggtree(tree) %<+% meta

# Plot with just organism names and source
p <- tree_data +
  geom_tiplab(aes(label = organism), size = 3) +
  geom_tippoint(aes(shape = source), size = 3, stroke = 0.5) +
  ggtitle("Annotated Phylogenetic Tree with Dataset Source") +
  theme_tree2() +
  theme(plot.title = element_text(hjust = 0.5))
p


# Save the normal plot
ggsave("phylogenetic_tree_large.pdf", 
       plot = p, 
       width = 25,      
       height = 50,    
       units = "in",
       limitsize = FALSE)  # Allows very large sizes


# Save the circular plot
p_circular <- p + layout_circular()
# Increase label size for circular tree
p_circular_detailed <- tree_data +
  geom_tiplab(aes(label = isoform_info, color = isoform_info), 
              size = 4,          # Larger text
              offset = 0.01) +   # Adjust spacing from tips
  geom_tippoint(aes(shape = source, fill = isoform_info), size = 4, stroke = 0.5) +
  scale_color_manual(values = isoform_colors) +
  scale_fill_manual(values = isoform_colors) +
  ggtitle("Annotated Phylogenetic Tree with Dataset Source") +
  theme_tree2() +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  guides(
    fill = "none",
    color = guide_legend(override.aes = list(shape = NA))
  ) +
  layout_circular()
p_circular_detailed
ggsave("tree_circular_ultra_detailed.pdf", 
       plot = p_circular_detailed, 
       width = 40, 
       height = 40, 
       units = "in",
       limitsize = FALSE)


