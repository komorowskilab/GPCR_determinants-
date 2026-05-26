# Data preparation code for alpha 2 adrenergic receptor
# (both fasta and metadata retrieval from GenPept record)


# ---- Packages, Libraries & Loading External Functions ----

# Install if not already installed
install.packages("rentrez")
install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("scales")
install.packages("rlang")


library(rentrez)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

# Load helper functions
source("helper_functions.R")


# ---- Parameters ----

db <- "protein"
chunk_size <- 50

# TM helix keywords (needed by parse_gp_record function)
tm_keywords <- paste0("TM helix ", 1:7)


# ---- Full Accession List ----
accessions_adrenergic <- readLines("accessions_adrenergic.txt")

# Check the number of samples 
cat("Loaded", length(accessions_adrenergic), "accessions\n")


# ---- Initialize Containers and Retrieve Data ----

all_rows <- list()
all_fasta <- list()

for (i in seq(1, length(accessions_adrenergic), by=chunk_size)) {
  this_chunk <- accessions_adrenergic[i:min(i + chunk_size - 1, length(accessions_adrenergic))]
  cat("Fetching accessions", i, "to", min(i + chunk_size - 1, length(accessions_adrenergic)), "\n")
  
  gp_chunk <- tryCatch(
    entrez_fetch(db = db, id = this_chunk, rettype = "gp", retmode = "text"),
    error = function(e) { warning(e); return(NULL) }
  )
  
  if (is.null(gp_chunk)) next
  
  records <- strsplit(gp_chunk, "//")[[1]]
  # Filter out empty or whitespace-only records
  records <- records[nchar(trimws(records)) > 0]
  
  for (rec in records) {
    parsed <- parse_gp_record(rec)
    if (is.null(parsed) || is.null(parsed$row)) next
    
    # Add metadata row
    all_rows[[length(all_rows) + 1]] <- parsed$row
    
    # Add sequence if present
    if (!is.null(parsed$sequence) && !is.na(parsed$sequence) && nchar(parsed$sequence) > 0) {
      fasta_entry <- paste0(">", parsed$row$accession, " ", parsed$row$organism, "\n", parsed$sequence)
      all_fasta[[length(all_fasta) + 1]] <- fasta_entry
    }
  }
  
  Sys.sleep(0.3)  
}



# ---- Save .csv and fasta ----

# Combine all metadata rows
adrenergic_metadata <- bind_rows(all_rows)

# Reset row names to be clean
rownames(adrenergic_metadata) <- NULL

# Export .csv metadata
write.csv(adrenergic_metadata, "adrenergic_metadata.csv", row.names = FALSE)

# Combine all FASTA entries into one string
fasta_lines_adrenergic <- unlist(all_fasta)

# Export fasta file
writeLines(fasta_lines_adrenergic, "adrenergic.fasta")



# ---- Results Visualization ----

# Summary table: number of species per genus 
genus_summary_adrenergic <- adrenergic_metadata %>%
  filter(!is.na(genus), !is.na(species)) %>%
  group_by(genus) %>%
  summarise(
    n_species = n_distinct(species),
    n_entries = n(),
    taxon_id = first(taxon_id),
    # Automatically add all phylo_info columns that exist in the data
    across(starts_with("phylo_info_"), ~first(.x), .names = "{.col}"),
    phylo_lineage = first(phylo_lineage)
  ) %>%
  arrange(desc(n_entries))

print(genus_summary_adrenergic)

# Optionally save as CSV
write.csv(genus_summary_adrenergic, "genus_summary_adrenergic.csv", row.names = FALSE)


# Count the isoforms
adrenergic_isoform_count <- adrenergic_metadata %>%
  dplyr::count(isoform_info, name = "isoform_count")


# Optionally save as CSV
write.csv(adrenergic_isoform_count, "adrenergic_isoform_count.csv", row.names = FALSE)


# TM helices completeness
ggplot(adrenergic_metadata, aes(x=all_7_TM_helices)) +
  geom_bar(color="black") +
  geom_text(stat="count",
            aes(label=..count..),
            position=position_stack(vjust=0.5),
            color="black", size=5) +
  scale_fill_manual(values=c("grey60", "firebrick"),) +
  labs(title="TM Helices Completeness (Adrenergic Dataset)",
       x="Presence of All 7 TM helices",
       y="Number of sequences") +
  theme_minimal(base_size = 14)



# TM completeness per helix
tm_counts_adrenergic <- colSums(adrenergic_metadata[, paste0("TM", 1:7)], na.rm = TRUE)
tm_df_adrenergic <- data.frame(
  TM = paste0("TM", 1:7),
  Count = tm_counts_adrenergic
)

ggplot(tm_df_adrenergic, aes(x=TM, y=Count)) +
  geom_bar(stat="identity", fill="grey60") +
  geom_text(aes(label=Count), vjust=-0.3, size=5, fontface="bold") +
  labs(title="Individual TM Helix Completeness (Adrenergic Dataset)",
       x="Transmembrane helix", y="Number of sequences containing helix") +
  theme_minimal(base_size = 14)


# See which entries are missing TM1,TM2,TM3....
absent_TM7_adrenergic <- adrenergic_metadata %>%
  filter(TM1 == FALSE)
# Preview
head(absent_TM7_adrenergic[, c("accession", "organism", "title")])


# Receptor-like isoforms
ggplot(adrenergic_metadata, aes(x=isoform_info)) +
  geom_bar() +
  scale_fill_manual(values=c("grey60", "firebrick"),) +
  geom_text(stat="count", aes(label=..count..), 
            position=position_stack(vjust=0.5), size=5, fontface="bold") +
  labs(title="Isoform Distribution (Adrenergic Dataset) ",
       x="Isoform info", y="Number of sequences") +
  theme_minimal(base_size = 14)

# Histograms of length_aa to spot unusually short/long sequences
#Basic Length Plot
ggplot(adrenergic_metadata, aes(x=length_aa)) +
  geom_histogram(binwidth=10, color="black", position="stack") +
  scale_fill_manual(values=c("grey60", "firebrick"),) +
  labs(title="Protein Sequence Lengths of Adrenergic Dataset",
       x="Sequence length (aa)", y="Count") +
  theme_minimal(base_size = 14)

#Detailed Plot
# Calculate statistics
mean_length <- mean(adrenergic_metadata$length_aa)
median_length <- median(adrenergic_metadata$length_aa)
min_length <- min(adrenergic_metadata$length_aa)
max_length <- max(adrenergic_metadata$length_aa)

# Plot with all statistics
ggplot(adrenergic_metadata, aes(x=length_aa)) +
  geom_histogram(binwidth=10, color="black", fill="grey60") +
  
  # Mean line
  geom_vline(aes(xintercept=mean_length), 
             color="firebrick", linetype="dashed", size=1) +
  
  # Median line
  geom_vline(aes(xintercept=median_length), 
             color="blue", linetype="dashed", size=1) +
  
  # Min line
  geom_vline(aes(xintercept=min_length), 
             color="darkgreen", linetype="dotted", size=1) +
  
  # Max line
  geom_vline(aes(xintercept=max_length), 
             color="purple", linetype="dotted", size=1) +
  
  # Labels
  annotate("text", x=mean_length, y=Inf, 
           label=paste("Mean =", round(mean_length, 1)), 
           vjust=1.5, hjust=-0.1, color="firebrick", size=3.5) +
  annotate("text", x=median_length, y=Inf, 
           label=paste("Median =", round(median_length, 1)), 
           vjust=3, hjust=-0.1, color="blue", size=3.5) +
  annotate("text", x=min_length, y=Inf, 
           label=paste("Min =", min_length), 
           vjust=4.5, hjust=-0.1, color="darkgreen", size=3.5) +
  annotate("text", x=max_length, y=Inf, 
           label=paste("Max =", max_length), 
           vjust=6, hjust=1.1, color="purple", size=3.5) +
  
  labs(title="Protein Sequence Lengths of Adrenergic Dataset",
       x="Sequence length (aa)", y="Count") +
  theme_minimal(base_size=14)

# Get species names for min and max
min_species <- adrenergic_metadata$organism[which.min(adrenergic_metadata$length_aa)]
max_species <- adrenergic_metadata$organism[which.max(adrenergic_metadata$length_aa)]
print(min_species)
print(max_species)