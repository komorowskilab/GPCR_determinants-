# helper_functions.R
# Helper functions for GenPept record parsing

parse_gp_record <- function(rec) {
  lines <- strsplit(rec, "\n")[[1]]
  if (length(lines) < 5) return(NULL)  # skip empty/malformed
  
  # Accession
  acc <- grep("^ACCESSION", lines, value=TRUE)
  acc <- if (length(acc) > 0) sub("ACCESSION +([A-Z0-9_.-]+).*", "\\1", acc) else NA
  
  # Skip if no accession found
  if (is.na(acc)) return(NULL)
  
  # Definition/Title
  def_line <- grep("^DEFINITION", lines, value=TRUE)
  title <- if (length(def_line)>0) sub("DEFINITION +", "", def_line) else NA
  
  # Organism
  org_line <- grep("^  ORGANISM", lines, value=TRUE)
  organism <- if (length(org_line)>0) sub("^  ORGANISM +", "", org_line) else NA
  
  # Genus / species
  genus <- species <- NA
  if (!is.na(organism)) {
    parts <- strsplit(organism, " +")[[1]]
    if (length(parts) >= 1) genus <- parts[1]
    if (length(parts) >= 2) species <- parts[2]
  }
  
  # Extract taxonomic lineage (phylo_info)
  phylo_info <- list()
  org_idx <- grep("^  ORGANISM", lines)
  if (length(org_idx) > 0) {
    # Find the end of the ORGANISM section (next major section or empty line)
    next_section <- grep("^[A-Z]", lines[(org_idx+1):length(lines)])[1]
    end_idx <- if (!is.na(next_section)) org_idx + next_section - 1 else length(lines)
    
    # Get all lines after ORGANISM until next section
    tax_lines <- lines[(org_idx+1):(end_idx-1)]
    # Remove leading spaces and combine
    tax_text <- paste(trimws(tax_lines), collapse=" ")
    # Split by semicolon and clean up
    tax_levels <- strsplit(tax_text, ";")[[1]]
    tax_levels <- trimws(tax_levels)
    tax_levels <- tax_levels[tax_levels != ""]
    
    # Store as phylo_info_1, phylo_info_2, etc.
    if (length(tax_levels) > 0) {
      for (i in seq_along(tax_levels)) {
        phylo_info[[paste0("phylo_info_", i)]] <- tax_levels[i]
      }
    }
  }
  
  # Vertebrate/Invertebrate status
  vertebrate_keywords <- c("Vertebrata", "Chordata") # expand as needed
  
  phylo_class <- "Invertebrate"
  if (length(phylo_info) > 0) {
    all_phylo <- paste(unlist(phylo_info), collapse=" ")
    if (any(sapply(vertebrate_keywords, function(x) grepl(x, all_phylo, ignore.case=TRUE)))) {
      phylo_class <- "Vertebrate"
    }
  }
  
  # Taxon ID (from db_xref, it can be searched on https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)
  taxon_id <- NA
  taxon_line <- grep("db_xref=\"taxon:", lines, value=TRUE)
  if (length(taxon_line) > 0) {
    taxon_id <- sub('.*db_xref="taxon:([0-9]+)".*', "\\1", taxon_line[1])
  }
  
  # Sex information
  sex <- NA
  sex_line <- grep("/sex=", lines, value=TRUE)
  if (length(sex_line) > 0) {
    sex <- sub('.*\\/sex="?([^"]+)"?.*', "\\1", sex_line[1])
  }
  
  # Tissue type
  tissue <- NA
  tissue_line <- grep("/tissue_type=", lines, value=TRUE)
  if (length(tissue_line) > 0) {
    tissue <- sub('.*\\/tissue_type="?([^"]+)"?.*', "\\1", tissue_line[1])
  }
  
  # Sequence length
  locus_line <- grep("^LOCUS", lines, value=TRUE)
  length_aa <- ifelse(length(locus_line)>0,
                      as.numeric(sub(".* ([0-9]+) aa.*", "\\1", locus_line)), NA)
  
  # Patent
  has_patent <- ifelse(!is.na(title) && grepl("patent", title, ignore.case=TRUE), "yes", "no")
  
  # Isoform / receptor-like info 
  isoform_info <- "no"
  
  if (!is.na(title)) {
    # Clean title: remove anything starting from [, remove trailing period, trim spaces
    clean_title <- gsub("\\s*\\[.*$", "", title)
    clean_title <- sub("\\.$", "", clean_title)
    clean_title <- trimws(clean_title)
    
    # First priority: isoform
    if (grepl("isoform", clean_title, ignore.case=TRUE)) {
      match <- regexpr("isoform.*", clean_title, ignore.case=TRUE)
      isoform_info <- substring(clean_title, match)
      
      # Second priority: receptor-like
    } else if (grepl("receptor-like", clean_title, ignore.case=TRUE)) {
      match <- regexpr("receptor-like.*", clean_title, ignore.case=TRUE)
      isoform_info <- substring(clean_title, match)
      
      # Third priority: other keywords
    } else {
      other_keywords <- c("beta", "like", "variant")
      if (any(grepl(paste(other_keywords, collapse="|"), clean_title, ignore.case=TRUE))) {
        isoform_info <- clean_title
      }
    }
    
    # Final trim (just in case)
    isoform_info <- trimws(isoform_info)
  }
  
  # TM helices (per helix + overall completeness) 
  tm_presence <- sapply(tm_keywords, function(k) grepl(k, rec, ignore.case=TRUE))
  names(tm_presence) <- paste0("TM", 1:7)
  
  all_7_TM <- ifelse(all(tm_presence), "yes", "no")
  
  # Sequence
  origin_start <- grep("^ORIGIN", lines)
  sequence <- NA
  if (length(origin_start) > 0) {
    seq_lines <- lines[(origin_start+1):length(lines)]
    seq_lines <- gsub("[0-9 ]", "", seq_lines)
    sequence <- toupper(paste(seq_lines, collapse=""))
  }
  
  # Build the base row
  row <- data.frame(
    accession = ifelse(length(acc) == 0, NA, acc),
    title = ifelse(length(title) == 0, NA, title),
    organism = ifelse(length(organism) == 0, NA, organism),
    genus = ifelse(length(genus) == 0, NA, genus),
    species = ifelse(length(species) == 0, NA, species),
    taxon_id = taxon_id,
    sex = sex,
    tissue_type = tissue,
    length_aa = ifelse(length(length_aa) == 0, NA, length_aa),
    has_patent = ifelse(length(has_patent) == 0, NA, has_patent),
    isoform_info = ifelse(length(isoform_info) == 0, NA, isoform_info),
    all_7_TM_helices = all_7_TM,
    stringsAsFactors = FALSE
  )
  
  # Add TM columns
  for (i in 1:7) {
    row[[paste0("TM", i)]] <- tm_presence[i]
  }
  
  # Add phylo_info columns dynamically
  if (length(phylo_info) > 0) {
    for (col_name in names(phylo_info)) {
      row[[col_name]] <- phylo_info[[col_name]]
    }
  }
  
  # Store phylo_info as a collapsed string for consistency
  row$phylo_lineage <- if (length(phylo_info) > 0) {
    paste(unlist(phylo_info), collapse="; ")
  } else {
    NA
  }
  
  return(list(row = row, sequence = sequence))
}

