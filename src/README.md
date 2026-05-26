# Doga Yalova Master's Thesis

In this project, we are studying determinants of two phylogenetically similar but pharmacologically different GPCRs using interpretable rule-based models.

The workflow is as follows:  

# insert image

## Code Organization
Each step of the workflow corresponds to specific R scripts and resources:

#### Data Retrieval Scripts
- The text files for both receptors with all accession IDs should be present.
- `data_preparation_adrenergic.R`
- `data_preparation_octopamine.R`
- `helper_functions.R`

#### Phylogenetic Tree of Data
- `phylogenetic_tree.R`

#### Multiple Sequence Alignment (MSA)
- Combined FASTA file (NCBI samples, mutants, chimera mutants, uncharacterized mutants) is submitted to
  Clustal-Omega MSA Tool: https://www.ebi.ac.uk/jdispatcher/msa/clustalo
  
#### Data Preparation, Train/Test Split
- `data_preparation.R`

#### Data Preparation Specific for the Model, MCFS, and ROSETTA

- Preparation of Distance Models from AF3 Predictions to Prepared Decision Table
- (See Appendix F)
- `cif2pdb.py`
- `pdb2tsv.py`
- `distance_pipeline.py` 
- `align_distances.py`

#### pLDDT Score Calculator and Visualization for Distance Model
- `plddt_score_calculator.py`
- `plddt_score_histogram.py`

- Preparation of Data Files
- **Sequence Model:** `sequence_model.R`
- **Physicochemical Model:** `physicochemical_model_prep.R`
- **Distance Model:** `distance_model_prep.R`

- Training and Testing 
- **Sequence Model:** `sequence_model.R`
- **Physicochemical Model:** `physicochemical_model.R`
- **Distance Model:** `distance_model.R`
- **Merged Model:** `merged_model.R`

- ROC Scores Overlay Plot
- `roc_overlay.R`

- Uncharacterized Mutant Plotting
- `uncharacterized_mutant_plot.py`

- Coloring the GPCR Residues According to Distances by Changing the B-Factor
- Column
- `changing_pdb_column_to_distance.py`

- Clustering Rules for Isoform Data
- Pipeline and code for sequence model was used for clustering the isoform data.

> **Note:** Code for sequence model `sequence_model.R` contains both preparation of the file and training/testing.
> **Note:** `wang.R`, `permutation_test.R`,`cluster_rules.R` are external script taken from the KomorowskiLab repository.
