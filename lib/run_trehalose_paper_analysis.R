# Run all functions necessary to generate trehalose paper figures & statstics. 

# Packages --------------------------------------------------------------------#
library(tidyverse)
library(ape)
library(survival)
library(seqinr)

# Functions -------------------------------------------------------------------# 
source("../lib/create_risk_scores_and_strata.R")
source("../lib/save_intermediate_data.R")
source("../lib/generate_descriptive_stats.R")
source("../lib/save_r_session_log.R")
source("../lib/plot_tree.R")
source("../lib/clogit.R")
source("../lib/variant_imputation.R")
source("../lib/save_four_gene_insertion_sequences.R")

# Generate results and figures for trehalose paper ----------------------------#
model_variables <- 
  read_tsv("../data/inputs/model_variables_plus_ribo_and_id.tsv")

# Create risk scores and strata based on KR logistic regression
create_risk_scores_and_strata(model_variables)

# Format additional input data necessary for analysis 
out_group <- "ERR232398"
keeper_path <- 
  "../data/inputs/key.csv"
tree_path <- "../data/inputs/input_tree.tree" 
snp_path <- "../data/inputs/SNP_matrix.tsv"
pan_path <- "../data/inputs/gene_presence_absence.Rtab"

model_path <- paste0("../data/outputs/severity_model.tsv")

save_data_for_tre_analysis(out_group, 
                           keeper_path,
                           tree_path, 
                           snp_path, 
                           pan_path,
                           model_path)

# Describe study cohort (ribotype, presence of trehalose variants, etc...)
metadata_path <- "../data/outputs/pre-analysis_trehalose_metadata.tsv"
generate_stats(metadata_path)

# Plot tree with trehalose variants on it 
updated_tree_path <- "../data/outputs/trehalose.tree"
plot_trehalose_tree(updated_tree_path, metadata_path) # FALSE to drop tip.labels

# Calculate conditional logistic regression and generate tables of results
run_clogit_models(metadata_path)

# Run imputation for trehalose variant presence on non-sequenced isolates
impute_variants(metadata_path, num_perm = 1000)

# Save sequence of four gene trehalose insertion
pan_seq_path <- "../data/inputs/pan_genome_reference.fa"
save_four_gene_insertion(pan_seq_path)

# Record R version and packages 
save_session_info()

# End of script ---------------------------------------------------------------#