# Goal: perform conditional logistic regression. 
# Predictor variable(s): trehalose variants. 
# Response variable: severe C. difficile infection outcome. 
# Conditioned on stratum, which is based on severe outcome risk score. 

# Process outline: 
# Generate presence/absence tables for each trehalose variant. 
# If table has any zeroes remove that varaint from consideration. 
# For each remaining variant perform a single, unadjusted clogit. 
# Perform multiple test correction (fdr)
# Examine statistics for each model. 

#' run_clogit_models
#' This function runs clogit models on all trehalose variants input and saves 
#'   the model results in a table as both pdf and a .tsv file. 
#' @param metadata_path Character. Path to model inputs, strata, etc...
#'
#' @noRd
run_clogit_models <- function(metadata_path){
  metadata <- read_tsv(metadata_path, 
                       col_names = TRUE)
  
  # Subset to only those samples that had WGS performed, were matched into a 
  #   stratum, and the stratum is complete.
  # Subset metadata to only IDs, predictor variables, response variable, 
  #   stratum, ribotype. 
  
  model_input <- 
    metadata %>% filter(WGS_performed == 1, Stratum_complete == 1)
  col_to_drop <- c("treA2", 
                   "ptsT", 
                   "treX", 
                   "treR2", 
                   "Cys171Ser",
                   "age", 
                   "gender..M.0.F.1.", 
                   "METS", 
                   "concurrentabx",
                   "Lowest_SBP", 
                   "highcreat", 
                   "highbili", 
                   "WBC", 
                   "Duplicated_Patient", 
                   "Missing_Model_Data", 
                   "Risk_Score",           
                   "Duplicated_Patient_No_Missing_Info", 
                   "Unmatched", 
                   "WGS_performed", 
                   "Stratum_complete")
  model_input <- 
    model_input[ , !(colnames(model_input) %in% col_to_drop), drop = FALSE]
  model_input <- 
    model_input %>%
    select(-ID, everything()) %>% 
    select(C171S_L172I_or_insertion, everything())
  not_variants <- c("Severe_Outcome", "Ribotype", "Stratum", "ID")
  num_variants <- sum(!colnames(model_input) %in% not_variants)
  good_variant <- NULL
  for (i in 1:num_variants) {
    temp_table <- 
      model_input %>% 
      select(Severe_Outcome, colnames(model_input)[i]) %>%
      table
    # If any of the squares in these tables are 0 remove them from the analysis 
    # because they will not be interpretable/informative.
    if (sum(temp_table == 0) == 0) {
      if (ncol(temp_table) == 2) {
        if (nrow(temp_table) == 2) {
          good_variant <- append(good_variant, i)
        }
      }
    }
  }
  # Now all tables have non-zero entries in each cell and are 2x2
  num_good_loci <- length(good_variant)
  unadjusted_models <- rep(list(NULL), num_good_loci)
  counter <- 1
  for (i in good_variant) {
    unadjusted_models[[counter]] <- 
      clogit(formula = Severe_Outcome ~ 
               model_input[ , i, drop = TRUE] + 
               strata(Stratum), 
             data = model_input)
    counter <- counter + 1
  }
  # Information stored in clogit model: 
  # coef = logs odd ratio
  # exp(coef) = odds ratio
  # 95%CI = exp(confint(model))
  model_results <- matrix(NA, ncol = 5, nrow = num_good_loci)
  colnames(model_results) <- 
    c("Variant", "OR", "95% CI (lower)", "95% CI (upper)", "P-value")
  for (m in 1:num_good_loci) {
    model_results[m, 1] <- colnames(model_input)[good_variant[m]]
    model_results[m, 2] <-
      format(round(summary(unadjusted_models[[m]])[]$coefficients[2], 2), 
             nsmall = 2)
    model_results[m, 3] <- 
      format(round(summary(unadjusted_models[[m]])[]$conf.int[3], 2), 
             nsmall = 2)
    model_results[m, 4] <- 
      format(round(summary(unadjusted_models[[m]])[]$conf.int[4], 2),
             nsmall = 2)
    model_results[m, 5] <- 
      format(round(coef(summary(unadjusted_models[[m]]))[ ,'Pr(>|z|)'], 2), 
             nsmall = 2)
  }
  
  # Apply multiple test correction to P-values. 
  model_results[ , 5] <- 
    format(round(p.adjust(model_results[ , 5], method = "fdr"), 2), nsmall = 2)

  # Save table 
  my_theme <- 
    ttheme_minimal(base_size = 11, 
                   padding = unit(c(6, 1), "mm"),
                   core = list(fg_params = list(hjust = 1, x = 0.9), 
                              bg_params = list(fill = c("grey95", "white"))),
                   rowhead = list(fg_params = list(hjust = 0, x = 0)), 
                   colhead = list(fg_params = list()))
  table_plot <- tableGrob(model_results, theme = my_theme)
  pdf(paste0("../figures/", 
             Sys.Date(), 
             "_clogit_all_trehalose_variants_table.pdf"))
  grid.draw(table_plot)
  dev.off()
  
  # Subset to just the variants of interest
  var_of_interest <- c("C171S_L172I_or_insertion",
                       "Leu172Ile", 
                       "four_gene_insertion")
  smaller_model <- model_results[model_results[ , 1] %in% var_of_interest, ]
  table_plot <- tableGrob(smaller_model, theme = my_theme)
  pdf(paste0("../figures/",
             Sys.Date(), 
             "_clogit_only_trehalose_variants_of_interest_table.pdf"))
  grid.draw(table_plot)
  dev.off()
  
  write_tsv(as_tibble(model_results), 
            paste0("../data/outputs/", 
                   Sys.Date(), 
                   "_clogit_all_trehalose_variants.tsv"))
  
  write_tsv(as_tibble(model_results), 
            paste0("../data/outputs/Supplementary_table_2.tsv"))
  
  write_tsv(as_tibble(smaller_model), 
            paste0("../data/outputs/", 
                   Sys.Date(), 
                   "_clogit_only_trehalose_variants_of_interest.tsv"))
} # end run_clogit_models()