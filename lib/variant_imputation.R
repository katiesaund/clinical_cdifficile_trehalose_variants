#' calculate_variant_ratio_by_ribo
#' For each ribotype in the WGS'd cohort count what percent of the isolates from 
#'   each ribotype has any of the three trehalose variants of interest (C171S, 
#'   L172I, or four gene insertion.) Return a tibble with the observed trehalose 
#'   utilization variant ratio for each ribotype. Save a table for supplement
#'   with information about ribotype-specific ratio for trehalose variant 
#'   presence. Calculate the ratio by ribotype for controls based only on data
#'   obtained from controls; similarly, only calculate the ratio by ribotype for
#'   cases based on data obtained from cases. This is essential to keep with the
#'   case-control study upon which this data is based.
#' @param variant_data Tibble. All of the relevant metadata. Rows correspond to
#'   individual isolates and columns contain, among other info, trehalose
#'   variants. 
#' @param case_or_ctrl String. Either "case" or "control."
#'
#' @return obs_tre_var_by_ribo. The fraction of samples that have the 
#'   C171S_L172I_or_insertion by ribotype. 
#' @noRd
calculate_variant_ratio_by_ribo <- function(variant_data, case_or_ctrl){
  is_case <- FALSE
  if (case_or_ctrl == "case") {
    is_case <- TRUE
  }
  obs_tre_var_by_ribo <- 
    variant_data %>%
    filter(!is.na(C171S_L172I_or_insertion)) %>%
    filter(Stratum_complete == TRUE) %>% 
    filter(Severe_Outcome == is_case) %>% 
    select(Ribotype, C171S_L172I_or_insertion) %>% 
    table() 
  obs_tre_var_by_ribo <- 
    cbind(obs_tre_var_by_ribo, rep(NA, nrow(obs_tre_var_by_ribo)) )
  colnames(obs_tre_var_by_ribo) <-
    c("tre_variant_absent", "tre_variant_present", "present_ratio")
  obs_tre_var_by_ribo <- as_tibble(obs_tre_var_by_ribo, rownames = "Ribotype")
  for (i in 1:nrow(obs_tre_var_by_ribo)) {
    temp_num_in_ribo <-
      sum(obs_tre_var_by_ribo$tre_variant_present[i] + 
            obs_tre_var_by_ribo$tre_variant_absent[i])
    obs_tre_var_by_ribo$present_ratio[i] <- 
      obs_tre_var_by_ribo$tre_variant_present[i] / temp_num_in_ribo
  }
  obs_tre_var_by_ribo <- obs_tre_var_by_ribo %>% select(Ribotype, present_ratio)
  
  save_version <- obs_tre_var_by_ribo
  save_version$present_ratio <- save_version$present_ratio * 100
  colnames(save_version)[2] <- 
    "% of isolates with trehalose utilization variant"
  write_tsv(save_version, 
            path = paste0("../data/outputs/", 
                          Sys.Date(), 
                          "_trehalose_prevalence_by_ribotype_",
                          case_or_ctrl, 
                          ".tsv"), 
            col_names = TRUE)
  return(obs_tre_var_by_ribo)
} # end calculate_variant_ratio_by_ribo()

#' drop_ribo_absent_from_matched_WGS
#' Remove those ribotypes from consideration which had no isolates that were 
#'   ever sequenced, because if there is no sequencing then the trehalose 
#'   utilization variant presence/absence cannot be calculated. If that's not 
#'   calculated then we cannot use those ribotypes for imputation of 
#'   presence/absence.
#'
#' @param variant_data Tibble. All of the relevant metadata. Rows correspond to
#'   individual isolates and columns contain, among other info, ribotypes.
#'
#' @return variant_data. Tibble. Same format as input variant_data, but the 
#'   ribotypes that are present in the non-sequenced cohort but absent in the 
#'   sequenced cohort are removed. 
#' @noRd
drop_ribo_absent_from_matched_WGS <- function(variant_data){
  ribo_present_in_WGS <- 
    variant_data %>% 
    filter(WGS_performed == 1, Stratum_complete == 1) %>% 
    select(Ribotype) %>% unique(.) %>% 
    unname() %>% 
    unlist()
  ribo_present_in_not_WGS <- 
    variant_data %>% 
    filter(WGS_performed == 0) %>% 
    select(Ribotype) %>% 
    unique(.) %>% 
    unname() %>% 
    unlist()
  ribo_absent_from_WGS <- setdiff(ribo_present_in_not_WGS, ribo_present_in_WGS)
  variant_data <- variant_data %>% filter(!(Ribotype %in% ribo_absent_from_WGS))
  return(variant_data)
} # end drop_ribo_absent_from_matched_WGS()

#' predict_num_with_var
#' Given the observed frequency with which trehalose utilization variants occur 
#'   by ribotype in the sequenced cohort and the number of not sequenced 
#'   isolates you have per ribotype, simply multiply these two numbers together 
#'   to get the total number isolates that should be assigned to "presence" 
#'   during the later imputation step. Only use control-dervied data for 
#'   controls; similarly, only used case-derived data for cases. This is
#'   essential to keep with the case-control study design upon which this data
#'   is based.
#'
#' @param variant_data Tibble. All of the relevant metadata. Rows correspond to
#'   individual isolates and columns contain, among other info, ribotypes.
#' @param tre_var_by_ribo Tibble. Each row is a different ribotype. Columns 
#'   include ribotype and the ratio of variant presence in that ribotype. 
#' @param case_or_ctrl String. Either "case" or "control."
#'
#' @return num_of_missing_data_by_ribotype. Tibble. The total number of isolates 
#'   that should get assigned to presence during imputation by ribotype. 
#' @noRd
predict_num_with_var <- function(variant_data, tre_var_by_ribo, case_or_ctrl){
  is_case <- FALSE
  if (case_or_ctrl == "case") {
    is_case <- TRUE
  }
  non_sequenced <- variant_data %>% 
    filter(WGS_performed == FALSE | Stratum_complete == FALSE) %>% 
    filter(Severe_Outcome == is_case)
  num_of_missing_data_by_ribotype <- 
    as.data.frame(non_sequenced %>% select(Ribotype) %>% table())
  colnames(num_of_missing_data_by_ribotype)[1] <- "Ribotype"
  num_of_missing_data_by_ribotype <- 
    left_join(num_of_missing_data_by_ribotype, tre_var_by_ribo, by = "Ribotype")
  num_of_missing_data_by_ribotype <- 
    num_of_missing_data_by_ribotype %>%
    mutate("assigned_to_1" = round(Freq * present_ratio, 0))
  return(num_of_missing_data_by_ribotype)
} # end predict_num_with_var()

#' generate_imputed_data
#' Generate a tibble with the imputed data. One row for each sample that
#'   wasn't sequenced. Each column, besides the ID, is a one run of the 
#'  imputation process. 
#'  Only generate imputed data on controls based on controls. 
#'  Similarly, only generate imputed data on cases based on cases. 
#' @param variant_data Tibble. Metadata. Isolates in rows. Features in columns.
#' @param n_perm Number. Number of times to run imputation step. 
#' @param n_per_ribo Tibble. Number of isolates to assign to trehalose 
#'   utilization variant presence per ribotype. 
#'
#' @return perm_imputed_data. Tibble with n_perm + 1 columns. 1 column is ID, 
#'   other columns correspond to imputations. Each row is an non-sequenced 
#'   isolate. 
#' @noRd
generate_imputed_data <- function(variant_data, n_perm, n_per_ribo, case_or_ctrl){
  is_case <- FALSE
  if (case_or_ctrl == "case") {
    is_case <- TRUE
  }
  
  non_sequenced <- variant_data %>% 
    filter(WGS_performed == FALSE | Stratum_complete == FALSE) %>% 
    filter(Severe_Outcome == is_case)
  
  n_per_ribo <- n_per_ribo %>% filter(!is.na(present_ratio)) 
  num_no_wgs_ribo <- length(unique(n_per_ribo$Ribotype))
  
  non_sequenced <- non_sequenced %>% filter(Ribotype %in% n_per_ribo$Ribotype)
  num_no_wgs_samples <- non_sequenced %>% nrow(.)
  
  
  perm_imputed_data <- matrix(0, nrow = num_no_wgs_samples, ncol = n_perm)
  row.names(perm_imputed_data) <- non_sequenced$ID
  set.seed(1)
  for (ribo in 1:num_no_wgs_ribo) {
    curr_ribo <- unique(n_per_ribo$Ribotype)[ribo]
    num_isolates <- 
      n_per_ribo %>% 
      filter(Ribotype == curr_ribo) %>% 
      select(assigned_to_1) %>% 
      unlist() %>%  
      unname()
    curr_ids <- 
      non_sequenced %>% 
      filter(Ribotype == curr_ribo) %>% 
      select(ID) %>% 
      unlist() %>% 
      unname()
    
    for (j in 1:n_perm) {
      ids_with_variants <- 
        sample(x = curr_ids, size = as.numeric(num_isolates), replace = FALSE)
      indices_with_variants <- 
        which(row.names(perm_imputed_data) %in% ids_with_variants) 
      perm_imputed_data[indices_with_variants, j] <- 1
    }
  }
  
  perm_imputed_data <- as_tibble(perm_imputed_data, rownames = "ID")
  return(perm_imputed_data)
} # end generate_imputed_data()

#' join_imputed_to_seq
#' This function joins the imputed variant data to the observed variant data. 
#'   Said another way, this function joins the non-matched data to the 
#'   matched data. It does this by replicating the matched variant data 
#'   n_perm times, then row-binding the imputed data to this replicated and
#'   matched data. At the end each column contains the matched & sequenced
#'   variant data and the imputed (non-matched) variant data. This format
#'   allows the analysis in the later functions to run tests on columns. 
#'
#' @param original_data Tibble. Metadata. Isolates in rows. The whole 1144 of 
#'   the original cohort. 
#' @param matched_seq_data Tibble. Metadata. Isolates in rows. Only data that 
#'   were in complete strata (matched) and sequenced. 
#' @param num_perm Number. Number of imputations performed.
#' @param all_data Tibble. Imputed data. Each rows is an isolate (names in ID 
#'   column). Other columns correspond to imputations. 
#'
#' @return var_data Tibble with imputed and observed variant data joined into 
#'   one tibble. 
#' @noRd
join_imputed_to_seq <- function(original_data, matched_seq_data, num_perm, imp_mat){
  temp_seq_var <- 
    matched_seq_data %>% 
    select(c(ID, C171S_L172I_or_insertion))
  temp_seq_var <- 
    cbind(temp_seq_var, 
          replicate(num_perm - 1, temp_seq_var$C171S_L172I_or_insertion))
  for (i in 2:ncol(temp_seq_var)) {
    colnames(imp_mat)[i] <- colnames(temp_seq_var)[i] <- paste0("imp", i - 1)
  }
  
  imp_and_seq_var <- rbind(imp_mat, temp_seq_var)
  all_data <- full_join(original_data, imp_and_seq_var, by = "ID")
  return(all_data)
} # end join_imputed_to_seq()

#' calculate_FE
#' 1. For the sequenced data calculate the FE for only the sequenced variants. 
#'   Save FE results. 2. For a tibble with the isolates in the rows and the 
#'   columns as the combo of real (observed) supplemented with imputed trehalose 
#'   utilization variants calculate the FE and return FE results. 
#' @param var_data Tibble. The isolates are in the rows and the columns are a 
#'   combination of real (observed) supplemented with imputed trehalose 
#'   utilization variants
#' @param n_perm Number. Number of imputations performed.
#'
#' @return FE_results. Matrix. Nrow = n_perm. Ncol = 4. Columns describe FE 
#'   results: "OR", "95% CI (lower)", "95% CI (upper)", and "P-value". 
#' @noRd
calculate_FE <- function(var_data, n_perm){
  # Sequenced, matched results
  seq_only <- var_data %>% filter(WGS_performed == 1, Stratum_complete == 1)
  seq_fe <- exact2x2(seq_only$Severe_Outcome, seq_only$C171S_L172I_or_insertion)
  seq_results <- matrix(NA, ncol = 5, nrow = 1)
  colnames(seq_results) <- c("OR", 
                             "95% CI (lower)", 
                             "95% CI (upper)",
                             "P-value", 
                             "Number Samples Included")
  seq_results[1, 1] <- format(round(seq_fe$estimate, 2), nsmall = 2)
  seq_results[1, 2] <- format(round(seq_fe$conf.int[1], 2), nsmall = 2)
  seq_results[1, 3] <- format(round(seq_fe$conf.int[2], 2), nsmall = 2)
  seq_results[1, 4] <- format(round(seq_fe$p.value, 2), nsmall = 2)
  seq_results[1, 5] <- nrow(seq_only)
  seq_results <- as_tibble(seq_results)
  write_tsv(seq_results, 
            path = paste0("../data/outputs/", 
                          Sys.Date(),
                          "_FE_only_sequenced_results.tsv"))
  
  # Imputed plus sequenced results
  var_data <- 
    var_data %>% 
    # select(-c(ID, 
    #           Ribotype, 
    #           WGS_performed, 
    #           Duplicated_Patient, 
    #           C171S_L172I_or_insertion))
    filter(!is.na(imp1))
  
  FE_results <- matrix(NA, ncol = 5, nrow = n_perm)
  colnames(FE_results)  <- c("OR", 
                             "95% CI (lower)", 
                             "95% CI (upper)",
                             "P-value", 
                             "Number Samples Included")
  for (i in 1:n_perm) {
    temp_col <- paste0("imp", i)
    temp_fe <- 
      exact2x2(var_data$Severe_Outcome, var_data[[as.character(temp_col)]])
    FE_results[i, 1] <- temp_fe$estimate
    FE_results[i, 2] <- temp_fe$conf.int[1]
    FE_results[i, 3] <- temp_fe$conf.int[2]
    FE_results[i, 4] <- temp_fe$p.value
    FE_results[i, 5] <- nrow(var_data)
    
  }
   
  write_tsv(as_tibble(FE_results), 
            path = paste0("../data/outputs/", 
                          Sys.Date(), 
                          "_FE_imputation_plus_sequenced_results.tsv"))
  
  return(FE_results)
} # end calculate_FE()

#' summarize_model_results
#' Generate summary statistics for FE test or logitistic regression from real + 
#'   imputed variants and save the results in a table. 
#'   
#' @param model_results model_results. Matrix. Nrow = n_perm. Ncol = 5. Columns
#'   describe model results: "OR", "95% CI (lower)", "95% CI (upper)", "P-value"
#'   and "Number Samples Included."
#' @param suffix Character. Description of test ("FE" or "logit"). Optional. 
#'
#' @return model_summary. Tibble. nrow = 1. Ncol = 6. Columns describe statistics 
#'   summarizing the larger matrix model_results. Columns described: "Median OR", 
#'   "Min. OR", "Max OR", "Median P-value", "Min. P-value", and "Max P-value". 
#' @noRd
summarize_model_results <- function(model_results, suffix = ""){
  model_results <- as_tibble(model_results)
  or_summary <- summary(model_results$OR)
  p_summary <- summary(model_results$`P-value`)
  model_summary <- matrix(NA, ncol = 7, nrow = 1)
  colnames(model_summary) <- c("Median OR",
                            "Min. OR",
                            "Max OR",
                            "Median P-value",
                            "Min. P-value",
                            "Max P-value", 
                            "Number Samples Included")
  model_summary <- as_tibble(model_summary)
  model_summary$`Median OR`[1] <- format(round(or_summary[3], 2), nsmall = 2)
  model_summary$`Min. OR`[1] <- format(round(or_summary[1], 2), nsmall = 2)
  model_summary$`Max OR`[1] <- format(round(or_summary[6], 2), nsmall = 2)
  model_summary$`Median P-value`[1] <- format(round(p_summary[3], 2), nsmall = 2)
  model_summary$`Min. P-value`[1] <- format(round(p_summary[1], 2), nsmall = 2)
  model_summary$`Max P-value`[1] <- format(round(p_summary[6], 2), nsmall = 2)
  model_summary$`Number Samples Included`[1] <- model_results[1, 5]
  
  write_tsv(model_summary, 
            path = paste0("../data/outputs/", 
                          Sys.Date(), 
                          "_imputation_plus_sequenced_results_summary", 
                          suffix, 
                          ".tsv"))
  return(model_summary)
} # end summarize_model_results()

#' describe_imputation_cohort
#' Save in a log file a basic summary of the size cohort included in the 
#'   imputation steps. 
#' @param final_data Tibble. Metadata of the isolates included in the 
#'   imputation steps. 
#'
#' @noRd
describe_imputation_cohort <- function(final_data){
  sequenced <- final_data %>% 
    filter(WGS_performed == 1, Stratum_complete == 1)
  imputed <- final_data %>% 
    filter(WGS_performed == 0 | Stratum_complete == 0) %>% 
    filter(!is.na(imp1))
  
  n_seq_sample <- nrow(sequenced)
  n_seq_ribo <- length(unique(sequenced$Ribotype))
  n_seq_case <- sequenced %>% filter(Severe_Outcome == 1) %>% nrow()
  n_seq_ctrl <- sequenced %>% filter(Severe_Outcome == 0) %>% nrow()
  
  n_imp_sample <- nrow(imputed)
  n_imp_ribo <- length(unique(imputed$Ribotype))
  n_imp_case <- imputed %>% filter(Severe_Outcome == 1) %>% nrow()
  n_imp_ctrl <- imputed %>% filter(Severe_Outcome == 0) %>% nrow()
  
  # n_unique_individ <- sum(!final_data$Duplicated_Patient)
  
  # Describe FE test -----------------------------------------------------------
  # Sequenced FE
  FE_seq_only <-
    final_data %>% 
    filter(WGS_performed == 1, Stratum_complete == 1)
  FE_seq_only_case <- 
    FE_seq_only %>% 
    filter(Severe_Outcome == 1) %>%
    nrow()
  FE_seq_only_ctrl <- 
    FE_seq_only %>% 
    filter(Severe_Outcome == 0) %>% 
    nrow()
  
  # Imputed + Sequenced FE
  # Imputed plus all sequenced results
  FE_all <- final_data %>% 
    filter(!is.na(imp1)) 
  FE_all_case <- 
    FE_all %>% 
    filter(Severe_Outcome == 1) %>% 
    nrow()
  FE_all_ctrl <-
    FE_all %>% 
    filter(Severe_Outcome == 0) %>% 
    nrow()
  
  # Describe logit -------------------------------------------------------------
  logit_input <- 
    final_data %>% 
    filter(!is.na(Risk_Score), !is.na(imp1))
  logit_case <- 
    logit_input %>% 
    filter(Severe_Outcome == 1) %>% 
    nrow()
  logit_ctrl <- 
    logit_input %>% 
    filter(Severe_Outcome == 0) %>% 
    nrow()

  sink(paste0("../data/outputs/", 
              Sys.Date(), 
              "_imputation_descriptive_stats", 
              suffix, 
              ".txt"))
  print("Imputation cohort")
  print(paste0("Sequenced & Match N = ", n_seq_sample))
  print(paste0("Sequenced & Match Ribotype N = ", n_seq_ribo))
  print(paste0("Sequenced & Match Case N = ", n_seq_case))
  print(paste0("Sequenced & Match Control N = ", n_seq_ctrl))
  
  print(paste0("Imputed N = ", n_imp_sample))
  print(paste0("Imputed Ribotype N = ", n_imp_ribo))
  print(paste0("Imputed Case N = ", n_imp_case))
  print(paste0("Imputed Control N = ", n_imp_ctrl))
  
  print("FE TEST ALL")
  print(paste0("FE Test N = ", nrow(FE_all)))
  print(paste0("FE Test Case N = ", FE_all_case))
  print(paste0("FE Test Control N = ", FE_all_ctrl))
  
  print("FE TEST Seq & Match only")
  print(paste0("FE Test N = ", nrow(FE_seq_only)))
  print(paste0("FE Test Case N = ", FE_seq_only_case))
  print(paste0("FE Test Control N = ", FE_seq_only_ctrl))
       
  print("Logit")
  print(paste0("Logit N = ", nrow(logit_input)))
  print(paste0("Logit Case N = ", logit_case))
  print(paste0("Logit Control N = ", logit_ctrl))
  
  
  # print(paste0("Unique inidividuals N = ", n_unique_individ))
  sink()
} # end describe_imputation_cohort()

#' calculate_logit
#' @description Perform a logistic regression. Response variable: severe 
#'   infection outcome. Predictors: Risk Score (as a continuous number) and the
#'   presence/absence of the three trehalose utilization variants as a single 
#'   variable (binary). 
#'
#' @param var_dat Tibble. All of the imputed + observed variant data. Rows are
#'   samples, columns are variants. Also contains other relevant model inputs. 
#'
#' @return model_results Output of the logistic regression.
#'   List of 6 data points: 
#'   Variant. Character. 
#'   OR. Numeric.
#'   95% CI (lower). Numeric. 
#'   95% CI (upper). Numeric. 
#'   P-value. Numeric. 
#'   Number Samples Included. Numeric. 
#' @export
#'
#' @examples
calculate_logit <- function(var_dat){
  # Rearrange code to have only the trehalose variants and the necessary columns
  # Remember that "imputed" columns have imputed data for non-sequenced results
  #  but have real data for sequenced results. 
  model_input <- var_dat %>% select(paste0("imp", 1:num_perm), 
                                    c("Risk_Score", 
                                      "Severe_Outcome", 
                                      "C171S_L172I_or_insertion"))
  # Drop isolates without a risk score or without variant data
  model_input <- model_input %>% filter(!is.na(Risk_Score), !is.na(imp1))
  print("dim of model input")
  print(dim(model_input))
  
  # Initialize data objects
  unadjusted_models <- rep(list(NULL), num_perm)
  imputation_index <- 1:num_perm
  
  # Compute logisitic regression based on variant and risk score 
  for (i in imputation_index) {
    unadjusted_models[[i]] <- 
      glm(formula = Severe_Outcome ~ model_input[ , i, drop = TRUE] + 
            Risk_Score, 
          data = model_input, 
          family = "binomial")
  }
  
  # Information stored in logit model: 
  # coef = logs odd ratio
  # exp(coef) = odds ratio
  # 95%CI = exp(confint(model))
  
  # Save logit results in a model results
  model_results <- matrix(NA, ncol = 6, nrow = num_perm)
  colnames(model_results) <- c("Variant",
                               "OR", 
                               "95% CI (lower)", 
                               "95% CI (upper)",
                               "P-value",
                               "Number Samples Included")
  for (m in imputation_index) {
    model_results[m, 1] <- colnames(model_input)[m]
    model_results[m, 2] <-
      format(round(exp(unadjusted_models[[m]][]$coefficients[2]), 2), 
             nsmall = 2)
    model_results[m, 3] <- 
      format(round(exp(confint(unadjusted_models[[m]])[2, 1]), 2), 
             nsmall = 2)
    model_results[m, 4] <- 
      format(round(exp(confint(unadjusted_models[[m]])[2, 2]), 2),
             nsmall = 2)
    model_results[m, 5] <- 
      format(round(coef(summary(unadjusted_models[[m]]))[2 ,'Pr(>|z|)'], 2), 
             nsmall = 2)
    model_results[m, 6] <- nrow(model_input)
  }
  
  # Order results by increasing OR
  model_results <- as_tibble(model_results)
  class(model_results$OR) <- "numeric"
  model_results <- model_results %>% arrange(OR)
  
  # Save table 
  write_tsv(x = model_results, 
            path = paste0("../data/outputs/", 
                          Sys.Date(), 
                          "_logit_with_score_imputation.tsv"))

  model_results <- as_tibble(model_results)
  model_results$OR <- as.numeric(model_results$OR)
  model_results$`95% CI (lower)` <- as.numeric(model_results$`95% CI (lower)`)
  model_results$`95% CI (upper)` <- as.numeric(model_results$`95% CI (upper)`)
  model_results$`P-value` <- as.numeric(model_results$`P-value`)

  # Get 10 representative model results
  # Evenly distributed results ranging from min to max point estimate
  # By 10% increments
  subset_row_index <- 
    round(nrow(model_results)/9) * seq(from = 1, to = 9, by  = 1)
  subset_row_index <- c(1, subset_row_index)
  if (max(subset_row_index) > nrow(model_results)) {
    subset_row_index[10] <- nrow(model_results)
  }
  
  subset_model_results <- subset_model_results[subset_row_index, , drop = FALSE]
   
  level_order <- 
    factor(subset_model_results$Variant, level = subset_model_results$Variant)
   
  CI_plot <- subset_model_results %>% 
    arrange(OR) %>% 
    ggplot(mapping = aes(x = level_order, y = OR)) + 
    geom_point() + 
    geom_errorbar(aes(ymax = `95% CI (upper)`, ymin = `95% CI (lower)`)) + 
    theme_bw() + 
    geom_abline(slope = 0, intercept = 1, color = "red") + 
    coord_flip() +
    xlab("Representative Imputations") + 
    ylab("Odds Ratio + 95% CI") + 
    ggtitle("OR, ordered from min(OR) to max(OR) by 10%") + 
    theme(axis.text.x = element_text(size = rel(3), colour = "black"), 
         axis.text.y = element_text(size = rel(0), color = "black"), 
         axis.title = element_text(size = rel(2)))
  
  ggsave(CI_plot,
         filename = paste0("../figures/", 
                           Sys.Date(),
                           "_logit_with_score_imputation_OR_and_CI.pdf"))
  
  pdf(paste0("../figures/", 
             Sys.Date(),
             "_logit_with_score_imputation_OR_hist.pdf"))
  hist(model_results$OR, 
       main = "Logit: Severe ~ Risk Score + Any Trehalose Variant",
       xlab = "OR", 
       breaks = 20)
  abline(v = 1.0, col = "red")
  dev.off()
  
  return(model_results)
} # end calculate_logit()

#' #' summarize_logit_results
#' #' @description Summarize the results of the logistic regression
#' #'
#' #' @param logit_results 
#' #' @param suffix 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' summarize_logit_results <- function(logit_results, suffix = ""){
#'   or_summary <- summary(logit_results$OR)
#'   class(logit_results$`P-value`) <- "numeric"
#'   p_summary <- summary(logit_results$`P-value`)
#'   logit_summary <- matrix(NA, ncol = 6, nrow = 1)
#'   colnames(logit_summary) <- c("Median OR",
#'                                "Min. OR",
#'                                "Max OR",
#'                                "Median P-value",
#'                                "Min. P-value",
#'                                "Max P-value")
#'   logit_summary <- as_tibble(logit_summary)
#'   logit_summary$`Median OR`[1] <- format(round(or_summary[3], 2), nsmall = 2)
#'   logit_summary$`Min. OR`[1] <- format(round(or_summary[1], 2), nsmall = 2)
#'   logit_summary$`Max OR`[1] <- format(round(or_summary[6], 2), nsmall = 2)
#'   logit_summary$`Median P-value`[1] <- 
#'     format(round(p_summary[3], 2), nsmall = 2)
#'   logit_summary$`Min. P-value`[1] <- format(round(p_summary[1], 2), nsmall = 2)
#'   logit_summary$`Max P-value`[1] <- format(round(p_summary[6], 2), nsmall = 2)
#'   
#'   write_tsv(logit_summary, 
#'             path = paste0("../data/outputs/", 
#'                           Sys.Date(), 
#'                           "_imputation_logit_results", 
#'                           suffix, 
#'                           ".tsv"))
#'   return(logit_summary)
#' } # end summarize_logit_results()


#' impute_variants
#' Wrapper function to perform trehalose utilization variant imputation on 
#'   isolates for which there is ribotype information but not sequence data. 
#' @param metadata_path Character. Path to metadata file. 
#' @param num_perm Number. Number of times to run imputation. 
#'
#' @noRd
impute_variants <- function(metadata_path, num_perm){
  variant_data <- read_tsv(metadata_path)
  # Select columns to work with
  variant_data <- 
    variant_data %>% 
    select(c(ID,
             Ribotype,
             Severe_Outcome, 
             WGS_performed,
             Duplicated_Patient, 
             C171S_L172I_or_insertion, 
             Stratum_complete, 
             Risk_Score))
  
  # Select isolates to work with
  # Remove isolates with Ribotypes that are "other", "Unique" or NA
  # Remove isolates that are from duplicated patients
  original_data <- variant_data
  matched_seq_data <- 
    variant_data %>% 
    filter(WGS_performed == 1, Stratum_complete == 1)
  variant_data <- 
    variant_data %>% 
    filter(Duplicated_Patient == 0,
           !(Ribotype %in% c("other", "Unique")), 
           !is.na(Ribotype))
  
  variant_data <- drop_ribo_absent_from_matched_WGS(variant_data)
  obs_tre_var_by_ribo_case <- 
    calculate_variant_ratio_by_ribo(variant_data, "case")
  obs_tre_var_by_ribo_ctrl <- 
    calculate_variant_ratio_by_ribo(variant_data, "control")
  num_predicted_per_ribo_case <-
    predict_num_with_var(variant_data, obs_tre_var_by_ribo_case, "case")
  num_predicted_per_ribo_ctrl <-
    predict_num_with_var(variant_data, obs_tre_var_by_ribo_ctrl, "control")
  imputed_var_case <- 
    generate_imputed_data(variant_data,
                          num_perm, 
                          num_predicted_per_ribo_case, 
                          "case")
  imputed_var_ctrl <- 
    generate_imputed_data(variant_data, 
                          num_perm, 
                          num_predicted_per_ribo_ctrl,
                          "ctrl")
  impute_var_all <- rbind(imputed_var_case, imputed_var_ctrl)
  
  imputed_and_matched <- 
    join_imputed_to_seq(original_data, matched_seq_data, num_perm, impute_var_all)
  fe_results <- calculate_FE(imputed_and_matched, num_perm)
  fe_summary <- summarize_model_results(fe_results, suffix = "_FE")
  logit_results <- calculate_logit(imputed_and_matched)
  logit_summary <- 
    summarize_model_results(logit_results, suffix = "_logit_risk_score")
  describe_imputation_cohort(imputed_and_matched)
} # end impute_variants()


