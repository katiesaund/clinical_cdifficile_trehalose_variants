#' Print statistics related to ribotype
#' @description This function prints to a log information about the number of 
#'   isolates from a specified ribotype in the sequenced cohort, sequenced 
#'   cases, and sequenced controls. 
#' @param metadata. Tibble with ribotype, sequencing status, and severity info. 
#' @param rbtyp Character. One of the ribotypes in the cohort. 
#'
#' @noRd
print_ribo_stats <- function(metadata, rbtyp){
  if (sum(metadata$WGS_performed & metadata$Stratum_complete) != nrow(metadata)) {
    stop("Assumes all included isolates were sequenced and from complete strata")
  }
  print(paste0("Num ", rbtyp, " in matched, case-control"))
  print(nrow(metadata %>% 
               filter(Ribotype == rbtyp)))
  print(paste0("Num ", rbtyp, " in matched cases"))
  print(nrow(metadata %>% 
               filter(Severe_Outcome == 1, 
                      Ribotype == rbtyp)))
  print(paste0("Num ", rbtyp, " in matched controls"))
  print(nrow(metadata %>% 
               filter(Severe_Outcome == 0, 
                      Ribotype == rbtyp)))
  
} # end print_ribo_stats()

#' Print phylogenetic distribution of this variant
#' @description This function prints to a log information about the number of 
#'   isolates with a specific variant and their ribotypes to gain insight into 
#'   the phylogenetic distribution of the variant.  
#' @param metadata Tibble with ribotype, sequencing status, and severity info. 
#' @param variant_name Character. Column name of the variant of interest found
#'   in metadata. 
#' @param ribotype_num Character. Ribotype name. 
#' @param num_ribo_in_matched Number. Number of ribotypes in the matched cohort
#'   that have the ribotype specified in ribotype_num. 
#'
#' @noRd
print_phylogenetic_distribution <- 
  function(metadata, variant_name, ribotype_num, num_ribo_in_matched){
  if (sum(metadata$WGS_performed & metadata$Stratum_complete) != nrow(metadata)) {
    stop("Assumes all included isolates were sequenced and from complete strata")
  }

  print(paste0(variant_name, " distribution"))
  print(paste0(ribotype_num, 
               " with ", 
               variant_name, 
               " / total ", 
               ribotype_num))
  variant_in_ribo  <- 
    metadata %>% 
    filter(Ribotype == ribotype_num) %>% 
    select(variant_name)
  print(paste0(sum(variant_in_ribo[ , 1]), "/", nrow(variant_in_ribo)))
  
  print(paste0("All isolates with ", variant_name, " / total isolates"))
  all_with_variant <- 
    metadata %>% 
    select(variant_name)
  print(paste0(sum(all_with_variant[ , 1]), "/", nrow(all_with_variant)))
  
  ribotypes_with_variant <- 
    metadata %>% 
    select(variant_name, Ribotype)
  ribotypes_with_variant <- 
    ribotypes_with_variant[ribotypes_with_variant[ ,  1] == 1, ]
  print(paste0("number of ribotypes with ",
               variant_name,
               " mutation / total"))
  print(table(ribotypes_with_variant$Ribotype))
  print(paste0(length(unique(ribotypes_with_variant$Ribotype)), 
               "/", 
               num_ribo_in_matched))
} # end print_phylogenetic_distribution()

#' Generate descriptve statistics of samples
#' @description This function generates the basic stats required to describe the 
#'   risk score matched data set. 
#'
#' @param metadata_path Character. Path to location of metadata file. 
#'
#' @noRd
generate_stats <- function(metadata_path){
  metadata <- read_tsv(metadata_path, 
                       col_names = TRUE)
  sink("../data/outputs/descriptive_stats.txt")
  # Study population: 
  # Number of CDI cases
  print("Number of C. difficile infections (original cohort)")
  print(nrow(metadata))
  
  # Number of individual patients
  print("Number of individual patients")
  print(nrow(metadata %>% filter(Duplicated_Patient == FALSE)))
  
  # Number of cases with complete data and from individual patients
  print("Number of cases with complete model data and from individual patient")
  print(nrow(metadata %>% filter(Duplicated_Patient_No_Missing_Info == FALSE)))
  
  # Number of total ribotypes (all 1144 cases)
  print("Number of ribotypes (N=1144)")
  print(length(unique(metadata$Ribotype))) 
  
  # Number of isolates given a propensity score
  print("Number of isolates with a propensity score")
  print(sum(!is.na(metadata$Risk_Score)))
  
  print("Number of isolates with placed in a stratum")
  print(sum(!is.na(metadata$Stratum)))
  
  # Number of isolates with WGS sequencing
  print("Number of isolates in a stratum that were sequenced")
  print(sum(metadata$WGS_performed))
  
  # Number of cases
  print("Total number of cases")
  print(nrow(metadata %>% 
               filter(WGS_performed == 1, metadata$Severe_Outcome == 1)))
  
  # Number of controls
  print("Total number of controls")
  print(nrow(metadata %>% 
               filter(WGS_performed == 1, metadata$Severe_Outcome == 0)))
  
  # Now only work with "good strata"
  print("Strata excluded because strata didn't have cases:")
  print(metadata %>% 
          filter(metadata$Stratum_complete == 0 & metadata$WGS_performed == 1) %>% 
          select(Stratum) %>% 
          unique() %>%
          unname() %>% 
          unlist())
  print("Isolates excluded because strata didn't have cases:")
  print(metadata %>% 
          filter(metadata$Stratum_complete == 0 & metadata$WGS_performed == 1) %>% 
          select(ID) %>% 
          unique() %>%
          unname() %>% 
          unlist())
  
  print("Below only with strata with at least 1 control and 1 case")
  metadata <- metadata %>% 
    filter(metadata$Stratum_complete == 1 & metadata$WGS_performed == 1)
  
  # Number of samples 
  print("Number of isolates sequenced, in complete strata, and therefore included in the matched analyses")
  print(nrow(metadata))
  
  print("Number of ribotypes for which samples were sequenced and in complete strata, and therefore included in the matched analyses")
  ribos <- unique(metadata$Ribotype)
  num_ribos <- length(ribos)
  print(num_ribos)
  
  print("All ribotypes in matched & sequenced")
  print(ribos)
  
  print("Number of ribotypes minus NA, unique, and other")
  ribos <- unique(metadata$Ribotype)
  ribos <- ribos[!ribos %in% c(NA, "Unique", "other")]
  print(length(unique(ribos)))
  
  # Number of cases
  print("Number of cases")
  print(nrow(metadata %>% 
               filter(metadata$Severe_Outcome == 1)))
  
  # Number of controls
  print("Number of controls")
  print(nrow(metadata %>% 
               filter(metadata$Severe_Outcome == 0)))
  
  # Controls / case stats
  case_ctrl_stat <-
    metadata %>% 
    select(Severe_Outcome, Stratum)
  
  ctrl_stat <- 
    case_ctrl_stat %>% 
    group_by(Stratum) %>% 
    mutate("case_in_stratum" = sum(Severe_Outcome == 1), 
           "ctrl_in_stratum" = sum(Severe_Outcome == 0)) %>% 
    filter(case_in_stratum == 1) %>% 
    select(-case_in_stratum)
  
  ctrls_per_stratum <- ctrl_stat %>% select(ctrl_in_stratum) 
  summary(ctrls_per_stratum$ctrl_in_stratum)
  
  print("Controls per stratum")
  print(summary(ctrls_per_stratum$ctrl_in_stratum))
  
  # Ribotype distribution
  print_ribo_stats(metadata, "014-020")
  print_ribo_stats(metadata, "015")
  print_ribo_stats(metadata, "017")
  print_ribo_stats(metadata, "027")
  print_ribo_stats(metadata, "053-160")
  print_ribo_stats(metadata, "078-126")
  
  # Phylogenetic distribution
  print_phylogenetic_distribution(metadata, 
                                  "Cys171Ser", 
                                  "017", 
                                  num_ribos)
  print_phylogenetic_distribution(metadata, 
                                  "Leu172Ile", 
                                  "027", 
                                  num_ribos)
  print_phylogenetic_distribution(metadata, 
                                  "four_gene_insertion",
                                  "078-126", 
                                  num_ribos)
  sink()
} # end generate_stats()
