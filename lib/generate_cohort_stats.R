#' print_ribo_stats
#' This function prints to a log information about the number of isolates from a 
#'   specified ribotype in the sequenced cohort, sequenced cases, and sequenced
#'   controls. 
#' @param metadata. Tibble with ribotype, sequencing status, and severity info. 
#' @param rbtyp Character. One of the ribotypes in the cohort. 
#'
#' @noRd
print_ribo_stats <- function(metadata, rbtyp){
  print(paste0("Num ",rbtyp, " in sequenced cohort"))
  print(nrow(metadata %>% 
               filter(WGS_performed == 1, Ribotype == rbtyp)))
  print(paste0("Num ", rbtyp, " in sequenced cases"))
  print(nrow(metadata %>% 
               filter(WGS_performed == 1, 
                      metadata$Severe_Outcome == 1, 
                      Ribotype == rbtyp)))
  print(paste0("Num ", rbtyp, " in sequenced controls"))
  print(nrow(metadata %>% 
               filter(WGS_performed == 1, 
                      metadata$Severe_Outcome == 0, 
                      Ribotype == rbtyp)))
  
} # end print_ribo_stats()

#' print_phylogenetic_distribution
#' This function prints to a log information about the number of isolates with a
#'   specific variant and their ribotypes to gain insight into the phylogenetic
#'   distribution of the variant.  
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
  print(paste0(variant_name, " distribution"))
  print(paste0(ribotype_num, 
               " with ", 
               variant_name, 
               " / total ", 
               ribotype_num))
  variant_in_ribo  <- 
    metadata %>% 
    filter(Ribotype == ribotype_num, 
           WGS_performed == TRUE) %>% 
    select(variant_name)
  print(paste0(sum(variant_in_ribo[ , 1]), "/", nrow(variant_in_ribo)))
  
  print(paste0("All isolates with ", variant_name, " / total isolates"))
  all_with_variant <- 
    metadata %>% 
    filter(WGS_performed == TRUE) %>% 
    select(variant_name)
  print(paste0(sum(all_with_variant[ , 1]), "/", nrow(all_with_variant)))
  
  ribotypes_with_variant <- 
    metadata %>% 
    filter(WGS_performed == TRUE) %>% 
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

#' generate_stats
#' This function generates the basic stats required to describe the propsenity
#'   score matched data set. 
#'
#' @param metadata_path Character. Path to location of metadata file. 
#'
#' @noRd
generate_stats <- function(metadata_path){
  metadata <- read_tsv(metadata_path, 
                       col_names = TRUE)
  sink(paste0("../data/outputs/", Sys.Date(), "_descriptive_stats.txt"))
  # Study population: 
  # Number of CDI cases
  print("Number of CDI cases")
  print(nrow(metadata))
  # Number of individual patients
  print("Number of individual patients")
  print(nrow(metadata %>% filter(Duplicated_Patient == FALSE)))
  # Number of cases with complete data and from individual patients
  print("Number of cases with complete model data and from individual patient")
  print(nrow(metadata %>% filter(Duplicated_Patient_No_Missing_Info == FALSE)))
  # Number of total ribotypes (all 1144 cases)
  print("Number of ribotypes (N=1144)")
  print(length(unique(metadata$Ribotype)) + 
                 sum(metadata$Ribotype == "Unique", na.rm = TRUE) + 
                 sum(is.na(metadata$Ribotype), na.rm = FALSE) -
          2) # Don't count NA and unique twice: also in length(unique(Ribotype))
  # Number of isolates given a propensity score
  print("Number of isolates with a propensity score")
  print(sum(!is.na(metadata$Propensity_Score)))
  
  # Number of isolates with WGS sequencing
  print("Size of final matched set")
  print(sum(metadata$WGS_performed))
  
  # Number of cases
  print("Num cases")
  print(nrow(metadata %>% 
               filter(WGS_performed == 1, metadata$Severe_Outcome == 1)))
  
  # Number of controls
  print("Num controls")
  print(nrow(metadata %>% 
               filter(WGS_performed == 1, metadata$Severe_Outcome == 0)))
  
  # Controls / case stats
  case_ctrl_stat <-
    metadata %>% 
    filter(WGS_performed == 1) %>% 
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

  # number of ribotypes in the sequenced, matched cohort
  matched_cohort <- 
    metadata %>% 
    filter(WGS_performed == 1)
  
  num_ribo_in_matched <- 
    length(unique(matched_cohort$Ribotype)) + 
    sum(matched_cohort$Ribotype == "Unique", na.rm = TRUE) + 
    sum(is.na(matched_cohort$Ribotype), na.rm = FALSE) -
    2 # Don't count NA and unique twice: also in length(unique(Ribotype))
  
  # Phylogenetic distribution
  print_phylogenetic_distribution(metadata, 
                                  "Cys171Ser", 
                                  "017", 
                                  num_ribo_in_matched)
  print_phylogenetic_distribution(metadata, 
                                  "Leu172Ile", 
                                  "027", 
                                  num_ribo_in_matched)
  print_phylogenetic_distribution(metadata, 
                                  "four_gene_insertion",
                                  "078-126", 
                                  num_ribo_in_matched)
  sink()
} # end generate_stats()
