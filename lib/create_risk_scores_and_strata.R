# Adapted Krishna Rao's script into a function and updated outputs. 

#' Create the risk scores and group isolates into strata 
#' @description Generate the severe outcome risk scores and matched strata from 
#'   model in Rao 2015 CID paper for trehalose analysis.  
#' @param model_data Matrix or dataframe with model data.  
#'
#' @return Dataframe with dataset's propensity scores, stratum (when 
#'   applicable), and original data that went into creating model. 
#' @noRd
create_risk_scores_and_strata <- function(model_data){
  model_param <- c("CDC_SEV_ATT",
                   "age", 
                   "gender..M.0.F.1.", 
                   "METS", 
                   "concurrentabx", 
                   "Lowest_SBP", 
                   "highcreat", 
                   "highbili", 
                   "WBC")
  model_data <- 
    model_data %>% 
    mutate("Duplicated_Patient" = duplicated(model_data$PATIENT_NUM))
  
  # Remove isolates with missing clinical data
  incomplete_ids <- which(colnames(model_data) %in% model_param) 
  missing_data <-
    model_data[!complete.cases(model_data[ , incomplete_ids]), , drop = FALSE]
  missing_data <- missing_data %>% mutate("Missing_Model_Data" = TRUE)
  complete_data <- 
    model_data[complete.cases(model_data[ , incomplete_ids]), , drop = FALSE]
  complete_data <- complete_data %>% mutate("Missing_Model_Data" = FALSE)
  
  # Generate logistic regression
  model <- glm(CDC_SEV_ATT ~ 
                 age + 
                 gender..M.0.F.1. + 
                 METS + 
                 concurrentabx + 
                 Lowest_SBP + 
                 highcreat + 
                 highbili + 
                 WBC, 
               family = binomial(logit), 
               data = complete_data)

  # Generate propensity scores
  complete_data$scores <- predict(model, type = c("response"))

  duplicated_ptx <- 
    complete_data[duplicated(complete_data$PATIENT_NUM), , drop = FALSE]
  duplicated_ptx <- 
    duplicated_ptx %>% mutate("Duplicated_Patient_No_Missing_Info" = TRUE)
  complete_data <- 
    complete_data[!duplicated(complete_data$PATIENT_NUM), , drop = FALSE]
  complete_data <- 
    complete_data %>% mutate("Duplicated_Patient_No_Missing_Info" = FALSE)
  
  # Matching  
  cases <- complete_data[complete_data$CDC_SEV_ATT == 1, , drop = FALSE]
  cases <- arrange(cases, desc(scores))
  
  controls <- complete_data[complete_data$CDC_SEV_ATT == 0, , drop = FALSE]
  controls <- arrange(controls, desc(scores))
  controls$diff <- as.numeric(NA)
  
  cases <- arrange(cases, desc(CDC_SEV_ATT), desc(scores))
  cases$stratum <- as.integer(NA)
  cases$diff <- as.numeric(NA)
  
  cases <- as.data.frame(cases)
  controls <- as.data.frame(controls)
  
  for (i in 1:length(cases$SAMPLE_ID)) {
    cases[i, "stratum"] <- i
    scores <- cases[i, "scores"]
    tempframe <- controls[abs(scores - controls$scores) < .1, , drop = FALSE]
    tempframe$diff <- abs(scores - tempframe$scores)
    tempframe <- arrange(tempframe, diff)
    if (length(tempframe$SAMPLE_ID) > 3) { 
      j <- 4 
    } else { 
      j <- length(tempframe$SAMPLE_ID) 
    }
    
    tempframe <- tempframe[1:j, ]
    tempframe$stratum <- i
    cases <- rbind(cases, tempframe)
    controls <- controls[!(controls$SAMPLE_ID %in% tempframe$SAMPLE_ID), ]
  }
  
  # Some cases will not have a matching control; let's eliminate those. 
  nomatch <- cases[is.na(cases$PATIENT_NUM), "stratum"]
  unmatched_cases <- cases %>% filter(stratum %in% nomatch, !is.na(SAMPLE_ID))
  
  temp <- cases[!(cases$stratum %in% nomatch), , drop = FALSE]
  temp <- temp %>% mutate("Unmatched" = FALSE)
  
  controls <- controls %>% mutate("stratum" = NA)
  unmatched <- rbind(unmatched_cases, controls)
  unmatched <- unmatched %>% mutate("Unmatched" = TRUE)
  
  # Subset data to add only important columns
  unmatched <- unmatched %>% select(-PATIENT_NUM, -stratum, -diff)
  missing_data <- missing_data %>% select(-PATIENT_NUM)
  duplicated_ptx <- duplicated_ptx %>% select(-PATIENT_NUM)
  
  # Include all important information in output table
  temp <- select(temp, c(-diff, -PATIENT_NUM))
  foo <- full_join(temp, unmatched)
  foo <- full_join(foo, duplicated_ptx)
  final <- full_join(foo, missing_data)

  colnames(final) <- c("ERIN_ID", 
                       "Severe_Outcome",
                       "age", 
                       "gender..M.0.F.1.", 
                       "METS", 
                       "concurrentabx", 
                       "Lowest_SBP", 
                       "highcreat", 
                       "highbili", 
                       "WBC", 
                       "Ribotype", 
                       "Duplicated_Patient",
                       "Missing_Model_Data", 
                       "Risk_Score", 
                       "Duplicated_Patient_No_Missing_Info", 
                       "Stratum", 
                       "Unmatched")
  
  write_tsv(final, path = "../data/outputs/severity_model.tsv")
} # end create_risk_scores_and_strata()