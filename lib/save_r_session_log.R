#' save_session_info
#' Save the R session info, loaded R libaries to text file. 
#'
#' @noRd
save_session_info <- function(){
  log <- utils::capture.output(utils::sessionInfo())
  write(log, 
        file = paste0("../data/outputs/", 
                      Sys.Date(), 
                      "_R_session_log.txt"))
} # end save_session_info()
