#' Save R session information
#' @description Save the R session info & libraries to text file. 
#'
#' @noRd
save_session_info <- function(){
  log <- utils::capture.output(utils::sessionInfo())
  write(log, 
        file = "../data/outputs/R_session_log.txt")
} # end save_session_info()
