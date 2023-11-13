#' Source all R scripts in a folder 
#'
#' Source all R scripts in a folder 
#'
#'
#' @export
#'
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}