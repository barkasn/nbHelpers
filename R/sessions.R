#' @importFrom utils read.table
#' @importFrom Matrix readMM
NULL

#' Save the session in the current working directory with a file name that
#' includes the time stamp and process id
#' @param prefix the prefix to use,for the filename
#' @return the name of the file that was saved in
#' @importFrom fastSave save.image.fast
#' @export preserve.state
preserve.state <- function(prefix='savepoint_') {
  file <- paste0(prefix,gsub(' ','_',Sys.time()),'_',Sys.getpid(),'.RDataF')
  save.image.fast(file)
  file
}
