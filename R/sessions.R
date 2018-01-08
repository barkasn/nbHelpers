#' @importFrom utils read.table
#' @importFrom Matrix readMM
NULL

#' Save the session in the current working directory with a file name that
#' includes the time stamp and process id
#' @param prefix the prefix to use,for the filename
#' @return the name of the file that was saved in
#' @export preserve.state
preserve.state <- function(prefix='savepoint_') {
  file <- paste0(prefix,gsub(' ','_',Sys.time()),'_',Sys.getpid(),'.RDataF')
  save.image.fast(file)
  file
}

#' @title save a R session image fast
#' @description saves an R image of the current session much faster and with
#' better compression that the built-in save.image function. Requires system lbzip2
#' utility
#' @param filename the file to save the image to
#' @param tmpfile optional temporary file name, using a ramdisk location will further accelerate saving
#' @param verbose logical verbosity level
#' @export save.image.fast
save.image.fast <- function(filename,tmpfile = NULL, verbose = FALSE) {
  if (is.null(filename)) {
    stop("filename argument is required");
  }
  if (is.null(tmpfile)) {
    tmpfile <- tempfile();
  }

  if (verbose) { cat('Saving to temp file ', tmpfile,'\n') }
  save.image(tmpfile, compress = FALSE);

  cmd <- paste0('lbzip2 -c ', tmpfile, ' > ', filename)

  if (verbose) { cat('Compressing...\n') };
  system(cmd);

  if (verbose) { cat('Deleting temporary file...') }
  unlink(tmpfile)
}

#' Loads an image generated with save.image.fast()
#' @description loads an image generated with the save.image.fast function. This function
#' requires that the system has the lbunzip2 command installed
#' @param filename filename to load
#' @param tmpfile temporary file to use, use of ramdisk file will accelerate loading
#' @param verbose logical verbosity level
#' @param envir enviroment in which to load the data, by default the calling environment
#' @export load.image.fast
load.image.fast <- function(filename, tmpfile = NULL, verbose=F, envir = parent.frame()) {
  if (is.null(filename)) {
    stop('filename argument is required');
  }
  if (is.null(tmpfile)) {
    tmpfile <- tempfile();
  }

  if (verbose) { cat('Decompressing into temporary file ', tmpfile, '...\n') }
  cmd <- paste0('lbunzip2 -c ', filename, ' > ', tmpfile);
  system(cmd);

  if  (verbose) { cat('Loading...\n') }
  load(tmpfile, envir = envir);

  if (verbose) { cat('Deleting temporary file...\n') }
  unlink(tmpfile);
}
