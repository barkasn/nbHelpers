#' @importFrom utils read.table
#' @importFrom Matrix readMM
NULL

#' Get all unique combinations of elements in a list
#' @param nms elements to get combinations
#' @return list of pairs of elements
#' @importFrom gtools combinations
#' @export uniqCombs
uniqCombs <- function(nms) {
  nms <- unique(nms)

  combs <- combinations(n = length(nms), r = 2, v = nms, repeats.allowed =F)
  ret <- split(t(combs), rep(1:nrow(combs), each=ncol(combs)))

  ret
}

#' Head and tail of a table combined
#' @param x the table object
#' @param n a single integer
#' @return a table made of the n first elemets of x and n last elements of x
headtail <- function(x, n=6) {
  rbind(head(x,n), tail(x,n))
}

#' @title Find minimum string distance strings in a vector
#' @description Find the minimum string distance between strings in a vector,
#' this is useful for things like finding the number of allowable barcode
#' mismatches
#' @importFrom stringdist stringdistmatrix
#' @export minStringDistance
minStringDistance <- function(strings) {
  ds <- stringdist::stringdistmatrix(strings,strings);
  min(ds[upper.tri(ds)])
}

#' @title get NW, NE, SW and SE corners of a dataframe or matrix
#' @description return a new dataframe or matrix with the elements at the
#' four corners of the one provided
#'
#' @param vals the dataframe or matrix of interest
#' @param n how many rows/cols to return from each end
#'
#' @return the new matrix
#'
#' @export tableCorners
tableCorners <- function(vals, n = 4) {
  if (is.null(vals)) error("vals is NULL");
  if (n < 1) error("Invalid value for n");

  d <- dim(vals)
  rows <- c(1:n,(d[1]-n):(d[1]))
  cols <- c(1:n,(d[2]-n):(d[2]))
  vals[rows, cols]
}

#' Generates a new vector from two vectors by combining them according to a third
#' @description generates a vector from two vectors depending of the logical value of
#' a third vector
#' @param logical vector to use to combine the vTrue and vFalse vectors
#' @param vTrue values to return if the logical is true
#' @param vFalse values to return if the logical is false
#' @return a new vector with some values from vTrue and some for vFalse
#' @export mixVectorsOnLogical
mixVectorsOnLogical <- function(logical,vTrue,vFalse) {
  mapply(function(v, vt, vf) {ifelse(v, vt, vf)}, logical, vTrue, vFalse);
}

#' Convert list of lists of numbers to array
#' @description Converts lists of lists of numbers to a numeric array fast
#' this has been optimised against other possible implementations and
#' been found to the the fastest
#' @param inputList a list of lists of numbers
#' @return an array where the top-level lists are rows and low-level lists are columns
#' @examples
#' mylist <- list()
#' mylist_ <- list()
#' for(i in 1:10) {
#'   for(j in 1:100) {
#'     mylist[[j]] <- i*j
#'   }
#'   mylist_[[i]] <- mylist
#' }
#' str(mylist_)
#' r3 <- ll2a.3(mylist_)
#' @export ll2a
ll2a <- function(inputList){
  f <- as.matrix(as.data.frame(sapply(inputList, function(x) (unlist(x)),simplify = F),stringsAsFactors = F,as.is=T))
  rownames(f) <- NULL
  colnames(f) <- NULL

  t(f)
}

#' Write a list of data.frames as individual csv files
#' @description Write a list of data.frames as individual csv files
#' @param x named list of data.frames
#' @param output.dir the directory to save the files in
#' @export writeListOfDFs
writeListOfDFs <- function(x, output.dir) {
  # TODO: Add checks: list is named; elements are data.frames; output.dir exists

  mapply(function(x,name) {
    write.table(x, paste0(output.dir,'/',name,'.csv'))
  }, deres, names(deres));

  invisible(NULL)
}

#' Get a vector of the names of an object named by the names themselves
#' @description Get a named vector of the names of an object. This is useful
#' with lapply when passing names of objects as it ensures that the output list
#' is also named
#' @param g an objects on which we can call names()
#' @export namedNames
namedNames <- function(g) {
  n <- names(g)
  names(n) <- n;
  n
}

