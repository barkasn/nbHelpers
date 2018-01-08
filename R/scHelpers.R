#' @importFrom utils read.table
#' @importFrom Matrix readMM
NULL

#' Merge expression matrices in a list, keeping only common genes
#' @param data.all a list of sparse matrices
#' @param prefix.names logical, denoting whether to prefix the cell names with the table name
#' @param prefix.sep character, separator to use between sample and cell name, defaults to _
#' @return a single merged matrix
#' @export mergeMatrices
mergeMatrices <- function(data.all, prefix.names=T, prefix.sep = '_') {
  common.genes <- Reduce(intersect, lapply(data.all, function(x) {rownames(x)}))

  data.all.common <- mapply(function(x, name) {
    y <- x[common.genes,]
    if(prefix.names) {
      colnames(y) <- paste0(name, prefix.sep,  colnames(y))
    }
    y
  }, data.all, names(data.all))
  x <- do.call(cbind, data.all.common)

  x
}

#' Get the colors of a pagoda1 serialised app object as a factor
#' @description get the colors of a pagoda1 serialised app object as a factor from the app$results$colcol list
#' @param filename the rds file to load
#' @param colcol.name the name of the colcol item to retrive
#' @return a named factor
getFactorFromP1rds <- function(filename, colcol.name) {
  app <- readRDS(filename)
  env <- app@.xData
  app2 <- as.list(env)
  colcol <- app2$results$colcol
  fac <- colcol[[colcol.name]]$data
  fac
}
