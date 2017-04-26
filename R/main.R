#' @author Nikolas Barkas

#'
#'
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'




#' @title  Read 10x matrix
#' @description This function reads a 10x matrix from the specified directory.
#' @param path location of 10x output
#' @return None
#' @importFrom utils read.table
#' @importFrom Matrix readMM
#' @export
read10xMatrix <- function(path) {
  matrixFile <- paste(path, 'matrix.mtx', sep='/');
  genesFile <- paste(path, 'genes.tsv', sep='/');
  barcodesFile <- paste(path, 'barcodes.tsv', sep='/');

  if (!file.exists(matrixFile)) { stop('Matrix file does not exist');  }
  if (!file.exists(genesFile)) { stop('Genes file does not exist'); }
  if (!file.exists(barcodesFile)) { stop('Barcodes file does not exist'); }

  x <- as.matrix(readMM(matrixFile));
  genes <- read.table(genesFile)
  rownames(x) <- genes[,2];
  barcodes <- read.table(barcodesFile);
  colnames(x) <- barcodes[,1]
  storage.mode(x) <- 'integer'
  invisible(x);
}
