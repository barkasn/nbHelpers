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


#' Get the top level structure of an object
#' @param x the object to examine
#' @export str1
str1 <- function(x) {str(x,1)}

#' Get the two top levels of structure of an object
#' @param x the object to examine
#' @export str2
str2 <- function(x) {str(x,2)}

#' Reload a package from the library
#' @param pkg the package to reload
#' @export rldpkg
rldpkg <- function(pkg) {
    require('devtools')
    reload(inst(pkg))
}

#' Load a sparse matrix from an HDF5 file ensuring that the
#' second dimension is sorted
#' @param path location of the HDF files
#' @return a Matrix sparse matrix
#' @export load.hca.matrix
load.hca.matrix <- function(path) {
    require('rhdf5')
    require('Matrix')

    cat('Loading HDF5 file...');
    genes <- h5read(path, 'GRCh38/genes')
    gene.names <- h5read(path, 'GRCh38/gene_names')
    barcodes <- h5read(path, 'GRCh38/barcodes')
    indptr <- h5read(path, 'GRCh38/indptr')
    shape <- h5read(path, 'GRCh38/shape')
    indices <- h5read(path, 'GRCh38/indices')
    data <- h5read(path, 'GRCh38/data')
    cat('done\n');

    cat('Sorting elements...');
    ord1 <- unlist(lapply(1:(length(indptr)-1), function(i) {
        rowStart <- indptr[i] + 1;
        ## indptr[i+1] is the start of the next one in 0 indexed
        rowEnd <- indptr[i+1];
        rowIndices <- indices[c(rowStart:rowEnd)]
        rowElementOrder <- order(rowIndices) + rowStart - 1;
        rowElementOrder
    }))
    cat('done\n')

    cat('Generating matrix...')
    res <- new('dgCMatrix',
        i=as.integer(indices[ord1]),
        p=as.integer(indptr),
        x=as.double(data[ord1]),
        Dim=as.integer(shape),
        Dimnames=list(gene.names,barcodes)
        )
    cat('done\n')

    res
}

#' Get the symmetrical set difference between two sets (XOR)
#' @param a first set
#' @param b sencond set
#' @return elements in one of the two sets but not in the intersection
#' @export symmsetdiff
symmsetdiff <- function(a,b) {c(setdiff(a,b),setdiff(b,a))}

#' Convert DESeq2 results table into a tibble
#' @import tibble
#' @return results table in tibble format
#' @export deseqRes2Tibble
deseqRes2Tibble <- function(res) {
    if(class(res) != 'DESeqResults') stop('res is not DESeq2 results');
    res <- as.data.frame(res)
    res$gene <- rownames(res)
    as.tibble(res)
}

#' Get the specified substring of string delimited by a character
#' @param x the strings to breakdown now (can be vector)
#' @param split the character to breakdown by
#' @param integer specifying which substing to get
#' @return extracted substrings
#' @export strpart
strpart <- function(x, split, n, fixed=FALSE) {
    sapply(strsplit(as.character(x),split,fixed=fixed),'[',n)
}

#' Set the display parameter for X11 plotting
#' @param display.number the display number to set
#' @param host the host to set, defaults to localhost
#' return NUll
#' @export setDisplay
setDisplay <- function(display.number = NULL, host='localhost') {
    if(is.null(display.number)) stop('display number not specified');
    Sys.setenv(DISPLAY=paste0(host,':',display.number));
    invisible(NULL)
}

#' Get the display paramters for X11 plotting
#' @return DISPLAY variable value
#' @export getDisplay
getDisplay <- function() {
    Sys.getenv('DISPLAY');
}

#' Return TRUE if the parameter is an error object
#' @param x the object to test
#' @export is.error
is.error <- function(x) {
  inherits(x, c("try-error", "error"))
}

#' get number of rows and cols to use when printing a n.items number of items
#' @param n.items number of items
#' @param square force number of columns and rows to be equal
#' @export getParMfrow
getParMfrow <- function(n.items, square = FALSE) {
  n <- ceiling(sqrt(n.items))
  if (square)  {
    c(n,n);
  } else {
    m <- ceiling(n.items/n)
    c(n,m)
  }
}

#' Get the number of commong elements among the top n elements
#' of two character vectors where n varies between a specified range
#' @param itemsA character vector A
#' @param itemsB character vector B
#' @param start start of overlap sequence
#' @param seq_end end of overlap sequence,
#' if NULL up to the size of the smallest of the two lists (default: NULL)
#' @return dataframe with two columns
#' @export vectorSeqOverlap
vectorSeqOverlap <- function(itemsA, itemsB, start=1, seq_end=NULL) {
    ## Process params
    if( class(itemsA) != 'character') {
        stop('itemsA is not a character vector');
    }
    if ( class(itemsB) != 'character') {
        stop('itemsB is not a character vector');
    }
    if(!is.numeric(start)) {
        stop('start is not a number')
    }
    if(is.null(seq_end)) seq_end <- min(length(itemsA), length(itemsB))
    if(start > seq_end) {stop('start is larger than end')}
    ## Get sequence to calculater overlap
    s <- seq(start,seq_end)
    names(s) <-s
    ## Calculate the overlap at different levels down the list
    overlap <- unlist(lapply(s, function(n) {
        length(intersect(head(itemsA,n=n),head(itemsB,n=n)))
        }))
    ## return
    data.frame(
        row.names=NULL,
        n=as.numeric(names(overlap)),
        overlap=as.numeric(overlap)
    )
}

#' Convert NA values in x to the value specified by val
#' @export NA2VALUE
NA2VALUE <- function(x, val) {
  x[is.na(x)] <- c(val); x
}

#' Convert NA values in x to FALSE
#' @export NA2FALSE
NA2FALSE <- function(x) {
  NA2VALUE(x,FALSE)
}

#' Get the number of processors allocated to the current job by
#' the Slurm scheduler from the environment. Return 1 otherwise
#' @export getSlurmCpus
getSlurmCpus <- function() {
  NA2VALUE(as.numeric(Sys.getenv()['SLURM_CPUS_ON_NODE']),1)
}

#' Get a vector with the levels of a factor named with their
#' own name. Useful for lapply loops over factor levels
#' @param f a factor
#' @return named character vector of factor levels
#' @export namedLevels
namedLevels <- function(f) {
  if(!is.factor(f)) {
    stop('f is not a factor')
  }
  r <- levels(f)
  names(r) <- r
  r
}

#' Like list.files but return the full path and the vector entries
#' are named only after the basename
#' @param path location to list
#' @export getNamedDirectoryListing
getNamedDirectoryListing <- function(path) {
  x <- list.files(path)
  z <- paste0(path,'/',x)
  names(z) <- x
  z
}

#' Parallel apply using the packages parallel, BiocParallel or regular lapply
#' @description From the pagoda package with some minor improvements originally called papply
#' @param ... parameters to the apply functions, do not include BPPARAM or mc.cores
#' @param n.cores numbere of cores to use, will automatically detect the number of physical cores if possible
#' @param mc.preschedule parameter for mclapply
#' @export plapply
plapply <- function(...,n.cores=detectCores(logical=FALSE), mc.preschedule=FALSE) {
  if(n.cores>1) {
    if(requireNamespace("parallel", quietly = TRUE)) {
      return(parallel::mclapply(...,mc.cores=n.cores,mc.preschedule=mc.preschedule))
    }

    if(requireNamespace("BiocParallel", quietly = TRUE)) {
      # It should never happen because parallel is specified in Imports
      return(BiocParallel::bplapply(... , BPPARAM = BiocParallel::MulticoreParam(workers = n.cores)))
    }
  }

  # fall back on lapply
  lapply(...)
}

#' Get the percent of entries that are filled in a sparse matrix
#' @param mat matrix of class dgCMatrix or dgRMatrix
#' @return percent of entries that are non-zero
#' @export sparseMatrixDensity
sparseMatrixDensity <- function(mat) {
  if(!class(mat) %in% c("dgCMatrix","dgRMatrix")) stop('mat needs to be of class dgCMatrix or dgRMatrix')
  dim0 <- dim(mat)
  sum(mat@x > 0) / (dim0[1] * dim0[2])
}

#' Transpose a dgRMatrix and simultaneously convert it to dgCMatrix
#' @param inmat input matrix in dgRMatrix format
#' @return A dgCMatrix that is the transposed dgRMatrix
#' @export transpose_dgRMatrix
transpose_dgRMatrix <- function(inmat) {
    if(class(inmat) != 'dgRMatrix')
        stop('inmat is not of class dgRMatrix')
    out <- new('dgCMatrix',
               i=inmat@j,
               p=inmat@p,
               x=inmat@x,
               Dim=rev(inmat@Dim),
               Dimnames=rev(inmat@Dimnames)
               )
    out
}                                 

#' Given two vectors return the number of elements
#' only in a, in their intersection and only in b
#' @param a first vector
#' @param b second vector
#' @return a vector with (element count only in a, element count in intersection, element count only in b)
#' @export vennCounts
vennCounts <- function(a,b) {
    i <- length(intersect(a,b))
    c(length(a) - i, i, length(b))
}
