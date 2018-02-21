#' @importFrom utils read.table
#' @importFrom Matrix readMM
NULL

#' Break down a factor returning the names of the elements
#' in each level as character vectors in a list
#' @param f a factor to breakdown
#' @return a list of factor levels with the names of elemements in them
#' @export factorBreakdown
factorBreakdown <- function(f) {
  if(!is.factor(f)) stop('not a factor!')
  lvls <- levels(f);
  names(lvls) <- lvls;
  lapply(lvls, function(l) {
    r <- names(f)[f == l]
    names(r) <- r
    r
  })
}

#' Break down a factor returning the names of the elements
#' in each level as character vectors in a list
#' @param f a factor to breakdown
#' @return a list of factor levels with the names of elemements in them
#' @export factorBreakdown.2
factorBreakdown.2 <- function(f) {tapply(names(f),f, identity) }

#' Convert a factor to a character vector, while preserving the names
#' @param f the factor to convert
#' @return a named character vector
#' @export factor2Char
factor2Char <- function(f) {
  r <- as.character(f);
  names(r) <- names(f);
  r
}

#' Combine two named factors when the new factor contains a breaking up of the
#' a subset of the original clusters
#' @param originalFactor the factor in which to breakup clusters, this is not modified
#' @param newFactor the new factor specifying new cluster identities for some of the elements
#' @return a new factor with the above combination of the factos
replaceClusterFactors <- function(originalFactor, newFactor, newPrefix=NULL) {
  if(is.null(newPrefix)) stop('newPrefix is null');
  wf <- factor2Char(originalFactor);
  nwf <- factor2Char(newFactor);
  if(any(!names(nwf) %in% names(wf))) stop('newFactor is not a subset of originalFactor');
  wf[names(nwf)] <- paste0(c(newPrefix), nwf);
  wf <- as.factor(wf);
  wf
}

#' Generate a ggplot2 plot comparing two factors for a common set of datapoints
#' @description Generates a ggplot barplot that compares two names two nameed factors
#' for the specified points. The x axis corresponds to factor B and the y axis
#' factor A.
#' @param factorA the factor the distibution of which with respect to factorB will be assessed (y axis)
#' @param factorB the factor to breakdown factorA by
#' @param points the names of the elemets of the factorsA and factorB to include for the
#' purposed of the plot
#' @examples
#' fA <- rep(c('A','B','C','D'), each = 4);
#' fB <- rep(c('W','Y'), each = 8)
#' n <- paste0('n',1:16)
#' names(fA) <- n
#' names(fB) <- n
#' plotFactorsPercent(fA, fB, points = names(fA)[1:9])
plotFactorsPercent <- function(factorA, factorB, points = NULL, factorA.name = 'factorA', factorB.name = 'factorB' ) {
  require(ggplot2)
  require(reshape2)

  if (is.null(points)) {
    if (!(all(names(factorA) %in% names(factorB)) & all(names(factorB) %in% names(factorA)))) {
      points <- intersect(names(factorA), names(factorB))
      warning("points is null and the two factors contain different elements: using only common elements")
    } else {
      points <- names(factorA)
    }
  }

  tmpdf <- data.frame(factorA=factorA[points], factorB = factorB[points]);
  tmpa <- acast(tmpdf, factorA ~ factorB, fun.aggregate=length, value.var='factorB')

  tmpa.norm <- sweep(tmpa, 2, colSums(tmpa), FUN='/')
  tmpa.norm.m <- melt(tmpa.norm)

  ggplot(tmpa.norm.m, aes(x=as.factor(Var2), y=value, fill=as.factor(Var1))) + geom_bar(stat='identity') +
    scale_x_discrete(name = factorB.name) + scale_fill_discrete(name = factorA.name) +
    scale_y_continuous(name=paste0('proportion of ',factorA.name))
}





#' Remove levels and elements from a factor that have fewer than min.size entries
#' @param f factor
#' @param min.size minimum size levels to keep
#' @return a new factor
#' @export removeLevelsBySize
removeLevelsBySize <- function(f, min.size = 3) {
    if (!is.factor(f)) stop('f is not a factor!')

    min.size <- 3
    f <- ngrp

    lvl.counts <- table(f)
    cl.remove <- names(lvl.counts)[unname(lvl.counts) < 3]

    removeLevelsFromFactor(f, cl.remove)
}

#' Remove all elements from specificed clusters from a factor
#' and drop the respective levels
#' @param f factor
#' @param cl.remove character vector f clusters to remove
#' @return a new factor
#' @export removeLevelsFromFactor
removeLevelsFromFactor <- function(f, cl.remove) {
    fc <- factor2Char(f)
    nfc <- fc[!fc %in% cl.remove]
    as.factor(nfc)
}




#' Merge a previously broken down factor
#' @param bdf a broken down factor
#' @return a factor
#' @export mergeFactorx
mergeFactor <- function(bdf) {
    ret <- unlist(unname(mapply(function(cells,cl.name) {
        x <- rep(as.character(cl.name), length(cells))
        names(x) <- cells
        x
    },bdf, names(bdf))))
    as.factor(ret)
}
