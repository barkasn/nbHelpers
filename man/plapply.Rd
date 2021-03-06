% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{plapply}
\alias{plapply}
\title{Parallel apply using the packages parallel, BiocParallel or regular lapply}
\usage{
plapply(..., n.cores = detectCores(logical = FALSE),
  mc.preschedule = FALSE)
}
\arguments{
\item{...}{parameters to the apply functions, do not include BPPARAM or mc.cores}

\item{n.cores}{numbere of cores to use, will automatically detect the number of physical cores if possible}

\item{mc.preschedule}{parameter for mclapply}
}
\description{
From the pagoda package with some minor improvements originally called papply
}
