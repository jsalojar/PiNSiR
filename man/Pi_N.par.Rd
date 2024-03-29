% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/picalc.R
\name{Pi_N.par}
\alias{Pi_N.par}
\title{Calculate Pi_N for loci with a high impact or moderate impact mutation using parallel computation}
\usage{
Pi_N.par(TH, exon.sel, exon.mut, nSNP = NA, datatype = "full", cores = 15)
}
\arguments{
\item{TH}{Site-wise theta calculated using ANGSD}

\item{exon.sel}{Selected gene models}

\item{exon.mut}{Coordinates of the non-synonymous positions.}

\item{nSNP}{If given, divides the total Pi with the given number of SNPs.}

\item{datatype}{"full" assumes that all positions are called Pi is a mean of observed values, "SNP" assumes only variant positions and the Pi is divided by the total number of non-synonymous positions in the selected gene models.}

\item{cores}{Number of cores for parallel computation}
}
\value{
a list with Pi, chromosome-wise cumulative sums and numbers of SNPs.
}
\description{
Calculate Pi_N for loci with a high impact or moderate impact mutation using parallel computation
}
