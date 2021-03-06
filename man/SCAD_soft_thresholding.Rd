% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FSCseq.R
\name{SCAD_soft_thresholding}
\alias{SCAD_soft_thresholding}
\title{SCAD soft thresholding function}
\usage{
SCAD_soft_thresholding(diff_beta, lambda, alpha)
}
\arguments{
\item{diff_beta}{difference in pairwise cluster betas (log2 means)}

\item{lambda}{numeric penalty parameter, lambda >= 0}

\item{alpha}{numeric penalty parameters, 0 <= alpha < 1}
}
\value{
diff_beta value, after SCAD thresholding rule is applied
}
\description{
Takes pairwise distances between cluster log2 means and
penalty parameters lambda and alpha, and imposes
SCAD thresholding rule
}
\examples{
beta1=1
beta2=2
SCAD_soft_thresholding(beta2-beta1,0,0)
SCAD_soft_thresholding(beta2-beta1,0.5,0.45)
SCAD_soft_thresholding(beta2-beta1,1,1)
}
