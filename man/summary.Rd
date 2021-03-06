% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\name{summary}
\alias{summary}
\title{Summarize FSCseq clustering results}
\usage{
summary(res, true_cls = NULL, true_disc = NULL)
}
\arguments{
\item{res}{Output of FSCseq clustering analysis}

\item{true_cls}{(optional) vector of true cluster labels to calculate Adjusted Rand Index (ARI)}

\item{true_disc}{(optional) TRUE/FALSE vector of true cluster-discriminatory gene status to calculate True Positive Rate (TPR) and False Positive Rate (FPR)}
}
\value{
Summary of FSCseq clustering results: #clusters (K), clusters,
ARI (if true_cls input), TPR and FPR (vs. true cluster-discriminatory genes if true_disc input)
}
\description{
Outputs relevant summary of FSCseq clustering results
}
\examples{
sim.dat = FSCseq::simulateData(B=1, g=10000, K=2, n=50, LFCg=1, pDEg=0.05, beta0=12, phi0=0.35, nsims=1, save_file=F)[[1]]
\dontrun{FSCseq_results = FSCseq_workflow(cts=sim.dat$cts, K_search=c(2:3), lambda_search=c(1.0, 1.5), alpha_search=c(0.1, 0.2))}
\dontrun{summary(FSCseq_workflow$results)}

}
\references{
\url{https://github.com/DavidKLim/FSCseq}
}
\author{
David K. Lim, \email{deelim@live.unc.edu}
}
