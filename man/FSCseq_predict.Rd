% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FSCseq.R
\name{FSCseq_predict}
\alias{FSCseq_predict}
\title{Prediction via FSCseq}
\usage{
FSCseq_predict(
  X = NULL,
  fit,
  cts_train = NULL,
  cts_pred,
  SF_train = NULL,
  SF_pred,
  maxit = 100,
  eps = 1e-04
)
}
\arguments{
\item{X}{(optional) design matrix of dimension n by p}

\item{fit}{FSCseq results object}

\item{cts_train}{Training count data matrix of dimension g by n_train}

\item{cts_pred}{Prediction/test count data matrix of dimension g by n_pred}

\item{SF_train}{Vector of length n_train, size factors for training subjects. Can be accessed from the fit.}

\item{SF_pred}{Vector of length n_pred, size factors for prediction subjects}

\item{maxit}{Maximum number of iterations to run if prediction batches are estimated, default is 100.}

\item{eps}{Tolerance for relative change in Q function convergence criterion if prediction batches are estimated, default is 1e-4.}
}
\value{
list containing outputs
final_clusters: vector of length n of resulting clusters,
wts: k by n matrix of E step weights
}
\description{
Performs prediction of cluster membership via trained model
}
\references{
\url{https://github.com/DavidKLim/FSCseq}
}
\author{
David K. Lim, \email{deelim@live.unc.edu}
}
