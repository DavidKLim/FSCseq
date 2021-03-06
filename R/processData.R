#' Pre-process the count matrix
#'
#' Estimate size factors using DESeq2, and pre-filter based on low count and
#' low variability genes using custom thresholds on median normalized count and
#' median absolute deviation (MAD) values.
#'
#' @param y integer, gene expression count matrix (output from simulateData)
#' @param geoMeans (optional) numeric vector of length g (number of rows of y): custom geometric means of gene counts. Used to estimate prediction set size factors based on training set
#' @param estimateSFtype string, input for estimateSizeFactors "type" argument. must be 'ratio', 'poscounts', or 'iterate'. See DESeq2 vignette.
#' @param med_filt logical, TRUE to filter low-count genes via median threshold
#' @param MAD_filt logical, TRUE to filter low-variable genes via MAD quantile threshold
#' @param med_thresh numeric, median threshold for pre-filtering low-count genes (default 100, i.e. pre-filters genes with median normalized count below 100)
#' @param MAD_quant_thresh numeric value between 0 to 100, quantile threshold imposed on MAD values for pre-filtering low-variable genes (default 50, i.e. pre-filters genes below 50th quantile of MAD values)
#'
#' @return list containing the following objects:
#' dds: DESeq2 output.
#' size_factors: numeric vector, estimated size factors, from DESeq2
#' norm_y: numeric normalized count matrix, corrected for differences in sequencing depth, from DESeq2
#' idx: logical vector, TRUE for inclusion after pre-filtering low-count and low-variable genes
#' row_medians: numeric vector, median normalized count for each gene
#' row_MADs: numeric vector, median absolute deviation (MAD) value of log(norm_y+0.1) for each gene
#'
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
#'
#' @examples
#' sim.dat = simulateData(B=1, g=10000, K=2, n=50, LFCg=1, pDEg=0.05, beta0=12, phi0=0.35, nsims=1, save_file=F)[[1]]
#' proc.dat = processData(sim.dat$cts)
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix varianceStabilizingTransformation rlogTransformation estimateSizeFactors sizeFactors counts
#' @importFrom matrixStats rowMedians
#'
#' @export
processData = function(y, geoMeans = NULL, estimateSFtype = "ratio",
                       med_filt = TRUE, MAD_filt = TRUE,
                       med_thresh = 100, MAD_quant_thresh = 50) {
  n = ncol(y)
  g = nrow(y)
  coldata <- data.frame(matrix(rep(1, n), nrow = n))
  rownames(coldata) <- colnames(y)
  colnames(coldata) <- "int_only"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = y,
                                        colData = coldata,
                                        design = ~ 1)

  if (is.null(geoMeans)) {
    dds = DESeq2::estimateSizeFactors(dds, type = estimateSFtype)
  } else{
    dds = DESeq2::estimateSizeFactors(dds, type = estimateSFtype, geoMeans =
                                        geoMeans)
  }


  size_factors <- DESeq2::sizeFactors(dds)
  norm_y <- DESeq2::counts(dds, normalized = TRUE)

  idx = rep(T, g)
  row_medians = matrixStats::rowMedians(norm_y)
  row_MADs = apply(log(norm_y + 0.1), 1, mad)

  if (med_filt) {
    idx = (idx & (row_medians >= med_thresh))
  }
  if (MAD_filt) {
    idx = (idx & (row_MADs >= quantile(row_MADs, MAD_quant_thresh / 100)))
  }
  results = list(
    dds = dds,
    size_factors = size_factors,
    norm_y = norm_y,
    idx = idx,
    row_medians = row_medians,
    row_MADs = row_MADs
  )

  return(results)
}



#' Summarize FSCseq clustering results
#'
#' Outputs relevant summary of FSCseq clustering results
#'
#' @param res Output of FSCseq clustering analysis
#' @param true_cls (optional) vector of true cluster labels to calculate Adjusted Rand Index (ARI)
#' @param true_disc (optional) TRUE/FALSE vector of true cluster-discriminatory gene status to calculate True Positive Rate (TPR) and False Positive Rate (FPR)
#'
#' @return Summary of FSCseq clustering results: #clusters (K), clusters,
#' ARI (if true_cls input), TPR and FPR (vs. true cluster-discriminatory genes if true_disc input)
#'
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
#'
#' @examples
#' sim.dat = FSCseq::simulateData(B=1, g=10000, K=2, n=50, LFCg=1, pDEg=0.05, beta0=12, phi0=0.35, nsims=1, save_file=F)[[1]]
#' \dontrun{FSCseq_results = FSCseq_workflow(cts=sim.dat$cts, K_search=c(2:3), lambda_search=c(1.0, 1.5), alpha_search=c(0.1, 0.2))}
#' \dontrun{summary(FSCseq_workflow$results)}
#'
#' @importFrom mclust adjustedRandIndex
#'
#' @export
summary = function(res, true_cls = NULL, true_disc = NULL) {
  # extract relevant quantities from FSCseq results
  K=res$K; cls=res$cls; disc = res$discriminatory

  # Output K
  cat(paste("K:", K,"\n"))

  # if true_cls input, compare cls to true_cls:
  ARI = NA
  if(!is.null(true_cls)){
    cat(paste( "True K:", length(unique(true_cls)),"\n" ))
    cat(paste( "ARI:", adjustedRandIndex(true_cls,cls),"\n" ))
    print(table(true_cls, cls))
  }
  cat("-----------------\n")

  # if true_disc input, compare disc to true_disc
  TPR = NA; FPR = NA
  if(!is.null(true_disc)){
    cat(paste( "TPR:", sum(true_disc & disc)/sum(true_disc), "\n" ))
    cat(paste( "FPR:", sum(!true_disc & disc)/sum(!true_disc), "\n" ))
    cat("-----------------\n")
  }

  # Output first 5 cluster labels, and first 5 genes' cluster-discriminatory status
  cat("cls (first 5 samples):\n")
  print(head(cls, n=5))
  cat("disc (first 5 genes):\n")
  print(head(disc, n=5))

  return(list(K=K, cls=cls, disc=disc,
              ARI=ARI, TPR=TPR, FPR=FPR))
}
