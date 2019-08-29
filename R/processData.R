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
#' @param calc_vsd TRUE if want to output count after variance-stabilizing transform
#' @param calc_rld TRUE if want to output count after rlog transform
#'
#' @return list containing the following objects:
#' dds: DESeq2 output.
#' size_factors: numeric vector, estimated size factors, from DESeq2
#' norm_y: numeric normalized count matrix, corrected for differences in sequencing depth, from DESeq2
#' idx: logical vector, TRUE for inclusion after pre-filtering low-count and low-variable genes
#' row_medians: numeric vector, median normalized count for each gene
#' row_MADs: numeric vector, median absolute deviation (MAD) value of log(norm_y+0.1) for each gene
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix varianceStabilizingTransformation rlogTransformation estimateSizeFactors sizeFactors counts
#' @importFrom matrixStats rowMedians
#'
#' @export
processData = function(y,geoMeans=NULL,estimateSFtype="ratio",
                       med_filt=TRUE, MAD_filt=TRUE,
                       med_thresh=100, MAD_quant_thresh=50){
  n=ncol(y)
  g=nrow(y)
  coldata<-data.frame(matrix(rep(1,n),nrow=n))
  rownames(coldata)<-colnames(y)
  colnames(coldata)<-"int_only"
  dds<-DESeq2::DESeqDataSetFromMatrix(countData = y,
                                      colData = coldata,
                                      design = ~ 1)

  if(is.null(geoMeans)){
    dds = DESeq2::estimateSizeFactors(dds,type=estimateSFtype)
  } else{
    dds = DESeq2::estimateSizeFactors(dds,type=estimateSFtype,geoMeans=geoMeans)
  }


  size_factors<-DESeq2::sizeFactors(dds)
  norm_y<-DESeq2::counts(dds,normalized=TRUE)

  idx=rep(T,g)
  row_medians=matrixStats::rowMedians(norm_y)
  row_MADs=apply(log(norm_y+0.1),1,mad)

  if(med_filt){
    idx = ( idx & (row_medians >= med_thresh) )
  }
  if(MAD_filt){
    idx = ( idx & (row_MADs >= quantile(row_MADs,MAD_quant_thresh/100)) )
  }
  results=list(dds=dds, size_factors=size_factors,norm_y=norm_y,
               idx=idx, row_medians=row_medians, row_MADs=row_MADs)

  return(results)
}
