#' Minimal workflow for FSCseq
#'
#' Full FSCseq workflow based on minimal working defaults
#'
#' @param cts integer matrix, count matrix of dimension g by n. Must be integers (counts)
#' @param ncores integer, number of cores (for parallel computing). Default is 1
#' @param batch vector of batch, to use as covariates. Default is one batch (NULL).
#' @param X optional input design matrix to specify p arbitrary covariates/confounders. Must be matrix of dimension n x p. If batch and X are both specified, then X is augmented to incorporate batch as covariates.
#' @param true_cls (optional) integer vector of true groups, if available, for diagnostic tracking.
#' @param true_disc (optional) logical vector of true discriminatory genes, if available, for diagnostic tracking.
#' @param method string, either "EM" or "CEM". Default is "EM"
#' @param n_rinits integer, number of additional random initializations (on top of Hierarchical and K-means) to be searched. Default is 20
#' @param med_filt integer, threshold for minimum median gene normalized count for pre-filtering. med_filt=0 pre-filters no genes via this criterion. Default is 500.
#' @param MAD_filt integer, value between 0 and 100. quantile threshold for gene log MAD of normalized count. MAD_filt=0 pre-filters no genes via this criterion. Default is 50.
#' @param K_search integer vector, values of K (number of clusters) to be searched. Default is 2:6
#' @param lambda_search numeric vector, values of lambda to be searched. Default is seq(0.25,3,0.25)
#' @param alpha_search numeric vector, values of alpha to be searched. Default is c(0.01,seq(0.05,0.50,0.05))
#' @param OS_save logical, TRUE: saves progress of computationally costly warm starts (multiple initializations). Default is TRUE
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
#' @param trace.prefix (optional) string, prefix of file name to store trace output.
#' @param nMB integer, number of minibatches to use in M step. Default is 5
#' @param dir_name string, name of directory specified for saved results (if OS_save = TRUE) and diagnostics (if trace = TRUE)
#'
#' @return list with K, cls, discriminatory, and fit
#'
#' @export
FSCseq_workflow = function(cts, ncores = 1, batch = NULL, X = NULL, true_cls = NULL, true_disc = NULL,
                           method = "CEM", n_rinits = 1, med_filt = 500, MAD_filt = 50, K_search = c(2:6),
                           lambda_search = seq(0.25, 5, 0.25), alpha_search = c(0.01, seq(0.05, 0.5, 0.05)),
                           OS_save = TRUE, trace = F, trace.prefix = "", nMB = 5, dir_name = "Saved_Results") {
  ifelse(!dir.exists(dir_name),
         dir.create(dir_name, recursive = T),
         FALSE)
  if (trace) {
    ifelse(!dir.exists(sprintf("%s/Diagnostics", dir_name)),
           dir.create(sprintf("%s/Diagnostics", dir_name), recursive = T),
           FALSE)
  }
  # input true_cls or true_disc for clustering/FS diagnostics tracking

  if (!trace) {
    trace.file = NULL
  }

  if (is.null(batch)) {
    cat("No input batch. All samples from the same batch\n")
    batch = rep(1, ncol(cts))
  }

  # X1: Design matrix for batch
  if (length(unique(batch)) != 1) {
    # more than one batch
    X1 = matrix(nrow = ncol(cts), ncol = length(unique(batch)))
    for (i in 1:length(unique(batch))) {
      X1[, i] = (batch == unique(batch)[i]) ^ 2    # cell-means coding of batch
    }
  } else{
    # all the same batch: no need to adjust
    X1 = NULL
  }

  if(!is.null(X)){
    X = ifelse(is.null(X1), X, cbind(X,X1))
  } else{
    X = X1
  }

  cat("Computing size factors...\n")
  processed.cts = FSCseq::processData(cts, med_thresh = med_filt, MAD_quant_thresh = MAD_filt)
  norm_y = processed.cts$norm_y
  SF = processed.cts$size_factors
  idx = processed.cts$idx
  processed.dat = processed.cts[-1] # remove dds object (large)
  rm("processed.cts")

  mb_size = floor(sum(idx) / nMB)

  # warm starts for each value of K
  cat("Initializing warm starts...\n")
  list_res = list()
  for (c in 1:length(K_search)) {
    cat(paste("K =",K_search[c],"... "))
    fname = sprintf("%s/OS%d.out", dir_name, K_search[c])
    if (trace) {
      trace.file = sprintf("%s/Diagnostics/OS%d.txt", dir_name, K_search[c])
    }

    if (file.exists(fname)) {
      load(fname)
      warning(
        "Previous initialization results found in directory. Loading these initializations. Note: an error may occur if previous results correspond to a different dataset!"
      )
      list_res[[c]] = res
    } else{
      res = FSCseq::FSCseq(ncores = ncores, X = X, y = cts[idx, ], k = K_search[c],
                           lambda = 0.05, alpha = 0.01, size_factors = SF, norm_y = norm_y[idx, ],
                           true_clusters = true_cls, true_disc = true_disc[idx], init_parms = FALSE,
                           init_coefs = NULL, init_phi = NULL, init_cls = NULL, init_wts = NULL,
                           n_rinits = n_rinits, method = method, trace = trace, trace.file = trace.file,
                           mb_size = mb_size)
      list_res[[c]] = res
      if (OS_save) {
        save(res, file = fname)
      }
    }
    cat("done.\n")
  }

  # tuning (K, lambda, alpha) jointly
  n_tune = length(K_search) * length(alpha_search) * length(lambda_search)
  BICs = matrix(nrow = n_tune, ncol = 4)
  index = 1
  list_res_tune = list()

  cat("Tuning parameters...\n")
  for (c in 1:length(K_search)) {
    cat(paste("K =",K_search[c],"... "))
    for (a in 1:length(alpha_search)) {
      for (l in 1:length(lambda_search)) {
        if (trace) {
          trace.file = sprintf("%s/Diagnostics/JOINT%d_%f_%f.txt",
                               dir_name, K_search[c], lambda_search[l], alpha_search[a])
        }

        if (l == 1) {
          init_coefs = list_res[[c]]$coefs
          init_phi = list_res[[c]]$phi
          init_cls = list_res[[c]]$clusters
          init_wts = list_res[[c]]$wts
        } else if (l > 1) {
          init_coefs = list_res_tune[[index - 1]]$coefs
          init_phi = list_res_tune[[index - 1]]$phi
          init_cls = list_res_tune[[index - 1]]$clusters
          init_wts = list_res_tune[[index - 1]]$wts
        }
        res = FSCseq::FSCseq(ncores = ncores, X = X, y = cts[idx, ], k = K_search[c],
                             lambda = lambda_search[l], alpha = alpha_search[a],
                             size_factors = SF, norm_y = norm_y[idx, ],
                             true_clusters = true_cls, true_disc = true_disc[idx],
                             init_parms = TRUE, init_coefs = init_coefs, init_phi = init_phi,
                             init_cls = init_cls, init_wts = init_wts, trace = trace,
                             trace.file = trace.file, mb_size = sum(idx))

        list_res_tune[[index]] = res
        BICs[index, ] = c(K_search[c], alpha_search[a], lambda_search[l], res$BIC)
        if (trace) {
          print(BICs[index, ])
        }
        index = index + 1
      }
    }
    cat("done.\n")
  }

  optim_id = which.min(BICs[, 4])[1]   # if more than one result with minimum BIC, select first one
  optim_res = list_res_tune[[optim_id]]

  K = length(unique(optim_res$clusters))
  cls = optim_res$clusters
  discriminatory = optim_res$discriminatory

  results = list(
    K = K,
    cls = cls,
    discriminatory = discriminatory,
    fit = optim_res
  )

  return(list(processed.dat = processed.dat, results = results))

}

#' Minimal workflow for FSCseq_predict
#'
#' Full FSCseq workflow based on minimal working defaults
#'
#' @param X covariates (optional)
#' @param fit FSCseq results object. Accessed by $results from FSCseq_workflow object
#' @param cts Training set count matrix, dimension g by n. Used as pseudo-reference to calculate prediction set size factors
#' @param cts_pred integer matrix, count matrix of dimension g by n_pred. Must be integers (counts)
#' @param idx boolean vector: TRUE if gene passed pre-filtering step
#'
#' @return list with processed.dat.pred (processed prediction data), and prediction results
#'
#' @export
FSCseq_predict_workflow = function(X = NULL, fit, cts, cts_pred, idx) {
  geoMeans = exp(rowMeans(log(cts)))   # input custom geometric means (relative to training set)
  cat("Computing size factors...\n")
  processed.dat.pred = FSCseq::processData(y = cts_pred, geoMeans = geoMeans,
                                           med_filt = FALSE, MAD_filt = FALSE)
  SF_pred = processed.dat.pred$size_factors

  cat("Computing predictive posterior probabilities...\n")
  res_pred = FSCseq_predict(X = NULL, fit = fit, cts_pred = cts_pred[idx, ], SF_pred = SF_pred)

  return(list(processed.dat.pred = processed.dat.pred, results = res_pred))
}
