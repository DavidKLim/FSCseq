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
#' @param method string, either "EM" or "CEM". Default is "CEM"
#' @param n_rinits integer, number of additional random initializations (on top of Hierarchical and K-means) to be searched. Default is 1
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
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
#'
#' @examples
#' sim.dat = FSCseq::simulateData(B=1, g=10000, K=2, n=50, LFCg=1, pDEg=0.05, beta0=12, phi0=0.35, nsims=1, save_file=F)[[1]]
#' \dontrun{FSCseq_results = FSCseq_workflow(cts=sim.dat$cts, K_search=c(2:3), lambda_search=c(1.0, 1.5), alpha_search=c(0.1, 0.2))}
#' \dontrun{summary(FSCseq_workflow$results)}
#'
#' @export
FSCseq_workflow = function(cts, ncores = 1, batch = NULL, X = NULL, true_cls = NULL, true_disc = NULL,
                           method = "CEM", n_rinits = 1, med_filt = 500, MAD_filt = 50, K_search = c(2:6),
                           lambda_search = seq(0.25, 5, 0.25), alpha_search = c(0.01, seq(0.05, 0.5, 0.05)),
                           OS_save = T, tune_save=F, trace = F, trace.prefix = "", nMB = 5, dir_name = "Saved_Results") {
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
    ##more than one batch
    X1 = matrix(nrow = ncol(cts), ncol = length(unique(batch)))
    for (i in 1:length(unique(batch))) {
      X1[, i] = (batch == unique(batch)[i]) ^ 2    # cell-means coding of batch
    }
    # if(min(batch)!=0){batch = batch - min(batch)}
    # X1 = matrix(nrow = ncol(cts), ncol = length(unique(batch))-1)  # reference-cell coding. batch=0 is the reference: to do this, I have to estimate an intercept in M-step!!
    # for(i in 2:(length(unique(batch)))){
    #   X1[, i-1] = (batch == sort(unique(batch))[i])^2
    # }
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
        fname = sprintf("%s/tune%d_%f_%f.out", dir_name, K_search[c],alpha_search[a],lambda_search[l])
        if(file.exists(fname)){
          load(fname)
        }else{
          res = FSCseq::FSCseq(ncores = ncores, X = X, y = cts[idx, ], k = K_search[c],
                             lambda = lambda_search[l], alpha = alpha_search[a],
                             size_factors = SF, norm_y = norm_y[idx, ],
                             true_clusters = true_cls, true_disc = true_disc[idx],
                             init_parms = TRUE, init_coefs = init_coefs, init_phi = init_phi,
                             init_cls = init_cls, init_wts = init_wts, trace = trace,
                             trace.file = trace.file, mb_size = sum(idx))
          if(tune_save){
            save(res,file=fname)
          }
        }

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

  if(OS_save){
    cat("Removing saved interim results...\n")
    for(c in 1:length(K_search)){
      fname = sprintf("%s/OS%d.out", dir_name, K_search[c])
      file.remove(fname)
    }
  }
  if(tune_save){
    cat("Removing saved interim results...\n")
    for(c in 1:length(K_search)){for(a in 1:length(alpha_search)){for(l in 1:length(lambda_search)){
      fname = sprintf("%s/tune%d_%f_%f.out", dir_name, K_search[c],alpha_search[a],lambda_search[l])
      file.remove(fname)
    }}}
  }
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
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
#'
#' @examples
#' sim.dat = simulateData(B=1, g=10000, K=2, n=50, LFCg=1, pDEg=0.05, beta0=12, phi0=0.35, nsims=1, save_file=F)[[1]]
#' \dontrun{FSCseq_results = FSCseq_workflow(cts=sim.dat$cts, K_search=c(2:3), lambda_search=c(1.0, 1.5), alpha_search=c(0.1, 0.2))}
#' \dontrun{proc.dat = FSCseq_results$processed.dat}
#' \dontrun{pred_results = FSCseq_predict_workflow(fit=FSCseq_results$results$fit, cts=sim.dat$cts, cts_pred=sim.dat$cts_pred, idx=proc.dat$idx)}
#'
#' @export
# FSCseq_predict_workflow = function(X = NULL, fit, cts, cts_pred, idx) {
  # add cts_train, batch_train, batch_pred, SF_train into input
  # this function will calculate SF_pred
  # create X_train and X_pred
  # incorporating covariates (train set) and covariates (pred set)
  # need to incorporate option not batch-adjust (in FSCseq.R too: make sure X_train=NULL and X_pred=NULL works)
  # in FSCseq.R: incorporate covariates in X_train and X_pred along with batch

  ############################################ OLD ###################################################

  # geoMeans = exp(rowMeans(log(cts)))   # input custom geometric means (relative to training set)
  # cat("Computing size factors...\n")
  # processed.dat.pred = FSCseq::processData(y = cts_pred, geoMeans = geoMeans,
  #                                          med_filt = FALSE, MAD_filt = FALSE)
  # SF_pred = processed.dat.pred$size_factors
  #
  # cat("Computing predictive posterior probabilities...\n")
  # res_pred = FSCseq_predict(X = NULL, fit = fit, cts_pred = cts_pred[idx, ], SF_pred = SF_pred)
  #
  # return(list(processed.dat.pred = processed.dat.pred, results = res_pred))

  ############################################ NEW ###################################################
FSCseq_predict_workflow = function(res, X_covar_train = NULL, cts_train, SF_train=NULL, batch_train=NULL,
                                   X_covar_pred = NULL, cts_pred, batch_pred=NULL) {
  #### fit = straight from FSCseq_workflow output
  #### X_train and X_pred are covariates
  # idx: ids filtered by rowMedians and MAD
  n_train = ncol(cts_train)
  n_pred = ncol(cts_pred)
  n = n_train + n_pred

  if(!is.null(SF_train)){if(length(SF_train) != n_train){stop("SF_train must be length of ncol(cts_train).")}}
  if(!is.null(batch_train)){if(length(batch_train) != n_train){stop("batch_train must be length ncol(cts_train).")}}
  if(!is.null(batch_pred)){if(length(batch_pred) != n_pred){stop("batch_pred must be length ncol(cts_pred).")}}

  filt_idx = res$processed.dat$idx

  geoMeans = exp(rowMeans(log(cts_train)))
  # Computing size factors of just prediction set (using geoMeans)
  processed.dat.pred = FSCseq::processData(y = cts_pred, geoMeans = geoMeans,
                                           med_filt = FALSE, MAD_filt = FALSE)
  SF_pred = processed.dat.pred$size_factors

  # Design matrix for covariates (non-batch)
  if(is.null(X_covar_train) & is.null(X_covar_pred)) {
    X_covar = NULL
  } else if(!is.null(X_covar_train) & is.null(X_covar_pred)){
    stop("Covariates in training, but no covariates in prediction. Need to input relevant prediction set covariates")
    X_covar = X_covar_train
  } else if(is.null(X_covar_train) & !is.null(X_covar_pred)){
    stop("No covariates in training, but covariates in prediction. Must have adjusted for covariates in training to adjust in prediction")
    X_covar = X_covar_pred
  } else{
    if(ncol(X_covar_train) != ncol(X_covar_pred)){
      stop("Number of columns in X_covar_train and X_covar_pred must be the same. Did you include batch in these design matrices?")
      # users may include batch in these. Put in vignette that batch info will be concatenated into the design matrix, and shouldn't be included here
    }
    X_covar = rbind(X_covar_train, X_covar_pred)
  }
  p_covar = if(!is.null(X_covar)){ncol(X_covar)}else{0}

  # Design matrix for batch
  if(is.null(batch_train) & is.null(batch_pred)){
    # no training batch info, no prediction batch info
    B_train=0; B_pred=0
    cat("No batch info in data (training or prediction). Not adjusting for batch effects in prediction...\n")
    X_batch = NULL
  }else{
    if(is.null(batch_train) & !is.null(batch_pred)){
      # no training batch info, yes prediction batch info
      warning("No training batch info in data, but prediction batch info input.")
      if(length(unique(batch_pred))==1){
        warning("Assuming one training batch and one prediction batch.")
        batch = c(rep(1,n_train),rep(1,n_pred))
      } else{
        warning("Assuming training samples from one separate batch")
        batch_pred = batch_pred - min(batch_pred) + 2   # prediction batch starts from 2, .... B_pred+1
        batch = c(rep(1,n_train), batch_pred)
      }
    }else if(!is.null(batch_train) & is.null(batch_pred)){
      # yes training batch info, no prediction batch info
      warning("No prediction batch info in data. Assuming prediction batch from one separate batch. Was this intended?")
      batch_train = batch_train - min(batch_train) + 1 # training batch starts from 1, ..., B_train
      batch = c(batch_train, rep(length(unique(batch_train))+1,n_pred))
    }else{
      # yes training batch info, yes prediction batch info
      if(any(batch_pred %in% batch_train)){
        warning("Found training batches in prediction batch. Is this intended?")
      }
      # rare case where batch_pred is poorly numbered, i.e.: batch_train=1,...,B and batch_pred = (B+3),....,(B+11)
      if(min(batch_pred) > max(batch_train)+1){ batch_pred = batch_pred - (min(batch_pred) - max(batch_train) - 1) }
      batch = c(batch_train, batch_pred)
    }
    X_batch = matrix(0, nrow=n, ncol=length(unique(batch)))
    for (i in 1:length(unique(batch))) { X_batch[, i] = (batch == unique(batch)[i]) ^ 2 }
  }

  if(is.null(X_batch) & is.null(X_covar)){
    X=NULL
  }else if(!is.null(X_batch) & is.null(X_covar)){
    X=X_batch
  }else if(is.null(X_batch) & !is.null(X_covar)){
    X=X_covar
  }else{
    X=cbind(X_covar, X_batch)
  }

  cat("Computing predictive posterior probabilities...\n")

  #fit=res$results$fit; cts_train=cts; batch_train=batch
  res_pred = FSCseq::FSCseq_predict(X=X, p_covar=p_covar,
                            fit=res$results$fit, cts_train=cts_train[filt_idx,], #batch_train=batch_train, batch_pred=batch_pred,
                            cts_pred=cts_pred[filt_idx,], SF_train=SF_train, SF_pred=SF_pred)
  covariates = list(batch_train = batch_train,
                    batch_pred = batch_pred,
                    X_batch = X_batch,
                    X_covar_train = X_covar_train,
                    X_covar_pred = X_covar_pred)

  res_pred$covariates=covariates
  return(res_pred)
}
