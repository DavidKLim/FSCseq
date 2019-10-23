#' Minimal workflow for FSCseq
#'
#' Full FSCseq workflow based on minimal working defaults
#'
#' @param cts integer matrix, count matrix of dimension g by n. Must be integers (counts)
#' @param ncores integer, number of cores (for parallel computing). Default is 1
#' @param batch vector of batch, to use as covariates. Default is no batches.
#' @param true_cls (optional) integer vector of true groups, if available, for diagnostic tracking.
#' @param true_disc (optional) logical vector of true discriminatory genes, if available, for diagnostic tracking.
#' @param method string, either "EM" or "CEM". Default is "EM"
#' @param n_rinits integer, number of additional random initializations (on top of Hierarchical and K-means) to be searched. Default is 20
#' @param med_filt numeric, threshold for minimum median gene normalized count for pre-filtering. med_filt=0 pre-filters no genes via this criterion. Default is 500.
#' @param MAD_filt numeric, value between 0 and 100. quantile threshold for gene log MAD of normalized count. MAD_filt=0 pre-filters no genes via this criterion. Default is 50.
#' @param K_search integer vector, values of K (number of clusters) to be searched. Default is 2:6
#' @param lambda_search numeric vector, values of lambda to be searched. Default is seq(0.25,3,0.25)
#' @param alpha_search numeric vector, values of alpha to be searched. Default is c(0.01,seq(0.05,0.50,0.05))
#' @param OS_save logical, TRUE: saves progress of computationally costly warm starts (multiple initializations). Default is TRUE
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
#' @param trace.file (optional) string, file name to store trace output. For example: 'test.txt'
#' @param nMB integer, number of minibatches to use in M step. Default is 5
#'
#' @return list with K, cls, discriminatory, and fit
#'
#' @export
FSCseq_workflow = function(cts,ncores=1,batch=rep(1,ncol(cts)),true_cls=NULL,true_disc=NULL,
                           method="EM",n_rinits=20,med_filt=500,MAD_filt=50,
                           K_search=c(2:6),lambda_search=seq(0.25,3,0.25),alpha_search=c(0.01,seq(0.05,0.5,0.05)),
                           OS_save=TRUE,trace=F,trace.file=NULL,nMB=5){

  dir_name="Saved_Results"
  if(OS_save){
    ifelse(!dir.exists(dir_name),
           dir.create(dir_name),
           FALSE)
  }
  # input true_cls or true_disc for clustering/FS diagnostics tracking

  if(length(unique(batch))==1){
    X=NULL
  } else{
    X=matrix(nrow=ncol(cts),ncol=length(unique(batch)))
    for(i in 1:length(unique(batch))){
      X[,i] = (batch==unique(batch)[i])^2    # cell-means coding of batch
    }
  }
  mb_size = floor(nrow(cts)/nMB)

  processed.cts=FSCseq::processData(cts,med_thresh=med_filt,MAD_quant_thresh=MAD_filt)
  norm_y=processed.cts$norm_y
  SF=processed.cts$size_factors
  idx=processed.cts$idx

  # warm starts for each value of K
  list_res=list()
  for(c in 1:length(K_search)){
    fname=sprintf("%s/OS%d.out",dir_name,K_search[c])
    if(file.exists(fname)){
      load(fname)
      list_res[[c]]=res
    } else{
      res=FSCseq::FSCseq(ncores=ncores,X=X, y=cts[idx,], k=K_search[c],
                                     lambda=lambda_search[1],alpha=alpha_search[1],
                                     size_factors=SF,norm_y=norm_y[idx,],
                                     true_clusters=true_cls, true_disc=true_disc,
                                     init_parms=FALSE,init_coefs=NULL,init_phi=NULL,init_cls=NULL,
                                     n_rinits=n_rinits,method=method,
                                     trace=trace,trace.file=trace.file,
                                     mb_size=mb_size)
      list_res[[c]]=res
      if(OS_save){
        save(res,file=sprintf("%s/OS%d.out",dir_name,K_search[c]))
      }
    }
  }

  # tuning (K, lambda, alpha) jointly
  n_tune=length(K_search)*length(alpha_search)*length(lambda_search)
  BICs=matrix(nrow=ntune,ncol=4)
  index=1
  list_res_tune = list()
  for(c in 1:length(K_search)){for(a in 1:length(alpha_search)){for(l in 1:length(lambda_search)){
    if(a==1 & l==1){
      res = list_res[[c]]
    } else {
      if(a>1 & l==1){
        init_coefs=list_res[[c]]$coefs; init_phi=list_res[[c]]$phi; init_cls=list_res[[c]]$clusters
      } else if(a>=1 & l>1){
        init_coefs=list_res_tune[[index-1]]$coefs; init_phi=list_res_tune[[index-1]]$phi; init_cls=list_res_tune[[index-1]]$clusters
      }
      res = FSCseq::FSCseq(ncores=ncores,X=X, y=cts[idx,], k=K_search[c],
                         lambda=lambda_search[l],alpha=alpha_search[a],
                         size_factors=SF,norm_y=norm_y[idx,],
                         true_clusters=true_cls, true_disc=true_disc,
                         init_parms=TRUE,init_coefs=init_coefs,init_phi=init_phi,init_cls=init_cls,
                         trace=trace,trace.file=trace.file,
                         mb_size=sum(idx))   # minibatching disabled after warm start
    }

    list_res_tune[[index]]=res
    BICs[index,]=c(K_search[c],alpha_search[a],lambda_search[l],res$BIC)
    index=index+1
  }}}

  optim_id = which.min(BICs[,4])[1]   # if more than one result with minimum BIC, select first one
  optim_res = list_res_tune[[optim_id]]

  K=length(unique(optim_res$clusters))
  cls=optim_res$clusters
  discrimiantory=optim_res$discriminatory

  results=c(K=K,
            cls=cls,
            discriminatory=discriminatory,
            fit=optim_res)

  return(results)

}
