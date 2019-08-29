#' Simulate RNA-seq bulk gene expression count data
#'
#' Simulate data based on input simulation parameters. Size factors are custom input
#' or simulated from N(1,0.25)
#'
#' @param n integer, sample size
#' @param k integer, number of clusters
#' @param g integer, number of genes
#' @param pi numeric vector, mixture proportions, or proportion in each cluster. Equal sizes by default
#' @param coefs numeric matrix of simulated log2 means. Must have g rows and k columns, for each gene and cluster. If NULL, must input beta, p_disc, and LFC
#' @param beta numeric or numeric vector, gene ``baseline" log2 mean. If not vector, then all genes will be simulated with same baseline
#' @param p_disc numeric value from 0 to 1, proportion of genes that are discriminatory (that LFC is applied).
#' @param LFC numeric value greater than 0, log2 fold change to be applied on discriminatory genes. For each gene, one random cluster simulated to be either up-regulated (add LFC) or down-regulated (subtract LFC)
#' @param disp string, either "gene" or "cluster". "gene" simulates gene-level overdispersion parameters, and "cluster" simulates cluster-level overdispersion parameters
#' @param phi numeric, overdispersion parameters. must be vector of length g if disp=='gene', or matrix of g rows and k columns if disp=='cluster'
#' @param size_factors numeric vector of simulated size factors for each subject. By default, simulated from N(1,0.25)
#' @param batch_effects numeric vector of effects corresponding to each level of 'batch'. Must have length equal to the number of batches specified. Effects are applied in order (i.e. first element of batch_effects refers to the effect of batch 1, etc)
#' @param batch integer vector of simulated batch of each sample. Must be numeric 1 to (number of batches).
#'
#' @return list containing the following objects:
#' y: count matrix of g by n
#' clusters: cluster indices for each of the n samples
#' batch: input known batch of each sample (all 1 by default)
#' batch_effects: input batch effects for each batch (0 by default)
#' batch_g: genes on which batch effects were applied (on random half of genes)
#' size_factors: simulated or input size factors (by default, simulated from N(1,0.25))
#' coefs: matrix of log2 means for each gene and cluster, either input matrix or constructed based on beta, p_disc, and LFC
#' phi: numeric overdispersion parameters, either vector or matrix based on dispersion scheme
#' true_disc: logical vector of length g, TRUE if gene is simulated discriminatory (not simulated by same beta across clusters), FALSE if not
#'
#' @export
simulateData <- function(n,k,g,pi=rep(1/k,k),
                         coefs=NULL,beta,p_disc,LFC,
                         disp="gene",phi,
                         size_factors=NULL,
                         batch_effects=0,batch=rep(1,n)){

  if(length(pi)!=k){stop("pi must have length k")}
  if(any(pi<=0) | any(pi>1)){stop("pi must all be between 0 and 1")}
  if(!is.null(coefs)){
    if(!is.matrix(coefs)){stop("coefs must be a matrix")}
    if(nrow(coefs) != g | ncol(coefs) != k){stop("coefs must have g rows and k columns")}
  }
  if(disp=="gene"){
    if(!(length(phi) %in% c(1,g))){stop("phi must have length 1 or g")}
  } else if(disp=="cluster"){
    if(!is.matrix(phi)){stop("phi must be a matrix")}
    if(nrow(phi) != g | ncol(phi) != k){stop("phi must have g rows and k columns")}
  } else{ stop("disp must be 'gene' or 'cluster'") }

  # simulate size factors from N(1,0.25) if custom size factors are not input
  if(is.null(size_factors)){
    size_factors = simulateSizeFactors(n)
  }else{
    if(length(size_factors)!=n){stop("input size_factors must have length n")}
    if(any(size_factors<=0)){stop("input size_factors must be >0 (1=no adjustment)")}
  }
  if(length(batch)!=n){stop("batch must be specified for all n samples")}
  if(length(batch_effects) != length(unique(batch))){stop("Number of batch effects, and number of batches specified do not match.")}


  if(is.null(coefs)){
    message("Custom coefs matrix not specified. Constructing matrix based on input beta, p_disc, and LFC")
    if(!(length(beta) %in% c(1,g))){stop("beta must be a vector of length 1 (all genes same value) or g")}
    if(p_disc < 0 | p_disc > 1){stop("p_disc must be between 0 and 1")}
    coefs=matrix(beta,nrow=g,ncol=k)                      # initialize coefs matrix
    disc_ids=sample(1:g,ceiling(p_disc*g),replace=F)      # sample appropriate number of genes to be disc
    for(j in disc_ids){
      LFC_applied=(rbinom(1,1,0.5)-0.5)*2*LFC             # Random up/down regulation by LFC
      coefs[j,sample(1:k,1)] =+ LFC_applied     # LFC applied to random cluster
    }
  }

  true_disc=apply(coefs,1,function(x) {(max(x)-min(x))!=0})     # sees which genes have range >0

  batch_eff = batch_effects[batch]    # batch effect applied for each sample
  y<-matrix(rep(0,times=g*n),nrow=g)  # initialize count matrix
  z = rmultinom(n,1,pi)          # sample cluster membership from multinomial of mixture proportions
  batch_g = sample(1:g, 0.5*g)   # Apply batch on 50% of genes

  cl = rep(0,n)
  for(c in 1:k){
    cl[z[c,]==1] = c
  }
  if(length(phi)==1){phi=rep(phi,g)}

  if(disp=="gene"){
    for(j in 1:g){
      for(i in 1:n){
        if(j %in% batch_g){
          y[j,i] = rnbinom( 1, size = 1/phi[j], mu = size_factors[i]*2^(coefs[j,cl[i]] + batch_eff[i]))
        } else{
          y[j,i] = rnbinom( 1, size = 1/phi[j], mu = size_factors[i]*2^(coefs[j,cl[i]]))
        }
      }
    }
  } else if(disp=="cluster"){
    for(j in 1:g){
      for(i in 1:n){
        if(j %in% batch_g){
          y[j,i] = rnbinom( 1, size = 1/phi[j,cl[i]], mu = size_factors[i]*2^(coefs[j,cl[i]] + batch_eff[i]))
        } else{
          y[j,i] = rnbinom( 1, size = 1/phi[j,cl[i]], mu = size_factors[i]*2^(coefs[j,cl[i]]))
        }
      }
    }
  }

  result<-list(y=y,
               clusters=cl,
               batch=batch,
               batch_effects=batch_effects,
               batch_g=batch_g,       # genes on which batch effect was applied
               size_factors=size_factors,
               coefs=coefs,
               phi=phi,
               true_disc=true_disc)
  return(result)
}

#' Simulate size factors from N(1,0.25)
simulateSizeFactors <- function(n){
  size_factors=rnorm(n,mean=1,sd=0.25)
  size_factors[size_factors<0]=abs(size_factors[size_factors<0])    # replace negative value with its absolute value
  size_factors[size_factors==0]=1        # if size_factors somehow =0, then set it to 1 (no adjustment)
  return(size_factors)
}



