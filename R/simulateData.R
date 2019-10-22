#' Simulate data with gene-wise dispersion parameters
simulate_counts=function(K,B,n,g,
                         cls,SF,
                         beta,phi,LFCg_mat,#noise_mat,
                         batch,LFCb,sigma_b,sigma_g,DEb_ID,disp){     # batch of 0: no batch effect, batch of 1: yes effect
  # center batch at 0: lower batch gets downregulated, higher batch gets upregulated
  # if B=1, batch_eff = 0 with some noise by sigma_b
  batch_eff = (batch-(B-1)/2)*LFCb+rnorm(n,0,sigma_b)
  cts <- matrix(rep(0,times=g*n),nrow=g)

  b = beta + LFCg_mat # + noise_mat

  if(disp=="gene"){
    phi_mat = matrix(phi,nrow=g,ncol=K)      # construct matrix with same phi values for all cls if g-disp
  }
  for(j in 1:g){
    for(i in 1:n){
      if(DEb_ID[j]){      # if this is a gene that is differentially expressed due to batch
        cts[j,i] = rnbinom( 1, size = 1/phi_mat[j,cls[i]], mu = SF[i]*2^(b[j,cls[i]] + rnorm(1,0,sigma_g) + batch_eff[i]))
      } else{
        cts[j,i] = rnbinom( 1, size = 1/phi_mat[j,cls[i]], mu = SF[i]*2^(b[j,cls[i]] + rnorm(1,0,sigma_g)))
      }
    }
  }

  return(cts)
}

#' Simulate RNA-seq bulk gene expression count data
#'
#' Simulate data based on input simulation parameters. Size factors are custom input
#' or simulated from N(1,0.25)
#'
#' @param K integer, number of clusters
#' @param B integer, number of batches
#' @param g integer, number of genes
#' @param n integer, number of samples
#' @param pK vector of length K (optional): proportion of samples in each cluster
#' @param pB vector of length B (optional): proportion of samples in each batch
#' @param LFCg numeric, LFC for cluster-discriminatory genes
#' @param pDEg numeric, proportion of genes that are cluster-discriminatory
#' @param sigma_g numeric, Gaussian noise added from N(0,sigma_g). Default is 0.1
#' @param LFCb numeric, LFC for genes that are differentially expressed across batch. Default is 1.
#' @param pDEb numeric, proportion of genes that are differentially expressed across batch. Default is 0.5.
#' @param sigma_b numeric, Gaussian noise added to each batch (turned off, set to 0).
#' @param beta0 numeric, baseline log2 expression for each gene before LFC is applied
#' @param phi0 numeric, baseline overdispersion for each gene
#' @param SF vector of length n (optional), custom size factors from DESeq2. If NULL, simulated from N(1,0.25)
#' @param nsims integer, number of datasets to simulate given the input conditions. Default is 25.
#' @param disp string, either 'gene' or 'cluster' to simulate gene-level or cluster-level dispersions. Default is gene-level. Input phi must be g x K matrix if disp='cluster'
#' @param n_pred integer, number of samples in simulated prediction dataset. Default is 25
#' @param save_dir string (optional): directory to save files. Default: 'Simulations/<sigma_g>_<sigma_b>/B<B>'
#' @param save_file string (optional): prefix of file name to save simulated data to. Default: '<K>_<n>_<LFCg>_<pDEg>_<beta0>_<phi0>'
#'
#' @return saved file in '<save_dir>/<save_file>_sim<1:nsims>_data.RData'
#'
#' @export
simulateData<-function(K, B=1, g=10000, n, pK=NULL, pB=NULL,
                        LFCg, pDEg, sigma_g=0.1,
                        LFCb=1, pDEb=0.5, sigma_b=0,
                        beta0, phi0, SF=NULL,
                        nsims=25, disp="gene", n_pred=25, save_dir=NULL,save_file=NULL){

  if(is.null(save_dir)){
    dir_name=sprintf("Simulations/%f_%f/B%d",sigma_g,sigma_b,B)
  }else{dir_name=save_dir}
  ifelse(!dir.exists(dir_name),
         dir.create(dir_name),
         FALSE)

  if(is.null(save_file)){
    file_name=sprintf("%s/%d_%d_%f_%f_%f_%f",
                      dir_name,K,n,LFCg,pDEg,beta0[1],phi0[1])   # just the first elt of beta0 and phi0
  }else{file_name=save_file}

  match=match.call()
  # add N(0,sigma) noise to each gene/cluster
  # introduced sigma_g: draw LFC for each DE gene from N(LFCg,sigma_g). same w/ sigma_b for batch
  # introduced sigma: small Gaussian noise N(0,sigma) introduced to beta for all gene/cluster combinations
  #### this would make performance worse throughout
  # SF: defaults to be simulated N(1,0.25), cutoff at bottom 0.25, cutoff at top by 2 (based on BRCA sample estimates)

  # all the diff inputs of beta:
  # 1) scalar --> change to matrix(beta,nrow=g,ncol=K)
  # 2) vector of length g --> change to matrix(beta,nrow=g,ncol=K)
  # 3) matrix of dim 1x1 --> change to matrix(beta,nrow=g,ncol=K)
  # 4) matrix of dim gx1 or 1xg --> change to matrix(beta, nrow=g,ncol=K)
  # 5) matrix of dim gxK --> matrix(beta, nrow=g,ncol=K)
  if(is.null(pK)){     # custom pN (probability of being in each cl K) must be of length K, and add up to 1
    pK=rep(1/K,K)
    pB=rep(1/B,B)
  }
  if(disp=="gene"){
    if(length(phi0)==1){
      phi=rep(phi0,g)
    } else if(length(phi) != g){
      stop("phi must be either length 1 or g")
    }
  } else if(disp=="cluster"){
    if(length(phi0)==K){
      phi=matrix(phi0,nrow=g,ncol=K,byrow=TRUE)
    } else if(length(phi0)!=(g*K)){
      stop("phi must be gxK matrix")
    }
  } else{
    stop("specify 'gene' or 'cluster' for disp")
  }

  beta=matrix(beta0,nrow=g,ncol=K)

  for(i in 1:nsims){

    if(is.null(SF)){
      SF=rnorm(n,1,0.25)
      SF[SF<0.25]=0.25
      SF[SF>2]=2

      SF_pred=rnorm(n_pred,1,0.25)   # n_pred < n
      SF_pred[SF_pred<0.25]=0.25
      SF_pred[SF_pred>2]=2
    } else if(length(SF)!=n){
      stop("Custom size factors must be of length n")
    } else{
      SF_pred=sample(SF,n_pred,replace=T)
    }

    cls=sample(1:K,n,replace=T,prob=pK)
    batch=sample(1:B,n,replace=T,prob=pB)-1

    DEg_ID=rep(F,g)                                 # DEg across cls: just first g*pDEg genes
    DEg_ID[1:floor(pDEg*g)]=T
    DEb_ID=rep(F,g)                                 # select DEb across batch genes randomly
    DEb_ID[sample(1:g,floor(pDEb*g),replace=F)]=T

    LFCg_mat = matrix(0,nrow=g,ncol=K)
    #noise_mat = matrix(rnorm(g,0,sigma_g),nrow=g,ncol=K)


    for(j in 1:floor(pDEg*g)){
      signLFC=(j<floor(pDEg*g/2))*2-1           # +1 for j in 1:(DEgenes/2), -1 for rest
      LFCg_mat[j,sample(1:K,1,replace=T)] = signLFC*LFCg
    }

    cts<-simulate_counts(K=K,B=B,n=n,g=g,
                         cls=cls,SF=SF,
                         beta=beta,phi=phi,LFCg_mat=LFCg_mat,#noise_mat=noise_mat,
                         batch=batch,LFCb=LFCb,sigma_b=sigma_b,sigma_g=sigma_g,DEb_ID=DEb_ID,disp=disp) # gene-specific disp param

    cls_pred=sample(1:K,n_pred,replace=T,prob=pK)
    batch_pred=rep(0,n_pred)     # other batches are centered around 0. prediction batch assumed to be in the middle --> batch=0
    cts_pred <- simulate_counts(K=K,B=1,n=n_pred,g=g,
                                cls=cls_pred,SF=SF_pred,
                                beta=beta,phi=phi,LFCg_mat=LFCg_mat,#noise_mat=noise_mat,
                                batch=batch_pred,LFCb=LFCb,sigma_b=sigma_b,sigma_g=sigma_g,DEb_ID=DEb_ID,disp=disp) # gene-specific disp param

    sim_params=list(K=K,B=B,g=g,n=n,n_pred=n_pred,pK=pK,pB=pB,
                    LFCg=LFCg,pDEg=pDEg,sigma_g=sigma_g,
                    LFCb=LFCb,pDEb=pDEb,sigma_b=sigma_b,
                    beta=beta,phi=phi,disp=disp,
                    LFCg_mat=LFCg_mat,
                    DEb_ID=DEb_ID)
    sim.dat=list(cts=cts, cts_pred=cts_pred,
                 cls=cls, cls_pred=cls_pred, batch=batch,
                 SF=SF, SF_pred=SF_pred,
                 DEg_ID=DEg_ID, sim_params=sim_params)

    # order params in file name: K, n, LFCg, pDEg, beta, phi (use universally)
    file.name=sprintf("%s/%s_sim%d_data.RData",
            dir_name,file_name,i)
    save(sim.dat,file=file.name)
  }

  file.name2=sprintf("%s/%s_dataParams.txt",dir_name,file_name)
  sink(file=file.name2)

  print(match)
  cat(paste("\nB=",B," ,g=",g,", sigma_g=",sigma_g,", sigma_b=",sigma_b,"\n",
            sep=""))
  cat(paste("LFCb=",LFCb,", pDEb=",pDEb,", disp=",disp,", n_pred=",n_pred,"\n",
            sep=""))
  cat("\npK:\n")
  write.table(pK,quote=F,col.names=F)
  cat("\npB:\n")
  write.table(pB,quote=F,col.names=F)
  cat("\nlast sim SF:\n")
  write.table(SF,quote=F,col.names=F)
  cat("\nlast sim cls:\n")
  write.table(cls,quote=F,col.names=F)
  cat("\nlast sim SF_pred:\n")
  write.table(SF_pred,quote=F,col.names=F)
  cat("\nlast sim cls_pred:\n")
  write.table(cls_pred,quote=F,col.names=F)

  sink()
}
