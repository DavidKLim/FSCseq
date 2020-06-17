#' Simulate data with gene-wise dispersion parameters
simulate_counts=function(K,B,n,g,
                         cls,SF,
                         beta,phi,LFCg_mat,
                         batch,batch_effects,sigma_b,sigma_g,DEb_ID,disp){     # batch of 0: no batch effect, batch of 1: yes effect
  # center batch at 0: lower batch gets downregulated, higher batch gets upregulated
  # if B=1, batch_eff = 0 with some noise by sigma_b

  #batch_eff = (batch-(B-1)/2)*LFCb+rnorm(n,0,sigma_b)  # previous parametrization: -0.5*LFCb or +0.5*LFCb for B=2, pred had +0*LFCb
  #batch_effects = c(-1,1) ## FOR B=2. Need general way --> user input.
  # batch_effects should be the same length as the number of batches, or length(unique(batch))
  batch_eff = batch_effects[batch+1] # batch goes from 0, 1, ... for batch_pred, just 0 input for batch --> 1 batch.

  cts <- matrix(rep(0,times=g*n),nrow=g)

  b = beta + LFCg_mat # + noise_mat

  if(disp=="gene"){
    phi_mat = matrix(phi,nrow=g,ncol=K)      # construct matrix with same phi values for all cls if g-disp
  }

  sigma_mat = matrix(rnorm(n*g,0,sigma_g),nrow=g,ncol=n)
  mu_mat = matrix(0,nrow=g,ncol=n)
  for(i in 1:n){
    mu_mat[,i] = SF[i]*2^(b[,cls[i]] + sigma_mat[,i] + DEb_ID*batch_eff[i])  # this was how data was simulated, and is consistent with DESeq2 paper
    #mu_mat[,i] = 2^(SF[i] + b[,cls[i]] + sigma_mat[,i] + DEb_ID*batch_eff[i])     ##### this shouldn't be right, but empirically giving better results?! why?
    cts[,i] = rnbinom( g, size = 1/phi_mat[,cls[i]], mu = mu_mat[,i] )
  }

  return(list(cts=cts,mu_mat=mu_mat,sigma_mat=sigma_mat))
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
#' @param sigma_g numeric, Gaussian noise added to each gene/sample N(0,sigma_g). Default is 0.1
#' @param LFCb numeric, LFC for genes that are differentially expressed across batch. Default is 1.
#' @param pDEb numeric, proportion of genes that are differentially expressed across batch. Default is 0.5.
#' @param sigma_b numeric, batch-specific Gaussian noise (default 0).
#' @param beta0 numeric, baseline log2 expression for each gene before LFC is applied
#' @param phi0 numeric, baseline overdispersion for each gene
#' @param SF vector of length n (optional), custom size factors from DESeq2. If NULL, simulated from N(1,0.25)
#' @param nsims integer, number of datasets to simulate given the input conditions. Default is 25.
#' @param disp string, either 'gene' or 'cluster' to simulate gene-level or cluster-level dispersions. Default is gene-level. Input phi must be g x K matrix if disp='cluster'
#' @param n_pred integer, number of samples in simulated prediction dataset. Default is 25
#' @param sim_batch_pred boolean: FALSE (no batch effect for prediction samples) or TRUE (batch effect)
#' @param LFCb_pred LFCb for batch-affected genes in prediction set. By default, = max(batch_effects) + LFCb/2.
#' @param save_file boolean: TRUE (save each set of simulations)
#' @param save_dir string (optional): directory to save files. Default: 'Simulations/<sigma_g>_<sigma_b>/B<B>'
#' @param save_pref string (optional): prefix of file name to save simulated data to. Default: '<K>_<n>_<LFCg>_<pDEg>_<beta0>_<phi0>'
#'
#' @return if save_file=TRUE, then saved file in '<save_dir>/<save_pref>_sim<1:nsims>_data.RData'. Otherwise, list of length 'nsims', with a sim.dat list object for each simulation
#'
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
#'
#' @examples
#' sim.dat = FSCseq::simulateData(B=1, g=10000, K=2, n=50, LFCg=1, pDEg=0.05, beta0=12, phi0=0.35, nsims=1, save_file=F)[[1]]
#'
#' @export
simulateDataBatch<-function(K=2, B=1, g=10000, n=50, pK=NULL, pB=NULL,
                            LFCg=1, pDEg=0.05, sigma_g=0.1,
                            LFCb=1, pDEb=0.5, sigma_b=0,
                            beta0=12, phi0=0.35, SF=NULL,
                            nsims=25, disp="gene", n_pred=25, sim_batch_pred=FALSE, LFCb_pred=NULL,
                            save_file=TRUE, save_dir=NULL, save_pref=NULL){
  if(B==1){LFCb=0}

  if(save_file){
    if(is.null(save_dir)){
      dir_name1=sprintf("Simulations/%f_%f",sigma_g,sigma_b)
      dir_name=sprintf("%s/B%d_LFCb%d",dir_name1,B,LFCb)
      ifelse(!dir.exists(dir_name1),
             dir.create(dir_name1),
             FALSE)
    }else{
      dir_name=save_dir
    }
    ifelse(!dir.exists(dir_name),
           dir.create(dir_name),
           FALSE)

    if(is.null(save_pref)){
      file_name=sprintf("%d_%d_%f_%f_%f_%f",
                        K,n,LFCg,pDEg,beta0[1],phi0[1])   # just the first elt of beta0 and phi0
    }else{file_name=save_pref}
  }
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
      stop("phi must be length 1 or g")
    } else{
      phi=phi0     # specifying specific gene-level disp's for all genes
    }
  } else if(disp=="cluster"){
    if(length(phi0)==K){
      phi=matrix(phi0,nrow=g,ncol=K,byrow=TRUE)
    } else if(length(phi0)!=(g*K)){
      stop("phi must be gxK matrix")
    } else{
      phi=phi0     # specifying specific cluster-level disps for all genes/clusters
    }
  } else{
    stop("specify 'gene' or 'cluster' for disp")
  }

  if(length(beta0) %in% c(1,g)){
    beta=matrix(beta0,nrow=g,ncol=K)    # if beta0 is length g, then will repeat for each cluster
  } else{
    stop("beta0 must be length 1 (common baseline) or g (gene-level baselines)")
  }
  all_sim.dats = list()
  for(i in 1:nsims){

    # order params in file name: K, n, LFCg, pDEg, beta, phi (use universally)
    if(save_file){
      # skip files that are already saved (already simulated)
      file.name=sprintf("%s/%s_sim%d_data.RData",
                        dir_name,file_name,i)
      if(file.exists(file.name)){warning("Data with same file name already simulated. Remove existing file to replace. Skipping simulation of existing file"); next}else{print(file.name)}
    }

    if(is.null(SF)){
      # SF's are drawn from N(1,0.25), truncated at (0.25, 2.00)
      SF=rnorm(n,1,0.25)
      SF[SF<0.25]=0.25
      SF[SF>2]=2

      SF_pred=rnorm(n_pred,1,0.25)   # n_pred < n
      SF_pred[SF_pred<0.25]=0.25
      SF_pred[SF_pred>2]=2
    } else if(length(SF)!=n){
      stop("Custom size factors must be of length n")
    } else{
      # if custom SF's input, then SF_pred are sampled from the custom SF's randomly
      SF_pred=sample(SF,n_pred,replace=T)
    }

    cls=sample(1:K,n,replace=T,prob=pK)
    batch=sample(1:B,n,replace=T,prob=pB)-1

    DEg_ID=rep(F,g)                                 # DEg across cls: first g*pDEg genes
    DEg_ID[1:floor(pDEg*g)]=T
    DEb_ID=rep(F,g)                                 # select DEb across batch genes randomly
    DEb_ID[sample(1:g,floor(pDEb*g),replace=F)]=T

    LFCg_mat = matrix(0,nrow=g,ncol=K)

    for(j in 1:floor(pDEg*g)){
      signLFC=(j<floor(pDEg*g/2))*2-1           # +1 for j in 1:(DEgenes/2), -1 for rest
      LFCg_mat[j,sample(1:K,1,replace=T)] = signLFC*LFCg
    }

    # make batch_effects be centered around LFCb, +- 0.5 between each consecutive batch

    batch_effects = LFCb*( c(1:B)-mean(c(1:B)) )

    fit<-simulate_counts(K=K,B=B,n=n,g=g,
                         cls=cls,SF=SF,
                         beta=beta,phi=phi,LFCg_mat=LFCg_mat,#noise_mat=noise_mat,
                         batch=batch,batch_effects=batch_effects,sigma_b=sigma_b,sigma_g=sigma_g,DEb_ID=DEb_ID,disp=disp) # gene-specific disp param
    cts = fit$cts; mu_mat = fit$mu_mat; sigma_mat=fit$sigma_mat

    cls_pred=sample(1:K,n_pred,replace=T,prob=pK)

    batch_pred = rep(B, n_pred)     # train batch: 0, 1, .. B-1. pred_batch: B. right now: ONLY ONE BATCH simulated for prediction set
    batch_effects_pred =
      if(sim_batch_pred){
        if(is.null(LFCb_pred)){
          max(batch_effects) + LFCb/2 # if batch effects were (-1,1), then this batch effect is 2: "slightly greater" effect. Is this even possible to correct for??
        } else{LFCb_pred}
      }else{ 0 }
    fit_pred <- simulate_counts(K=K,B=1,n=n_pred,g=g,
                                cls=cls_pred,SF=SF_pred,
                                beta=beta,phi=phi,LFCg_mat=LFCg_mat,#noise_mat=noise_mat,
                                batch=rep(0,n_pred), batch_effects = batch_effects_pred,sigma_b=sigma_b,sigma_g=sigma_g,DEb_ID=DEb_ID,disp=disp) # gene-specific disp param

    cts_pred = fit_pred$cts; mu_mat_pred = fit_pred$mu_mat; sigma_mat_pred = fit_pred$sigma_mat

    sim_params=list(K=K,B=B,g=g,n=n,n_pred=n_pred,pK=pK,pB=pB,
                    LFCg=LFCg,pDEg=pDEg,sigma_g=sigma_g,
                    LFCb=LFCb,pDEb=pDEb,sigma_b=sigma_b,
                    batch_effects=batch_effects,batch_effects_pred=batch_effects_pred,
                    beta=beta,phi=phi,disp=disp,
                    LFCg_mat=LFCg_mat,
                    DEb_ID=DEb_ID, mu_mat=mu_mat, mu_mat_pred=mu_mat_pred, sigma_mat=sigma_mat, sigma_mat_pred=sigma_mat_pred)
    sim.dat=list(cts=cts, cts_pred=cts_pred,
                 cls=cls, cls_pred=cls_pred, batch=batch, batch_pred=batch_pred,
                 SF=SF, SF_pred=SF_pred,
                 DEg_ID=DEg_ID, sim_params=sim_params)

    if(save_file){
      save(sim.dat,file=file.name)
    }else{all_sim.dats[[i]]=sim.dat}

  }
  if(save_file){
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

    sink()
  }else{
    return(all_sim.dats)
  }
}
