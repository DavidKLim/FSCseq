#' Logsumexp technique to prevent underflow
#'
#' This function takes in a vector, and uses the logsumexp
#' technique to prevent underflow issues. Used on
#' log_pi + loglike to prevent it being rounded to 0
#'
#' @param v a vector or matrix of numbers
#'
#' @return log(sum(exp(v))), without underflow issues
#'
#' @examples
#' log_pi=log(c(0.1,0.02,0.07))
#' l=matrix(c(-2000,-2100,-2050,-1000,-1050,-1100),nrow=3,ncol=2)
#' logsumexpc(log_pi+l)
#' log(sum(exp(log_pi+l)))
#'
#' @export
logsumexpc=function(v){
  if(any(is.infinite(v))){
    warning("infinite value in v\n")
  }
  if(length(v)==1){ return(v[1]) }
  sv = sort(v, decreasing=TRUE)
  res = sum(exp(sv[-1] - sv[1]))
  lse = sv[1] + log(1+res)
  lse
}

#' SCAD soft thresholding function
#'
#' Takes pairwise distances between cluster log2 means and
#' penalty parameters lambda and alpha, and imposes
#' SCAD thresholding rule
#'
#' @param diff_beta difference in pairwise cluster betas (log2 means)
#' @param lambda numeric penalty parameter, lambda >= 0
#' @param alpha numeric penalty parameters, 0 <= alpha < 1
#'
#' @return diff_beta value, after SCAD thresholding rule is applied
#'
#' @examples
#' beta1=1
#' beta2=2
#' SCAD_soft_thresholding(beta2-beta1,0,0)
#' SCAD_soft_thresholding(beta2-beta1,0.5,0.45)
#' SCAD_soft_thresholding(beta2-beta1,1,1)
#'
#' @export
SCAD_soft_thresholding=function(diff_beta,lambda,alpha){
  a=3.7
  #if(abs(diff_beta)<=2*lambda*alpha){
  if(abs(diff_beta)<=(alpha/(1-alpha))+lambda*alpha ){
    if(abs(diff_beta)<=alpha/(1-alpha)){
      return(0)
    } else{
      return(sign(diff_beta)*(abs(diff_beta)-alpha/(1-alpha)))
    }
  }else if(abs(diff_beta)>(alpha/(1-alpha))+lambda*alpha & abs(diff_beta)<=a*lambda*alpha){
    omega = ((a-1)*diff_beta)/(a-1-1/(lambda*(1-alpha)))
    if(abs(omega)-(a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha))) <= 0){
      return(0)
    } else{
      return(sign(omega)*(abs(omega)-(a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha)))) )
    }
  }else{
    return(diff_beta)
  }
}

#' Computes E step weights
#'
#' Takes current parameter estimates and likelihoods, and
#' computes E step weights. Also updates CEM temperature Tau,
#' if applicable
#'
#' @param wts matrix of dimension k by n, current E step weights. Used as template matrix to update
#' @param l matrix of dimension k by n, log-likelihood summed over all genes, for each k clusters and n subjects
#' @param pi vector of k mixture proportions
#' @param CEM logical, FALSE for EM and TRUE for CEM
#' @param Tau numeric, current temperature of CEM. Tau=1 for EM
#'
#' @return list containing updated wts: matrix of E step weights,
#' keep: logical matrix of whether a sample is included in cluster baseline calculation,
#' Tau: numeric: resultant CEM temperature (after annealing),
#' CEM: logical: whether CEM is continued to be run
#'
#' @export
E_step<-function(wts,l,pi,CEM,Tau){
  k=length(pi)
  n=ncol(l)

  if(!CEM){
    # update on weights
    logdenom = apply(log(pi) + l, 2,logsumexpc)
    for(c in 1:k){
      wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
    }
  } else if(CEM){
    # CEM update on weights
    logdenom = apply((1/Tau)*(log(pi)+l),2,logsumexpc)
    for(c in 1:k){
      wts[c,]<-exp((1/Tau)*(log(pi[c])+l[c,])-logdenom)
    }
    if(Tau>1){
      Tau = 0.9*Tau
    } else{
      Tau=1       # after Tau hits 1 --> fix at 1
      CEM=F       # after Tau is fixed at 1 --> no more stochastic sampling
    }
  }

  # UB and LB on weights
  for(i in 1:n){
    for(c in 1:k){
      if(is.na(wts[c,i])){
        wts[c,i]=1E-50
      } else if(wts[c,i]<1E-50){
        wts[c,i]=1E-50
      } else if(wts[c,i]>(1-1E-50)){
        wts[c,i]=1-1E-50
      }
    }
  }

  if(CEM){
    # CEM
    draw_wts=wts                 # initialize
    for(i in 1:n){
      set.seed(i)
      draw_wts[,i] = rmultinom(1,1,wts[,i])
    }
    seed_mult=1
    while(any(rowSums(draw_wts)==0)){
      #if(trace){cat("Drawing again",seed_mult,"\n")}
      for(i in 1:n){
        set.seed(seed_mult*n+i)
        for(c in 1:k){
          if(wts[c,i]<=(1E-50*10^seed_mult) & seed_mult<=48){
            wts[c,i]=1E-50*10^seed_mult
          } else if(wts[c,i]>=(1-(1E-50*10^seed_mult)) & seed_mult<=48){
            wts[c,i]=1-1E-50*10^seed_mult
          }
        }
        draw_wts[,i] = rmultinom(1,1,wts[,i])
      }
      seed_mult=seed_mult+1
      if(seed_mult>250){
        draw_wts[,n]=rep(1/k,k)
        break
      }
    }
    wts=draw_wts
  }          # Keep drawing until at least one in each cluster

  # Input in M step only samples with PP's > 0.001
  keep = (wts>0.001)^2      # matrix of 0's and 1's, dimensions k x n

  return(list(wts=wts,keep=keep,Tau=Tau,CEM=CEM))
}

#' Initializes parameter estimates
#'
#' Uses glm() and phi_ml_g() functions to initialize
#' estimates at the beginning of EM algorithm.
#'
#' @param j integer, gene index j
#' @param y_j vector of length n, expression of gene j
#' @param XX matrix of dimension nk by k+p, design matrix
#' @param k integer, number of clusters
#' @param offsets numeric vector of length n, subject-specific offsets
#' @param wts matrix of dimension k by n, E step weights
#' @param keep logical matrix of dimension k by n, whether a sample is included in calculation for cluster k
#'
#' @return list containing outputs coefs_j: vector of length k corresponding to each cluster log2 baseline,
#' phi_g: Overdispersion estimate
#'
#' @importFrom stats glm
#' @importFrom MASS theta.ml
#'
#' @export
glm.init=function(j,y_j,XX,k,offsets,wts,keep){
  ids = c(t(keep==1))  # PP filtering

  init.fit = stats::glm(as.integer(rep(y_j,k))[ids]~0+XX[ids,]+offset(rep(offsets,k)[ids]),family=poisson(link="log"),weights=c(t(wts))[ids])
  coefs_j = log2(exp(init.fit$coefficients))         # change to log2 scale

  coefs_j[is.na(coefs_j)] = 0
  mu = 2^(XX %*% coefs_j + offsets)

  # phi_ml_g() is actually a bit slower than theta.ml(), but more stable
  phi_g = phi_ml_g(y=as.integer(rep(y_j,k))[ids],
                   mu=mu[ids],
                   wts=c(t(wts))[ids],
                   limit=25,
                   p0=0,
                   trace=0)

  # phi_g = 1/MASS::theta.ml(y=as.integer(rep(y[j,],k))[ids],mu=mu[ids],weights=c(t(wts))[ids],limit=25,trace=F)


  results=list(coefs_j=coefs_j,
               phi_g=phi_g)
  return(results)
}


#' Parallelized glm.init
#'
#' Uses glm() and phi_ml_g() functions to initialize
#' estimates at the beginning of EM algorithm.
#'
#' @param j integer, gene index j
#'
#' @return list containing outputs coefs_j: vector of length k corresponding to each cluster log2 baseline,
#' phi_g: Overdispersion estimate
#'
#' @importFrom stats glm
#' @importFrom MASS theta.ml
#'
#' @export
glm.init_par=function(j){
  ids = c(t(keep==1))  # PP filtering

  init.fit = stats::glm(as.integer(rep(y[j,],k))[ids]~0+XX[ids,]+offset(rep(offsets,k)[ids]),family=poisson(link="log"),weights=c(t(wts))[ids])
  coefs_j = log2(exp(init.fit$coefficients))         # change to log2 scale

  coefs_j[is.na(coefs_j)] = 0
  mu = 2^(XX %*% coefs_j + offsets)

  # phi_ml_g() is actually a bit slower than theta.ml(), but more stable
  phi_g = phi_ml_g(y=as.integer(rep(y[j,],k))[ids],
                   mu=mu[ids],
                   wts=c(t(wts))[ids],
                   limit=25,
                   p0=0,
                   trace=0)

  # phi_g = 1/MASS::theta.ml(y=as.integer(rep(y[j,],k))[ids],mu=mu[ids],weights=c(t(wts))[ids],limit=25)

  results=list(coefs_j=coefs_j,
               phi_g=phi_g)
  return(results)
}

#' Parallelized M step
#'
#' Runs M_step() for gene j. Used to compute in parallel on ncores (specified in wrapper)
#'
#' @param j gene index
#'
#' @return Same output as M_step() function
#'
#' @export
M_step_par = function(j){
  if(Tau<=1 & a>6){
    if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){
      res=par_X[[j]]
      return(res)
    }}
  res = M_step(X=XX, y_j=as.numeric(rep(y[j,],k)), p=p, j=j, a=a, k=k,
               all_wts=wts, keep=c(t(keep)), offset=rep(offsets,k),
               theta=theta_list[[j]],coefs_j=coefs[j,],phi_j=phi[j,],
               cl_phi=cl_phi,est_phi=est_phi[j],est_covar=est_covar[j],
               lambda=lambda,alpha=alpha,IRLS_tol=IRLS_tol,maxit_IRLS=maxit_IRLS,
               optim_method=optim_method)
  return(res)
}

#' Parallelized M step (mclapply)
#'
#' Runs M_step() for gene j. Used to compute in parallel on ncores (specified in wrapper)
#'
#' @param j gene index
#'
#' @return Same output as M_step() function
#'
#' @export
M_step_par2 = function(j, XX, y, p, a, k,
                      wts, keep, offsets,
                      theta_list, coefs, phi,
                      cl_phi, est_phi, est_covar,
                      lambda, alpha, IRLS_tol, maxit_IRLS,
                      optim_method, Tau, disc_ids_list, par_X){
  if(Tau<=1 & a>6){
    if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){
      res=par_X[[j]]
      return(res)
    }}
  res = M_step(X=XX, y_j=as.numeric(rep(y[j,],k)), p=p, j=j, a=a, k=k,
               all_wts=wts, keep=c(t(keep)), offset=rep(offsets,k),
               theta=theta_list[[j]],coefs_j=coefs[j,],phi_j=phi[j,],
               cl_phi=cl_phi,est_phi=est_phi[j],est_covar=est_covar[j],
               lambda=lambda,alpha=alpha,IRLS_tol=IRLS_tol,maxit_IRLS=maxit_IRLS,
               optim_method=optim_method)
  return(res)
}


#' Wrapper for main FSCseq function
#'
#' Determine and search HC, KM, and random initializations by BIC,
#' and then full EM/CEM run based on the optimal one.
#' Calls EM_run function to perform the clustering.
#'
#' @param X (optional) design matrix of dimension n by p
#' @param y count matrix of dimension g by n
#' @param k integer, number of clusters
#' @param lambda numeric penalty parameter, lambda >= 0
#' @param alpha numeric penalty parameters, 0 <= alpha < 1
#' @param size_factors numeric vector of length n, factors to correct for subject-specific variation of sequencing depth
#' @param norm_y count matrix of dimension g by n, normalized for differences in sequencing depth
#' @param true_clusters (optional) integer vector of true groups, if available, for diagnostic tracking
#' @param true_disc (optional) logical vector of true discriminatory genes, if available, for diagnostic tracking
#' @param init_parms logical, TRUE: custom parameter initializations, FALSE (default): start from scratch
#' @param init_coefs matrix of dimension g by k, only if init_parms = TRUE
#' @param init_phi vector of dimension g (gene-specific dispersions) or matrix of dimension g by k (cluster-specific dispersions), only if init_parms = TRUE
#' @param init_cls (optional) vector of length n, initial clustering. If NA (default), multiple initializations will be searched
#' @param init_wts (optional) matrix of dim k by n to denote initial clustering (allowing partial membership). If both init_cls and init_wts specified, init_wts will be ignored and init_cls used as initial clusters
#' @param init_method if searching over (n_rinits+2) initializations (random + HC and KM inits) --> how to choose optimal init to run full EM.
#' @param n_rinits integer, number of additional random initializations to be searched (default 50 for EM, 10 for CEM)
#' @param maxit_inits integer, maximum number of iterations for each initialization search (default 15 for EM, or until temperature anneals down to below 2 for CEM)
#' @param maxit_EM integer, maximum number of iterations for full EM/CEM run (default 100)
#' @param maxit_IRLS integer, maximum number of iterations for IRLS algorithm, in M step (default 50)
#' @param EM_tol numeric, tolerance of convergence for EM/CEM, default is 1E-8
#' @param IRLS_tol numeric, tolerance of convergence for IRLS, default is 1E-6
#' @param disp string, either "gene" (default) or "cluster"
#' @param optim_method string, three options "direct", "NR", or "GD". Direct, Newton-Raphson, or Gradient descent (fixed step size of 2)
#' @param method string, either "EM" (default) or "CEM"
#' @param init_temp numeric, default for CEM: init_temp = nrow(y), i.e. number of genes. temp=1 for EM
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
#' @param mb_size minibatch size: # of genes to include per M step iteration
#'
#' @return list containing outputs from EM_run() function
#'
#' @export
FSCseq<-function(ncores=1,X=NULL, y, k,
                 lambda=0,alpha=0,
                 size_factors=rep(1,times=ncol(y)),
                 norm_y=y,
                 true_clusters=NULL, true_disc=NULL,
                 init_parms=FALSE,
                 init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                 init_phi=matrix(0,nrow=nrow(y),ncol=k),
                 init_cls=NULL,init_wts=NULL,
                 init_method="max",n_rinits=if(method=="EM"){50}else{10},         # fewer searches for CEM to minimize computational cost
                 maxit_inits=if(method=="EM"){15}else{ceiling(log(2/nrow(y))/log(0.9))}, # for CEM, tolerates end temp (Tau) of 2 at end of initialization
                 maxit_EM=100,maxit_IRLS=50,EM_tol=1E-6,IRLS_tol=1E-6,
                 disp=c("gene","cluster"),optim_method="direct",
                 method=c("EM","CEM"),init_temp=sqrt(nrow(y)),trace=F,trace.file=NULL,
                 mb_size=NULL){

  # y: raw counts
  # k: #clusters
  # size_factors: SF's derived from DESeq2
  # norm_y: counts normalized for sequencing depth by DESeq2
  # true_clusters: if applicable. For diagnostics tracking of ARI
  # true_disc: if applicable. For diagnostics tracking of disc/nondisc genes
  # init_parms: TRUE if initial coefficient estimates/dispersion estimates are input
  # init_coefs & init_phi: Initial estimates, if applicable
  # disp = c(gene, cluster), depending on whether dispersions are gene-specific or cluster-specific
  # init_cls: Initial clustering
  # n_rinits: Number of initial clusterings searched with maxit=15. More initializations = more chance to attain global max

  n = ncol(y)
  g = nrow(y)
  cat(paste(method,"model with",disp,"level dispersions specified.\n"))

  if(is.null(X)){
    cat("No covariates specified. Running cluster-specific intercept-only model.\n")
  } else{
    if (class(X) != "matrix") {
      tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
      if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
    }
    if (storage.mode(X)=="integer") storage.mode(X) <- "double"
    if (ncol(y) != nrow(X)) stop("X and y do not have the same number of observations")
    if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y")
  }

  if(alpha < 0 | alpha >= 1){
    stop("alpha must be 0 <= alpha < 1")
  }
  if(lambda<0){
    stop("lambda must be greater than 0")
  }

  if(method=="EM"){
    CEM=F
  } else if(method=="CEM"){
    CEM=T
  } else{
    stop("method must be 'EM' or 'CEM'.")
  }

  if(!is.null(true_clusters)){
    if(length(true_clusters) != n){
      warning("Length of true clusters not equal to n, can't track diagnostics")
      true_clusters=NULL
    }
  }
  if(!is.null(true_disc)){
    if(length(true_disc) != g){
      warning("Length of true disc gene labels is not equal to g, can't track diagnostics")
      true_disc=NULL
    }
  }
  if(!is.null(init_cls)){
    if(length(init_cls) != n){stop("Length of initial clusters not equal to n")}
  }

  ## Diagnostic file
  if(trace){
    if(!is.null(trace.file)){
      sink(file=trace.file)
    }
  }

  if(trace){
    cat(paste(sprintf("n=%d, g=%d, k=%d, l=%f, alph=%f, ",n,g,k,lambda,alpha),"\n"))
    cat("True clusters:\n")
    write.table(head(true_clusters,n=3),quote=F,col.names=F)
    write.table(tail(true_clusters,n=3),quote=F,col.names=F)
  }

  init_Tau=1               # Tau=1 is classic EM. If CEM, then if k=1 or if there are initial clusters input
                           # then init temp = 1 (assuming initial clusters are good --> don't want to perturb it)

  if(k==1){
    init_cls=rep(1,n)
  }
  if(is.null(init_cls) & is.null(init_wts) & k>1){
    if(CEM){
      init_Tau=init_temp     # temperature Tau = g is default for CEM: only when there are no init cls and when K>1
    }

    # Initial Clusterings
    ## Hierarchical Clustering
    d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
    model<-hclust(d,method="average")       # hierarchical clustering
    cls_hc <- cutree(model,k=k)

    ## K-means Clustering
    cls_km <- kmeans(t(log(norm_y+0.1)),k)$cluster

    #TESTING RANDOM CLUSTERING
    r_it=n_rinits
    rand_inits = matrix(0,nrow=n,ncol=r_it)

    cls_cov_collinear = function(cls,X){
      if(is.null(X)){
        return(FALSE)
      }
      p=ncol(X)
      collinear = rep(NA,p)
      for(l in 1:p){
        tab = table(cls,X[,l])
        rowZeroes = rowSums(tab==0)
        if(all(rowZeroes==ncol(tab)-1)){          # If all categorical variables are same for each respective cluster
          collinear[l]=TRUE
        } else{collinear[l]=FALSE}
      }

      if(sum(collinear)==0){
        return(FALSE)
      } else{return(TRUE)}
    }

    for(r in 1:r_it){
      set.seed(r)
      rand_inits[,r] = sample(1:k,n,replace=TRUE)
      while(sum(1:k %in% rand_inits[,r]) < k | cls_cov_collinear(rand_inits[,r],X)){        # If no sample in one cluster OR covariate and clustering
        rand_inits[,r] = sample(1:k,n,replace=TRUE)                                         # are completely collinear, then resample
      }
    }
    colnames(rand_inits) = paste("rand",c(1:r_it),sep="")

    # Iterate through EM with each initialization
    all_init_cls <- cbind(cls_hc,cls_km,rand_inits)
    n_inits = ncol(all_init_cls)
    init_cls_BIC <- rep(0,times=n_inits)

    all_fits = list()

    for(i in 1:n_inits){
      if(trace){
        cat(paste("INITIAL CLUSTERING:",colnames(all_init_cls)[i],"\n"))
      }
      if(!is.null(X)){
        # Check for 100% collinearity between predictor and initial clusters --> will not work. skip this init
        if(cls_cov_collinear(all_init_cls[,i],X)){
          init_cls_BIC[i] = NA
          cat(paste("Initial clusters ",colnames(all_init_cls)[i]," is perfectly collinear with X.\n"))
          next
        }
      }

      fit = EM_run(ncores,X,y,k,lambda,alpha,size_factors,norm_y,true_clusters,true_disc,
                   init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi,disp=disp,
                   init_cls=all_init_cls[,i], CEM=CEM,init_Tau=init_Tau,maxit_EM=maxit_inits,
                   maxit_IRLS=maxit_IRLS,EM_tol=EM_tol,IRLS_tol=IRLS_tol,trace=trace,optim_method=optim_method,
                   mb_size=mb_size)
      all_fits [[i]] = fit
      init_cls_BIC[i] <- fit$BIC
    }

    ## Selecting maximum BIC model ##
    if(init_method=="max"){
      fit_id = which.min(init_cls_BIC)
      init_cls = all_fits[[fit_id]]$clusters
      init_parms = T
      init_coefs = all_fits[[fit_id]]$coefs
      init_phi = all_fits[[fit_id]]$phi
      if(trace){
        cat("FINAL INITIALIZATION:\n")
        cat(paste(colnames(all_init_cls)[fit_id],"\n"))
      }
    }

    if(init_method %in% c("emaEM","BIA")){
      #### Consensus methods ####
      # Create Z matrices and calculate weights
      list_Z = list()
      w_BIC_emaEM = rep(NA,n_inits)
      w_BIC_emaEM_logdenom = logsumexpc(-0.5*init_cls_BIC)

      max_BIC = max(init_cls_BIC,na.rm=T)
      w_BIC_BIA = rep(NA,n_inits)
      w_BIC_BIA_logdenom = logsumexpc(-0.5*(init_cls_BIC-max_BIC))
      for(m in 1:n_inits){
        list_Z[[m]] = matrix(0,nrow=k,ncol=n)
        for(c in 1:k){
          list_Z[[m]][c,]=(all_init_cls[,m]==c)^2
        }
        # w_BIC_emaEM[m]=if(is.na(init_cls_BIC[m])){0} else{exp(-0.5*init_cls_BIC[m]-w_BIC_emaEM_logdenom)}
        # w_BIC_BIA[m]=if(is.na(init_cls_BIC[m])){0} else{exp(-0.5*(init_cls_BIC[m]-max_BIC)-w_BIC_BIA_logdenom)}
        ##### EXPERIMENT #####
        temp=1     # temp=1 is the same as regular emaEM and BIA
        w_BIC_emaEM[m]=if(is.na(init_cls_BIC[m])){0} else{exp((1/temp)*(-0.5)*init_cls_BIC[m]-logsumexpc((-0.5*init_cls_BIC)*(1/temp)))}
        w_BIC_BIA[m]=if(is.na(init_cls_BIC[m])){0} else{exp((1/temp)*(-0.5)*(init_cls_BIC[m]-max_BIC)-logsumexpc((-0.5*(init_cls_BIC-max_BIC))*(1/temp)))}
      }

      ## emaEM method (Michael & Melnykov 2016) ##
      if(init_method == "emaEM"){
        # create y_m, A_m for all inits --> A --> HC --> replace init_cls
        A = matrix(0,n,n)
        for(m in 1:n_inits){
          # construct A_m matrix for all inits (y_m omitted to conserve memory. can obtain if you take out w_BIC_emaEM[m] = 1)
          # sum directly here to save memory
          A_m=matrix(NA,n,n)
          for(i in 1:n){
            for(ii in 1:n){
              A_m[i,ii] = w_BIC_emaEM[m]*all(list_Z[[m]][,i] == list_Z[[m]][,ii])^2   # if all wts for sample i and sample ii are the same --> y_m[i,ii]=1. otherwise, 0
            }
          }
          A = A + A_m
        }
        model = hclust(as.dist(A),method="average")
        init_cls <- cutree(hclust(as.dist(A),method="average"),k=k)
      }

      ## BIA method (Hagan & White 2018) ##
      if(init_method == "BIA"){
        # create Z* = sum(w*Z) --> replace init_wts
        init_wts=matrix(0,nrow=k,ncol=n)
        for(m in 1:n_inits){
          init_wts = init_wts + w_BIC_BIA[m]*list_Z[[m]]
        }
      }
    }

  }


  # if init_method="max": will have specified init_cls, init_parms=T, init_coefs, and init_phi
  # if init_method="emaEM": will have specified init_cls, init_parms=F. init_wts will remain NULL
  # if init_method="BIA": will have specified init_wts, init_parms=F. init_cls will remain NULL

  # CEM is set to F for final run to convergence. If inits were searched by
  # CEM, then temperature would have already reached 1 --> no need for CEM.
  # no inits searched for K=1 --> CEM and EM are equivalent for K=1.
  results=EM_run(ncores,X,y,k,lambda,alpha,size_factors,norm_y,true_clusters,true_disc,
                 init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi,disp=disp,
                 init_cls=init_cls,init_wts=init_wts,CEM=F,init_Tau=1,
                 maxit_EM=maxit_EM,maxit_IRLS=maxit_IRLS,EM_tol=EM_tol,IRLS_tol=IRLS_tol,trace=trace,optim_method=optim_method,
                 mb_size=mb_size)

  if(trace){
    if(!is.null(trace.file)){
      sink()
    }
  }
  return(results)
}

#' EM/CEM run for FSCseq
#'
#' Performs clustering, feature selection, and estimation of parameters
#' using a finite mixture model of negative binomials
#'
#' @param X design matrix of dimension n by p
#' @param y count matrix of dimension g by n
#' @param k integer, number of clusters
#' @param lambda numeric penalty parameter, lambda >= 0
#' @param alpha numeric penalty parameters, 0 <= alpha < 1
#' @param size_factors numeric vector of length n, factors to correct for subject-specific variation of sequencing depth
#' @param norm_y count matrix of dimension g by n, normalized for differences in sequencing depth
#' @param true_clusters (optional) integer vector of true groups, if available, for diagnostic tracking
#' @param true_disc (optional) logical vector of true discriminatory genes, if available, for diagnostic tracking
#' @param init_parms logical, TRUE: custom parameter initializations, FALSE (default): start from scratch
#' @param init_coefs matrix of dimension g by k, only if init_parms = TRUE
#' @param init_phi vector of dimension g (gene-specific dispersions) or matrix of dimension g by k (cluster-specific dispersions), only if init_parms = TRUE
#' @param init_cls vector of length n, initial clustering. If NA (default), multiple initializations will be searched. init_wts or init_cls must be initialized
#' @param init_wts matrix of dim k x n: denotes cluster memberships, but can have partial membership. init_wts or init_cls must be initialized
#' @param CEM logical, TRUE for EM (default), FALSE for CEM
#' @param init_Tau numeric, initial temperature for CEM. Default is 1 for EM and g for CEM
#' @param maxit_EM integer, maximum number of iterations for full EM/CEM run (default 100)
#' @param maxit_IRLS integer, maximum number of iterations for IRLS algorithm, in M step (default 50)
#' @param EM_tol numeric, tolerance of convergence for EM/CEM, default is 1E-8
#' @param IRLS_tol numeric, tolerance of convergence for IRLS, default is 1E-6
#' @param disp string, either "gene" (default) or "cluster"
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
#' @param optim_method string, three options "direct", "NR", or "GD". Direct, Newton-Raphson, or Gradient descent (fixed step size of 2)
#' @param mb_size integer, size of minibatch
#'
#' @return FSCseq object: list containing outputs
#' k: integer order,
#' pi: numeric vector of mixture proportions,
#' coefs: numeric matrix of dimensions g by (k+p) cluster log2 baselines and covariate effects,
#' Q: numeric vector of Q function values across iterations,
#' BIC: numeric BIC value,
#' discriminatory: logical vector of length g of whether a gene is discriminatory,
#' init_clusters: numeric vector of clusters at beginning of EM/CEM algorithm,
#' init_coefs: matrix of dimension g by k only output if initial parameter estimates were specified,
#' init_phi: vector of length g or matrix of dimension g by k only output if initial parameter estimates were specified,
#' final_clusters: vector of length n of resulting clusters,
#' phi: dispersion estimates,
#' logL: total log-likelihood,
#' wts: k by n matrix of E step weights,
#' time_elap: amount of time elapsed (in seconds),
#' lambda: input lambda,
#' alpha: input alpha,
#' size_factors: input size factors,
#' norm_y: input norm_y,
#' DNC: 1 if EM didn't converge in maxits and 0 if it did,
#' LFCs: matrix of dimension g by maxit_EM of LFCs (max-min)/(k-1)
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom parallel makeCluster clusterExport stopCluster clusterEvalQ parLapply mclapply
#'
#' @export
EM_run <- function(ncores,X=NA, y, k,
                   lambda=0,alpha=0,
                   size_factors=rep(1,times=ncol(y)) ,
                   norm_y=y,
                   true_clusters=NA, true_disc=NA,
                   init_parms=FALSE,
                   init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                   init_phi=matrix(0,nrow=nrow(y),ncol=k),
                   init_cls=NULL,init_wts=NULL,
                   CEM=F,init_Tau=1,
                   maxit_EM=100, maxit_IRLS = 50,EM_tol = 1E-6,IRLS_tol = 1E-6,disp,trace=F,optim_method="direct",
                   mb_size=NULL){

  start_time <- Sys.time()

  n<-ncol(y)         # number of samples
  g<-nrow(y)         # number of genes
  p<-ncol(X)         # number of covariates
  if(is.null(p)){p=0}

  cl_X = matrix(0,nrow=k*n,ncol=k)
  ident_k = diag(k)
  for(i in 1:k){
    cl_X[((i-1)*n+1):(i*n),] = matrix(rep(ident_k[i,],n),ncol=k,nrow=n,byrow=T)
  }

  if(is.null(X)){
    XX=cl_X
    covars=F
  } else{
    XX = do.call("rbind", replicate(k, X,simplify=FALSE))
    XX = cbind(cl_X,XX)
    covars=T
  }

  if(disp=="gene"){
    cl_phi=0
  } else if(disp=="cluster"){
    cl_phi=1
  }

  # Initialize parameters
  finalwts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  coefs<-matrix(rep(0,times=(g*(k+p))),nrow=g)
  pi<-rep(0,times=k)
  phi= matrix(0,nrow=g,ncol=k)    # initial gene-specific dispersion parameters for negative binomial
  # --> Poisson (phi = 0 = 1/theta)
  theta_list <- list()            # temporary to hold all K x K theta matrices across EM iterations
  temp_list <- list()             # store temp to see progression of IRLS
  phi_list <- list()              # store each iteration of phi to see change with each iteration of EM
  coefs_list <- list()
  MSE_coefs = rep(0,maxit_EM)
  LFCs = matrix(0,nrow=g,ncol=maxit_EM)

  offsets=log2(size_factors)

  Q<-rep(0,times=maxit_EM)

  if(!is.null(init_cls)){
    cls=init_cls
    current_clusters<-cls

    # Initialize weights
    wts<-matrix(0,nrow=k,ncol=n)
    for(c in 1:k){
      wts[c,]=(cls==c)^2
    }
  } else if(!is.null(init_wts)){
    wts=init_wts

    cls=apply(wts,2,which.max)
    current_clusters=cls
  }

  # For use in CEM in E step #
  Tau = init_Tau

  phi_g = rep(0,times=g)
  DNC=0

  disc_ids_list = list()
  disc_ids=rep(T,g)

  diff_phi=matrix(0,nrow=maxit_EM,ncol=g)

  est_phi=rep(1,g)                          # 1 for true, 0 for false
  est_covar = if(covars){rep(1,g)
      } else{
        rep(0,g)
      }

  keep = (wts>0.001)^2

  ## Manual turn-off of estimating phi/covariates
  # if(init_parms){
  #   est_phi=rep(0,g)
  #   est_covar = rep(0,g)
  # }

  # all_temp_list = list()
  # all_theta_list = list()

  lower_K=FALSE         # tracks whether a lower K is necessary

  cl_agreement = rep(NA,maxit_EM)
  par_X=rep(list(list()),g)         # store M_step results

  n_mb = if(!is.null(mb_size)){ceiling(g/mb_size)     # number of minibatches in g. default: 1. experiment with 5 (mb_size = g/5)
    }else{1}


  ########### M / E STEPS #########
  for(a in 1:maxit_EM){
    EMstart= as.numeric(Sys.time())

    prev_clusters = apply(wts,2,which.max)     # set previous clusters (or initial if a=1)
    if(a==1){         # Initializations for 1st EM iteration
      start=as.numeric(Sys.time())
      if(init_parms){
        coefs=init_coefs
        if(cl_phi==0){
          phi_g=init_phi
          phi = matrix(rep(phi_g,k),ncol=k)
        } else{
          phi=init_phi
        }
      } else{
        # Initialization
        if(ncores>1){
          # parallelized Poisson glm + phi_ml_g initialization
          if(Sys.info()[['sysname']]=="Windows"){
            ## use mclapply for Windows ##
            clust = makeCluster(ncores)
            clusterEvalQ(cl=clust,library(FSCseq))
            clusterExport(cl=clust,varlist=c("keep","y","XX","k","offsets","wts"),envir=environment())
            par_init_fit = parallel::parLapply(clust, 1:g, glm.init_par)
            stopCluster(clust)
          } else{
            ## use mclapply for others (Linux/Debian/Mac) ##
            par_init_fit = parallel::mclapply(1:g, mc.cores=ncores, FUN= function(j){
              glm.init(j,y[j,],XX,k,offsets,wts,keep)
            })
          }

          all_init_params=t(sapply(par_init_fit,function(x) {c(x$coefs_j,x$phi_g)}))  # k+p+1 cols
          coefs=matrix(all_init_params[,1:(k+p)],ncol=(k+p))
          phi_g=all_init_params[,ncol(all_init_params)]
          phi=matrix(rep(phi_g,k),ncol=k)
        } else{
          # regular for loop Poisson glm + phi_ml_g init
          for(j in 1:g){
            #j,y_j,XX,k,covars,p,offsets,wts,keep
            init_fit=glm.init(j,y[j,],XX,k,offsets,wts,keep)
            coefs[j,]=init_fit$coefs_j
            phi[j,]=rep(init_fit$phi_g,k)
            phi_g[j]=init_fit$phi_g
          }
        }
      }
      for(j in 1:g){
        theta<-matrix(rep(0,times=k^2),nrow=k)
        for(c in 1:k){
          pi[c]=mean(wts[c,])
          for(cc in 1:k){
            theta[c,cc]<-SCAD_soft_thresholding(coefs[j,c]-coefs[j,cc],lambda,alpha)
          }
        }
        theta_list[[j]] <- theta
      }

      end=as.numeric(Sys.time())
      if(trace){
        cat(paste("Parameter Estimates Initialization Time Elapsed:",end-start,"seconds.\n"))
        cat(paste("Initial coefs (top/bottom 3):\n"))
        write.table(head(coefs,n=3),quote=F)
        write.table(tail(coefs,n=3),quote=F)
        cat(paste("Initial phi (top/bottom 3):\n"))
        write.table(head(phi,n=3),quote=F)
        write.table(tail(phi,n=3),quote=F)
        if(!is.null(true_clusters) & !is.null(true_clusters)){
          cat(paste("Initial ARI:",adjustedRandIndex(prev_clusters,true_clusters),"\n"))
        }
      }
    }

    # M step
    Mstart=as.numeric(Sys.time())

    # minibatching starts at iteration 5 (let EM stabilize first)
    if(a<5){
      mb_genes = 1:g
    } else{
      if(!is.null(mb_size)){
        mb_genes = sample(1:g,mb_size,replace=F)
      }
    }

    if(ncores>1){
      # M_step parallelized across ncores
      if(Sys.info()[['sysname']]=="Windows"){
        ## use parLapply for Windows ##
        clust = makeCluster(ncores)
        clusterEvalQ(cl=clust,library(FSCseq))
        clusterExport(cl=clust,varlist=c("XX","y","p","a","k","wts","keep","offsets",
                                         "theta_list","coefs","phi","cl_phi","est_phi","est_covar",
                                         "lambda","alpha","IRLS_tol","maxit_IRLS","optim_method",
                                         "disc_ids_list","Tau","par_X"),envir=environment())
        par_X_mb = parallel::parLapply(clust, mb_genes, M_step_par)
        stopCluster(clust)
      } else{
        ## use mclapply for others (Linux/Debian/Mac) ##
        par_X_mb = parallel::mclapply(mb_genes, mc.cores=ncores, FUN= function(j){
          M_step_par2(j, XX, y, p, a, k,
                      wts, keep, offsets,
                      theta_list, coefs, phi,
                      cl_phi, est_phi, est_covar,
                      lambda, alpha, IRLS_tol, maxit_IRLS,
                      optim_method, Tau, disc_ids_list, par_X)
        })
      }
      par_X[mb_genes] = par_X_mb     # replace minibatch results
    } else{
      # regular M_step for loop across genes
      for(j in mb_genes){
        if(Tau<=1 & a>6){if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){next}}
        par_X[[j]] <- M_step(X=XX, y_j=as.numeric(rep(y[j,],k)), p=p, j=j, a=a, k=k,
                                     all_wts=wts, keep=c(t(keep)), offset=rep(offsets,k),
                                     theta=theta_list[[j]],coefs_j=coefs[j,],phi_j=phi[j,],
                                     cl_phi=cl_phi,est_phi=est_phi[j],est_covar=est_covar[j],
                                     lambda=lambda,alpha=alpha,IRLS_tol=IRLS_tol,maxit_IRLS=maxit_IRLS,
                                     optim_method=optim_method     # added optim_method: need to update M_step.cpp in package
        )
      }
    }
    Mend=as.numeric(Sys.time())
    if(trace){cat(paste("M Step Time Elapsed:",Mend-Mstart,"seconds.\n"))}


    for(j in mb_genes){
      if(Tau<=1 & a>6){if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){next}}
      coefs[j,] <- par_X[[j]]$coefs_j
      theta_list[[j]] <- par_X[[j]]$theta_j
      disc_ids[j]=any(theta_list[[j]]!=0)
      temp_list[[j]] <- if(p>0){cbind(par_X[[j]]$temp_beta, par_X[[j]]$temp_gamma)}else{par_X[[j]]$temp_beta}
      if(cl_phi==1){
        phi[j,] <- par_X[[j]]$phi_j
      } else if(cl_phi==0){
        phi_g[j] <- (par_X[[j]]$phi_j)[1]
        phi[j,] = rep(phi_g[j],k)
      }

      ### IF any coefs/phi unstable --> missing after M step, re-initialize via glm.nb() for just that gene
      # if(any(is.na(coefs[j,])) | any(is.na(phi[j,]))){
      #   if(trace){
      #     if(any(is.na(phi[j,]))){cat(paste("Phi for gene",j,"didn't converge in M step. Reinitializing with glm.nb().\n"))}
      #     if(any(is.na(coefs[j,]))){cat(paste("coefs for gene",j,"didn't converge in M step. Reinitializing with glm.nb().\n"))}
      #   }
      #   ###### re-initialize gene j with glm.nb() ######
      #   init_fit=init_fit=glm.init(j,y[j,],XX,k,offsets,wts,keep)
      #   coefs[j,]=init_fit$coefs_j
      #   phi[j,]=rep(init_fit$phi_g,k)
      #   phi_g[j]=init_fit$phi_g
      #   theta<-matrix(rep(0,times=k^2),nrow=k)
      #   for(c in 1:k){
      #     for(cc in 1:k){
      #       theta[c,cc]<-SCAD_soft_thresholding(coefs[j,c]-coefs[j,cc],lambda,alpha)
      #     }
      #   }
      #   theta_list[[j]] <- theta
      #   disc_ids[j]=any(theta_list[[j]]!=0)
      # }
      ### correct numerical inconsistency in gene-disp (fusion occurs, but slightly different cluster means)
      ### inconsistencies on the order of 1e-7 to 1e-8: won't add much bias to equate them
      # if(cl_phi==0){
      #   theta0=which(theta_list[[j]]==0,arr.ind=TRUE)
      #   fused_pairs=theta0[(theta0[,1]!=theta0[,2]),]
      #   fused_pairs=t(apply(fused_pairs,1,sort))
      #   fused_pairs[!duplicated(fused_pairs),]
      #   if(ncol(fused_pairs)!=0){
      #     for(fused_c in 1:nrow(fused_pairs)){                            # set coefs equal
      #       coefs[j,fused_pairs[fused_c,2]] = coefs[j,fused_pairs[fused_c,1]]
      #     }
      #   }
      # }

      # calculate LFCs (should just be 0 if k=1)
      LFCs[j,a] = max(coefs[j,1:k])-min(coefs[j,1:k])

      if(a>6){
        # set relative change in phi across 5 iterations
        if(cl_phi==1){
          diff_phi[a,j]=mean(abs(phi[j,]-phi_list[[a-5]][j,])/phi_list[[a-5]][j,])
        } else if(cl_phi==0){
          diff_phi[a,j]=abs(phi_g[j]-phi_list[[a-5]][j,1])/phi_list[[a-5]][j,1]
        }

        # if(is.infinite(diff_phi[a,j]) | is.na(diff_phi[a,j])){
        #   diff_phi[a,j]=1
        # }

        if(diff_phi[a,j]<0.01 & all(current_clusters==prev_clusters)){
          est_phi[j]=0
        } else{
          est_phi[j]=1
        }
        if(trace){
          cat(paste("#genes to continue est. phi next iter:",sum(est_phi),".\n"))
          cat(paste("Avg % diff in phi est (across 5 its) gene 1 = ",diff_phi[a,1],"\n"))
        }
      }
    }

    # all_temp_list[[a]] = temp_list
    # all_theta_list[[a]] = theta_list

    # Marker of all nondisc genes (T for disc, F for nondisc)
    if(trace){
      cat(paste("Disc genes:",sum(disc_ids),"of",g,"genes.\n"))
    }

    # save current objects in lists (to track)
    disc_ids_list[[a]] = disc_ids
    phi_list[[a]] <- phi
    #coefs_list[[a]] = coefs       # don't need coefs_list: save some time/memory


    if(a>1){
      #MSE_coefs[a]=mean(abs((coefs_list[[a]]-coefs_list[[a-1]])))
      MSE_coefs[a] = mean(abs(coefs-prev_coefs))
      if(trace){cat(paste("MSE_coef:",MSE_coefs[a],"\n"))}
      cl_agreement[a] = mean(current_clusters==prev_clusters)
    }

    prev_coefs = coefs    # store prev coefs for comparison in MSE_coefs next iteration

    # update on pi_hat, and UB & LB on pi
    for(c in 1:k){
      pi[c]=mean(wts[c,])
      if(pi[c]<1E-6){
        if(trace){warning(paste("cluster proportion", c, "close to 0"))}
        pi[c]=1E-6
      } # lowerbound for pi
      if(pi[c]>(1-1E-6)){
        if(trace){warning(paste("cluster proportion", c, "close to 1"))}
        pi[c]=(1-1E-6)
      } # upperbound for pi
    }

    # nb log(f_k(y_i))
    l<-matrix(rep(0,times=k*n),nrow=k)
    if(covars){
      covar_coefs = matrix(coefs[,-(1:k)],ncol=p)
      cov_eff = X %*% t(covar_coefs)         # n x g matrix of covariate effects
    } else {cov_eff=matrix(0,nrow=n,ncol=g)}
    for(i in 1:n){
      for(c in 1:k){
        # All genes. nondisc genes should be cancelled out in E step calculation
        if(cl_phi==1){
          l[c,i]<-sum(dnbinom(y[,i],size=1/phi[,c],mu=2^(coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))    # posterior log like, include size_factor of subj
        } else if(cl_phi==0){
          l[c,i]<-sum(dnbinom(y[,i],size=1/phi_g,mu=2^(coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))
        }
      }
    }

    # store and check Q function
    Q[a]<- (log(pi)%*%rowSums(wts)) + sum(wts*l)

    # break condition for EM
    if(a>n_mb){
      if(cl_agreement[a]==1){
        if(abs((Q[a]-Q[a-n_mb])/Q[a-n_mb])<EM_tol){
          # stop conditions: relative difference in Q within EM_tol and clusters stop changing
          # n_mb = g/mb_size, or number of minibatches that g can contain (rounded up/ceiling)
          # check Q function convergence across n_mb iterations. if mb_size = g/5, then check every 5 iters
          finalwts<-wts
          break
    }}}
    if(a==maxit_EM){
      finalwts<-wts
      if(trace){warning("Reached max iterations.")}
      DNC=1
      break
    }

    # E step
    start_E = as.numeric(Sys.time())
    Estep_fit=E_step(wts,l,pi,CEM,Tau)
    wts=Estep_fit$wts; keep=Estep_fit$keep; Tau=Estep_fit$Tau; CEM=Estep_fit$CEM
    end_E = as.numeric(Sys.time())
    if(trace){cat(paste("E Step Time Elapsed:",end_E-start_E,"seconds\n"))}

    # Diagnostics Tracking
    current_clusters = apply(wts,2,which.max)

    if(length(unique(current_clusters))<k & length(unique(prev_clusters))<k){
      finalwts=wts
      lower_K=TRUE
      warning(sprintf("EM iteration ended at iter%d, suboptimal order (choose higher K)",a))
      break
    }

    #print(current_clusters)
    if(trace){
      cat(paste("EM iter",a,"% of cls unchanged (from previous):",sum(current_clusters==prev_clusters)/n,"\n"))
      if(!is.null(true_clusters) & !is.null(true_clusters)){cat(paste("ARI =",adjustedRandIndex(true_clusters,current_clusters),"\n"))}
      cat(paste("Cluster proportions:",pi,"\n"))
      if(!is.null(true_disc)){
        if(sum(true_disc)==0){
          disc_gene=1
          cat("No discriminatory genes. Printing Gene1 instead\n")
        } else{ disc_gene = which(true_disc^2==1)[1] }
        if(sum(true_disc)==length(true_disc)){
          nondisc_gene=2
          cat("No nondiscriminatory genes. Printing Gene2 instead\n")
        } else{ nondisc_gene = which(true_disc^2==0)[1] }
        cat(paste("Disc Gene",disc_gene,": # of IRLS iterations used in M step:",nrow(temp_list[[disc_gene]][rowSums(temp_list[[disc_gene]])!=0,]),"\n"))
        cat(paste("coef:",coefs[disc_gene,],"\n"))
        if(cl_phi==1){
          cat(paste("phi:",phi[disc_gene,],"\n"))
        } else if(cl_phi==0){
          cat(paste("phi:",phi_g[disc_gene],"\n"))
        }
        cat(paste("Nondisc Gene",nondisc_gene,": # of IRLS iterations used in M step:",nrow(temp_list[[nondisc_gene]][rowSums(temp_list[[nondisc_gene]])!=0,]),"\n"))
        cat(paste("coef:",coefs[nondisc_gene,],"\n"))
        if(cl_phi==1){
          cat(paste("phi:",phi[nondisc_gene,],"\n"))
        } else if(cl_phi==0){
          cat(paste("phi:",phi_g[nondisc_gene],"\n"))
        }
      } else{
        cat(paste("Gene1: # of IRLS iterations used in M step:",nrow(temp_list[[1]][rowSums(temp_list[[1]])!=0,]),"\n"))
        cat(paste("coef:",coefs[1,],"\n"))
        if(cl_phi==1){
          cat(paste("phi:",phi[1,],"\n"))
        } else if(cl_phi==0){
          cat(paste("phi:",phi_g[1],"\n"))
        }
      }
      cat(paste("Samp1: PP:",wts[,1],"\n"))
      EMend = as.numeric(Sys.time())
      cat(paste("EM iter",a,"time elapsed:",EMend-EMstart,"seconds.\n"))
      cat("-------------------------------------\n")
    }

  }

  if(trace){cat("-------------------------------------\n")}
  num_warns=length(warnings())

  final_clusters = apply(finalwts,2,which.max)

  nondiscriminatory=rep(FALSE,times=g)
  if(lambda==0){
    m=rep(k,g) # if no penalty --> all disc, k params estimated per gene
  } else{
    m<-rep(0,times=g)
    for(j in 1:g){
      m[j]=length(unique(theta_list[[j]][1,]))
      # m_row=rep(0,k)
      # for(c in 1:k){
      #   m_row[c] <- sum(theta_list[[j]][c,]!=0) + 1         # of parameters estimated
      # }
      # m[j]=min(m_row)
      # if(any(coefs[j,] %in% c(-100,100))){
      #   m[j]=m[j]+(sum(coefs[j,] %in% c(-100,100))-1)
      # }
      if(m[j]==1){nondiscriminatory[j]=TRUE}
    }
  }

  num_est_coefs = sum(m)
  num_est_params =
    if(cl_phi==1){
      sum(m)+(k-1)+p*g + k*g                # p*g for covariates, sum(m) for coefs, k*g for cluster phis, k-1 for mixture proportions
    } else{ sum(m)+(k-1)+g+p*g }            # sum(m) for coef for each discriminatory clusters (cl_phi=1). sum(m) >= g
  # sum(m)+g for #coef+#phi (cl_phi=0)
  # (k-1) for mixture proportions

  log_L<-sum(apply(log(pi) + l, 2, logsumexpc))
  BIC = -2*log_L + log(n)*num_est_params

  if(lower_K){
    BIC=NA
    if(trace){
      print("Choose lower K. clusters identified: cls")
      print(unique(current_clusters))
    }
  }

  if(trace){
    cat(paste("total # coefs estimated =",num_est_coefs,"\n"))
    cat(paste("total # params estimated =",num_est_params,"\n"))
    cat(paste("-2log(L) =",-2*log_L,"\n"))
    cat(paste("log(n) =",log(n),"\n"))
    cat(paste("BIC =",BIC,"\n"))
    cat(paste("Tau =",Tau,"\n"))

    disc_stats=cbind(m,(!nondiscriminatory)^2,disc_ids^2)
    colnames(disc_stats) = c("#params","disc","disc_ids")
    write.table(head(disc_stats,n=10),quote=F)
    write.table(tail(disc_stats,n=10),quote=F)
    cat("-------------------------------------\n")
    cat("Coefs:\n")
    write.table(head(coefs,n=10),quote=F)
    write.table(tail(coefs,n=10),quote=F)
    cat("-------------------------------------\n")
    cat("Phi:\n")
    if(cl_phi==1){
      write.table(head(phi,n=10),quote=F)
      write.table(tail(phi,n=10),quote=F)
    } else if(cl_phi==0){
      write.table(head(phi_g,n=10),quote=F)
      write.table(tail(phi_g,n=10),quote=F)
    }
    cat("-------------------------------------\n")
  }

  end_time <- Sys.time()
  time_elap <- as.numeric(end_time)-as.numeric(start_time)

  if(cl_phi==0){
    phi = phi_g
  }

  result<-list(k=k,
               pi=pi,
               coefs=coefs,
               Q=Q[1:a],
               BIC=BIC,
               discriminatory=!(nondiscriminatory),
               init_clusters=init_cls,#init_coefs=init_coefs,init_phi=init_phi,
               clusters=final_clusters,
               phi=phi,
               #logL=log_L,
               wts=wts,
               time_elap=time_elap,
               lambda=lambda,
               alpha=alpha,
               size_factors=size_factors,norm_y=norm_y,DNC=DNC,n_its=a
  )
  return(result)
}

#' Prediction via FSCseq
#'
#' Performs prediction of cluster membership via trained model
#'
#' @param X (optional) design matrix of dimension n by p
#' @param fit FSCseq object
#' @param y_pred prediction/test count data matrix of dimension g by n_pred
#' @param size_factors_pred vector of length n_pred, size factors for prediction subjects
#'
#' @return list containing outputs
#' final_clusters: vector of length n of resulting clusters,
#' wts: k by n matrix of E step weights
#'
#' @export
FSCseq_predict <- function(X=NULL,fit,y_pred,size_factors_pred){
  # fit: Output of EM
  # y_pred: Data to perform prediction on
  # size_factors_pred: SF's of new data
  # offsets: Additional offsets per sample can be incorporated
  if(is.null(X)){
    cat("No covariates specified. Predicting on cluster-specific intercept-only model.\n")
  }else{
    cat("Predicting, adjusting for input X...\n")
    if (any(is.na(X))) {stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X")}
  }

  covars = !is.null(X)             # what should be done with missing values
  if(covars){
    p=ncol(X)
  } else{
    cat("No covariates specified")
    p=0
  }

  if(length(size_factors_pred) != ncol(y_pred)){stop("length of prediction size factors must be the same as the number of columns in prediction data")}

  # fit is the output object from the EM() function
  init_coefs=fit$coefs
  init_phi=fit$phi
  init_lambda=fit$lambda
  init_alpha=fit$alpha


  cl_phi=!is.null(dim(init_phi))  # dimension of phi is null when gene-wise (vector)

  init_size_factors = size_factors_pred
  offsets=log2(init_size_factors)
  n=ncol(y_pred)
  g=nrow(y_pred)
  k=fit$k

  # nb log(f_k(y_i))
  l<-matrix(0,nrow=k,ncol=n)
  for(i in 1:n){
    for(c in 1:k){
      if(covars){
        covar_coefs = matrix(init_coefs[,-(1:k)],ncol=p)
        cov_eff = X %*% t(covar_coefs)         # n x g matrix of covariate effects
      } else {cov_eff=matrix(0,nrow=n,ncol=g)}

      if(cl_phi){
        l[c,i]<-sum(dnbinom(y_pred[,i],size=1/init_phi[,c],mu=2^(init_coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))    # posterior log like, include size_factor of subj
      } else if(!cl_phi){
        l[c,i]<-sum(dnbinom(y_pred[,i],size=1/init_phi,mu=2^(init_coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))
      }
    }    # subtract out 0.1 that was added earlier
  }

  pi=fit$pi

  # E step
  # Estimate weights
  wts = matrix(0,nrow=k,ncol=n)
  logdenom = apply(log(pi) + l, 2,logsumexpc)
  for(c in 1:k){
    wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
  }

  final_clusters<-rep(0,times=n)
  for(i in 1:n){
    final_clusters[i]<-which.max(wts[,i])
  }

  return(list(clusters=final_clusters,wts=wts))
}
