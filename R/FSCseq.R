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

#' Theta estimation with penalty
#'
#' theta.ml() with small ridge penalty to stabilize estimates
#'
#' @param v a vector or matrix of numbers
#'
#' @return theta
#'
#'
#' @export
theta.ml2=function (y, mu, n = sum(weights), weights, limit = 10, eps = .Machine$double.eps^0.25,
                    trace = FALSE)
{
  lambda=1E-25
  score <- function(n, th, mu, y, w,lambda=1E-25) {sum(w * (digamma(th +
                                                                      y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +
                                                                                                                         th)/(mu + th))) + th*lambda}
  info <- function(n, th, mu, y, w,lambda=1E-25) {sum(w * (-trigamma(th +
                                                                       y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +
                                                                                                                            th)^2)) + lambda}
  if (inherits(y, "lm")) {
    mu <- y$fitted.values
    y <- if (is.null(y$y))
      mu + residuals(y)
    else y$y
  }
  if (missing(weights))
    weights <- rep(1, length(y))
  t0 <- n/sum(weights * (y/mu - 1)^2)
  it <- 0
  del <- 1
  if (trace)
    message(sprintf("theta.ml: iter %d 'theta = %f'",
                    it, signif(t0)), domain = NA)
  while ((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    if(trace){
      print(score(n, t0, mu, y, weights))
      print(info(n, t0,mu, y, weights))
    }
    del <- score(n, t0, mu, y, weights)/(i <- info(n, t0,
                                                   mu, y, weights))
    t0 <- t0 + del
    if (trace)
      message("theta.ml: iter", it, " theta =",
              signif(t0))
  }
  if (t0 < 0) {
    t0 <- 0
    warning("estimate truncated at zero")
    attr(t0, "warn") <- gettext("estimate truncated at zero")
  }
  if (it == limit) {
    warning("iteration limit reached")
    attr(t0, "warn") <- gettext("iteration limit reached")
  }
  test = sqrt(1/i)
  if(is.na(test)){return(NULL)}
  attr(t0, "SE") <- sqrt(1/i)
  t0
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
E_step<-function(wts,l,pi,CEM,Tau,PP_filt){
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
  if(!is.null(PP_filt)){
    keep = (wts > PP_filt)^2      # matrix of 0's and 1's, dimensions k x n
  }else{keep=matrix(1,nrow=nrow(wts),ncol=ncol(wts))}

  return(list(wts=wts,keep=keep,Tau=Tau,CEM=CEM))
}

#' Estimates dispersion parameter phi
#'
#' Uses theta.ml2 and theta.mm functions to tryCatch
#' any unstable estimations
#'
#' @param y vector: gene read counts for one fixed gene
#' @param mu vector: mean. must be same length as y
#' @param dfr numeric: degrees of freedom for theta.mm()
#' @param weights vector: must be same length as y
#' @param limit numeric: max number of iterations for theta.ml() and theta.mm()
#' @param trace logical: True/False to show trace output from theta.ml
#'
#' @return phi_g: estimate of phi
#'
#' @importFrom MASS theta.mm
#'
#' @export
phi.ml=function(y,mu,dfr,weights,limit,trace){
  theta_g = NULL
  try(theta_g <- theta.ml2(y=y, mu=mu, weights=weights, limit=limit, trace=trace),silent=TRUE)

  if(is.null(theta_g)){
    try(theta_g <- MASS::theta.mm(y=y, mu=mu, dfr=dfr, weights=weights, limit=limit))
  } else{
    if(theta_g<0.1){
      # if theta_g was estimated by theta.ml2() but phi_g >10, i.e. theta_g < 0.1 (overdispersion too large: unstable estimate) --> re-estimate with mm
      try(theta_g <- MASS::theta.mm(y=y, mu=mu, dfr=dfr, weights=weights, limit=limit))
    }
  }
  if(is.null(theta_g)){
    theta_g=Inf      # set to Poisson (phi=0) if theta_g didn't converge in either method
  }
  phi_g = 1/theta_g
  return(phi_g)
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
#'
#' @export
glm.init=function(j,y_j,XX,k,offsets,wts,keep){
  ids = c(t(keep==1))  # PP filtering

  init.fit = stats::glm(as.integer(rep(y_j,k))[ids]~0+XX[ids,]+offset(rep(offsets,k)[ids]),family=poisson(link="log"),weights=c(t(wts))[ids])
  coefs_j = log2(exp(init.fit$coefficients))         # change to log2 scale

  coefs_j[is.na(coefs_j)] = 0
  coefs_j[which(coefs_j< -50)] = -50
  coefs_j[which(coefs_j>50)] = 50
  mu = 2^(XX %*% coefs_j + offsets)

  phi_g=phi.ml(y=as.integer(rep(y_j,k))[ids],
               mu=mu[ids],
               dfr=sum(ids)-1,
               weights=c(t(wts))[ids],
               limit=25,trace=F)

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
#'
#' @export
glm.init_par=function(j){
  ids = c(t(keep==1))  # PP filtering

  init.fit = stats::glm(as.integer(rep(y[j,],k))[ids]~0+XX[ids,]+offset(rep(offsets,k)[ids]),family=poisson(link="log"),weights=c(t(wts))[ids])
  coefs_j = log2(exp(init.fit$coefficients))         # change to log2 scale

  coefs_j[is.na(coefs_j)] = 0
  coefs_j[which(coefs_j< -50)] = -50
  coefs_j[which(coefs_j>50)] = 50

  mu = 2^(XX %*% coefs_j + offsets)

  phi_g=phi.ml(y=as.integer(rep(y[j,],k))[ids],
               mu=mu[ids],
               dfr=sum(ids)-1,
               weights=c(t(wts))[ids],
               limit=25,trace=F)

  results=list(coefs_j=coefs_j,
               phi_g=phi_g)
  return(results)
}

#' Parallelized M step (parLApply)
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
               cl_phi=cl_phi,est_covar=est_covar[j],
               lambda=lambda,alpha=alpha,
               IRLS_tol=IRLS_tol,maxit_IRLS=maxit_IRLS,
               CDA_tol=CDA_tol,maxit_CDA=maxit_CDA)

  if(est_phi[j]==1){
    ids = (c(t(keep))==1)
    mu = 2^(XX %*% res$coefs_j + offsets)

    phi_g_temp=phi.ml(y=as.integer(rep(y[j,],k))[ids],
                      mu=mu[ids],
                      dfr=sum(ids)-1,
                      weights=c(t(wts))[ids],
                      limit=25,trace=F)
    res$phi_j = rep(phi_g_temp,k)
  } else{ res$phi_j=phi[j,]}
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
                       lambda, alpha, IRLS_tol, CDA_tol, maxit_IRLS, maxit_CDA,
                       Tau, disc_ids_list, par_X){
  if(Tau<=1 & a>6){
    if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){
      res=par_X[[j]]
      return(res)
    }}
  res = M_step(X=XX, y_j=as.numeric(rep(y[j,],k)), p=p, j=j, a=a, k=k,
               all_wts=wts, keep=c(t(keep)), offset=rep(offsets,k),
               theta=theta_list[[j]],coefs_j=coefs[j,],phi_j=phi[j,],
               cl_phi=cl_phi,est_covar=est_covar[j],
               lambda=lambda,alpha=alpha,
               IRLS_tol=IRLS_tol,CDA_tol=CDA_tol,
               maxit_IRLS=maxit_IRLS,maxit_CDA=maxit_CDA)

  if(est_phi[j]==1){
    ids = (c(t(keep))==1)
    mu = 2^(XX %*% res$coefs_j + offsets)

    phi_g_temp=phi.ml(y=as.integer(rep(y[j,],k))[ids],
                      mu=mu[ids],
                      dfr=sum(ids)-1,
                      weights=c(t(wts))[ids],
                      limit=25,trace=F)
    res$phi_j = rep(phi_g_temp,k)
  } else{res$phi_j = phi[j,]}
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
#' @param n_rinits integer, number of additional random initializations to be searched (default 50 for EM, 10 for CEM)
#' @param maxit_inits integer, maximum number of iterations for each initialization search (default 15 for EM, or until temperature anneals down to below 2 for CEM)
#' @param maxit_EM integer, maximum number of iterations for full EM/CEM run (default 100)
#' @param maxit_IRLS integer, maximum number of iterations for IRLS algorithm, in M step (default 50)
#' @param maxit_CDA integer, maximum number of iterations for CDA loop (default 50)
#' @param EM_tol numeric, tolerance of convergence for EM/CEM, default is 1E-6
#' @param IRLS_tol numeric, tolerance of convergence for IRLS, default is 1E-4
#' @param CDA_tol numeric, tolerance of convergence for IRLS, default is 1E-4
#' @param method string, either "EM" (default) or "CEM"
#' @param init_temp numeric, default for CEM: init_temp = nrow(y), i.e. number of genes. temp=1 for EM
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
#' @param mb_size minibatch size: # of genes to include per M step iteration
#' @param PP_filt numeric between (0,1), threshold on PP for sample/cl to be included in M step estimation. Default is 1e-3
#'
#' @return list containing outputs from EM_run() function
#'
#' @export
FSCseq<-function(ncores=1,X=NULL, y, k,
                 lambda=0,alpha=0,
                 size_factors,norm_y,
                 true_clusters=NULL, true_disc=NULL,
                 init_parms=FALSE,init_coefs=NULL,init_phi=NULL,init_cls=NULL,init_wts=NULL,
                 n_rinits=if(method=="EM"){20}else if(method=="CEM"){1},         # fewer searches for CEM to minimize computational cost
                 maxit_inits=if(method=="EM"){15}else{ceiling(log(2/nrow(y))/log(0.9))}, # for CEM, tolerates end temp (Tau) of 2 at end of initialization
                 maxit_EM=100,maxit_IRLS=50,maxit_CDA=50,EM_tol=1E-6,IRLS_tol=1E-4,CDA_tol=1E-4,
                 disp="gene", # disp is commented out. just left for simulations
                 method="EM",init_temp=nrow(y),
                 trace=F,trace.file=NULL,
                 mb_size=NULL,PP_filt=1e-3){
  # disp = "gene". disp="cluster" turned off (for now)
  disp="gene"
  # init_method = c("max","emaEM","BIA"). Turned off --> max: minimum BIC
  # optim_method = c("direct","NR","GD"). Direct, Newton-Raphson, or Gradient descent (fixed step size of 2)). All equivalent
  # BIC_penalty = c("n","n2","ng","ng2,"en","en2","eng","eng2")     # turned off --> set to "ng"
  # gamma=NULL: EBIC hyperparameter (EBIC turned off).
  # set PP_filt = NULL or 0 < PP_filt < 1e-50: turn off PP_filt in M step (not recommended due to computation time)


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

  start_FSC = as.numeric(Sys.time())

  n = ncol(y)
  g = nrow(y)
  p<-if(is.null(X)){0}else{ncol(X)}         # number of covariates
  if(is.null(init_parms)){cat(paste(method,"model with",disp,"level dispersions specified.\n"))}

  # check inputs #
  if(is.null(X)){
    if(is.null(init_parms)){cat("No covariates specified. Running cluster-specific intercept-only model.\n")}
  } else{
    if (class(X) != "matrix") {
      tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
      if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
    }
    if (storage.mode(X)=="integer") storage.mode(X) <- "double"
    if (ncol(y) != nrow(X)) stop("X and y do not have the same number of observations")
    if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y")
  }
  if(alpha < 0 | alpha >= 1){stop("alpha must be 0 <= alpha < 1")}
  if(lambda<0){stop("lambda must be greater than 0")}
  if(method=="EM"){CEM=F} else if(method=="CEM"){CEM=T} else{stop("method must be 'EM' or 'CEM'.")}
  if(CEM){if(init_temp<0){stop('init_temp must be greater than 0')}}
  if(length(size_factors) != n){stop("size_factors must be of length n")}
  if(any(size_factors <= 0)){stop("size_factors must be > 0")}
  if(any(y < 0) | any(norm_y < 0)){stop("y and norm_y must all be >= 0")}
  if(ncores%%1!=0 | k%%1!=0 | n_rinits%%1!=0 | maxit_inits%%1!=0 | maxit_EM%%1!=0 | maxit_IRLS%%1!=0 | maxit_CDA%%1!=0){stop('ncores, k, n_rinits, maxit_inits, maxit_EM, maxit_IRLS, maxit_CDA must all be positive integers')}
  if(EM_tol<0 | IRLS_tol<0 | CDA_tol<0 | EM_tol>1 | IRLS_tol>1 | CDA_tol>1){stop('EM_tol, IRLS_tol, and CDA_tol must be between 0 and 1')}
  if(PP_filt<0 | PP_filt>=1){stop('PP_filt must be 0 <= PP_filt < 1')}
  if(!is.null(mb_size)){if(mb_size>g | mb_size<=0){stop('mb_size must be 0 < mb_size <= g')}}
  if(init_parms){if(is.null(init_coefs) | is.null(init_phi)){stop("init_coefs (gxk matrix) and init_phi (vector of length g) must be input if init_parms=TRUE")}}
  if(!is.null(init_cls)){
    if(length(unique(init_cls))>k){stop("Too many cluster levels in init_cls (less than specified value of k)")}
    if(length(init_cls) != n){stop("Length of initial clusters not equal to n")}
  }
  if(!is.null(init_wts)){
    if (class(init_wts) != "matrix") {
      tmp <- try(init_wts <- model.matrix(~0+., data=init_wts), silent=TRUE)
      if (class(tmp)[1] == "try-error") stop("init_wts must be a matrix or able to be coerced to a matrix")
    }
    if(nrow(init_wts)!=k | ncol(init_wts)!=n){stop('init_wts must be matrix of k rows and n columns')}
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

  n_mb = if(!is.null(mb_size)){ceiling(g/mb_size)     # number of minibatches in g. default: 1. experiment with 5 (mb_size = g/5)
  }else{1}
  maxit_EM = maxit_EM*n_mb
  maxit_inits = maxit_inits*n_mb

  if(k==1){
    init_cls=rep(1,n)
  }

  # # if initial clustering is provided (turned off: can now let fewer initial clusters than k)
  # if(!is.null(init_cls)){if(length(unique(init_cls))<k){ # if previous clustering labels have fewer than k cls (warm starts)
  #   if(trace){
  #     if(!is.null(trace.file)){
  #       sink()
  #     }
  #   }
  #   return(list(k=k,pi=NA,coefs=init_coefs,Q=NA,BIC=NA,discriminatory=NA,init_clusters=init_cls,clusters=init_cls,
  #               phi=init_phi,wts=NA,time_elap=0,lambda=lambda,alpha=alpha,size_factors=size_factors,norm_y=norm_y,DNC=T,n_its=0,lower_K=T))     # NA values for BIC/n_its/time_elap. don't run FSCseq
  # }}

  # if initial clustering is not provided
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
      # init clusters weren't specified, init_coefs/phi were not either
      fit = EM_run(ncores,X,y,k,lambda,alpha,size_factors,norm_y,true_clusters,true_disc,
                   init_parms=F,disp=disp,
                   init_cls=all_init_cls[,i], CEM=CEM,init_Tau=init_Tau,
                   maxit_EM=maxit_inits,maxit_IRLS=maxit_IRLS,maxit_CDA=maxit_CDA,
                   EM_tol=EM_tol,IRLS_tol=IRLS_tol,CDA_tol=CDA_tol,
                   trace=trace,mb_size=mb_size,PP_filt=PP_filt)

      all_fits [[i]] = fit
      init_cls_BIC[i] <- fit$BIC
    }

    ## Selecting maximum BIC model ##
    fit_id = which.min(init_cls_BIC)[1]
    init_cls = all_fits[[fit_id]]$clusters
    init_parms = T
    init_coefs = all_fits[[fit_id]]$coefs
    init_phi = all_fits[[fit_id]]$phi
    if(trace){
      cat("FINAL INITIALIZATION:\n")
      cat(paste(colnames(all_init_cls)[fit_id],"\n"))
    }
    # if(length(unique(init_cls))<k){
    #   warning("best initialization yielded smaller k. returning result from best initialization (EM to full convergence not run)")
    #   if(trace){
    #     if(!is.null(trace.file)){
    #       sink()
    #     }
    #   }
    #   return(all_fits[[fit_id]])
    # }

    # if(init_method %in% c("emaEM","BIA")){
    #   #### Consensus methods ####
    #   # Create Z matrices and calculate weights
    #   list_Z = list()
    #   w_BIC_emaEM = rep(NA,n_inits)
    #   w_BIC_emaEM_logdenom = logsumexpc(-0.5*init_cls_BIC)
    #
    #   max_BIC = max(init_cls_BIC,na.rm=T)
    #   w_BIC_BIA = rep(NA,n_inits)
    #   w_BIC_BIA_logdenom = logsumexpc(-0.5*(init_cls_BIC-max_BIC))
    #   for(m in 1:n_inits){
    #     list_Z[[m]] = matrix(0,nrow=k,ncol=n)
    #     for(c in 1:k){
    #       list_Z[[m]][c,]=(all_init_cls[,m]==c)^2
    #     }
    #     # w_BIC_emaEM[m]=if(is.na(init_cls_BIC[m])){0} else{exp(-0.5*init_cls_BIC[m]-w_BIC_emaEM_logdenom)}
    #     # w_BIC_BIA[m]=if(is.na(init_cls_BIC[m])){0} else{exp(-0.5*(init_cls_BIC[m]-max_BIC)-w_BIC_BIA_logdenom)}
    #     ##### EXPERIMENT #####
    #     temp=1     # temp=1 is the same as regular emaEM and BIA
    #     w_BIC_emaEM[m]=if(is.na(init_cls_BIC[m])){0} else{exp((1/temp)*(-0.5)*init_cls_BIC[m]-logsumexpc((-0.5*init_cls_BIC)*(1/temp)))}
    #     w_BIC_BIA[m]=if(is.na(init_cls_BIC[m])){0} else{exp((1/temp)*(-0.5)*(init_cls_BIC[m]-max_BIC)-logsumexpc((-0.5*(init_cls_BIC-max_BIC))*(1/temp)))}
    #   }
    #
    #   ## emaEM method (Michael & Melnykov 2016) ##
    #   if(init_method == "emaEM"){
    #     # create y_m, A_m for all inits --> A --> HC --> replace init_cls
    #     A = matrix(0,n,n)
    #     for(m in 1:n_inits){
    #       # construct A_m matrix for all inits (y_m omitted to conserve memory. can obtain if you take out w_BIC_emaEM[m] = 1)
    #       # sum directly here to save memory
    #       A_m=matrix(NA,n,n)
    #       for(i in 1:n){
    #         for(ii in 1:n){
    #           A_m[i,ii] = w_BIC_emaEM[m]*all(list_Z[[m]][,i] == list_Z[[m]][,ii])^2   # if all wts for sample i and sample ii are the same --> y_m[i,ii]=1. otherwise, 0
    #         }
    #       }
    #       A = A + A_m
    #     }
    #     model = hclust(as.dist(A),method="average")
    #     init_cls <- cutree(hclust(as.dist(A),method="average"),k=k)
    #   }
    #
    #   ## BIA method (Hagan & White 2018) ##
    #   if(init_method == "BIA"){
    #     # create Z* = sum(w*Z) --> replace init_wts
    #     init_wts=matrix(0,nrow=k,ncol=n)
    #     for(m in 1:n_inits){
    #       init_wts = init_wts + w_BIC_BIA[m]*list_Z[[m]]
    #     }
    #   }
    # }

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
                 maxit_EM=maxit_EM,maxit_IRLS=maxit_IRLS,maxit_CDA=maxit_CDA,
                 EM_tol=EM_tol,IRLS_tol=IRLS_tol,CDA_tol=CDA_tol,
                 trace=trace,mb_size=mb_size,PP_filt=PP_filt)

  end_FSC = as.numeric(Sys.time())
  results$total_time_elap = end_FSC-start_FSC

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
#' @param maxit_IRLS integer, maximum number of iterations for IRLS loop, in M step (default 50)
#' @param maxit_CDA integer, maximum number of iterations for CDA loop (default is 50)
#' @param EM_tol numeric, tolerance of convergence for EM/CEM, default is 1E-6
#' @param IRLS_tol numeric, tolerance of convergence for IRLS, default is 1E-4
#' @param CDA_tol numeric, tolerance of convergence for CDA, default is 1E-4
#' @param disp string, either "gene" (default) or "cluster"
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
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
                   maxit_EM=100, maxit_IRLS = 50,maxit_CDA=50,EM_tol = 1E-6,IRLS_tol = 1E-4,CDA_tol=1E-4,disp,trace=F,
                   mb_size=NULL,PP_filt){

  start_time <- Sys.time()

  n<-ncol(y)         # number of samples
  g<-nrow(y)         # number of genes

  p<-if(is.null(X)){0}else{ncol(X)}         # number of covariates

  n_mb = if(!is.null(mb_size)){ceiling(g/mb_size)     # number of minibatches in g. default: 1. experiment with 5 (mb_size = g/5)
  }else{1}

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

  # if PP_filt is not set to NULL --> keep only obs/cl in M step > PP_filt threshold
  if(!is.null(PP_filt)){
    keep = (wts > PP_filt)^2
  }else{keep = matrix(1,nrow=nrow(wts),ncol=ncol(wts))}

  # if any cluster has 0 samples to estimate with, set keep to all samples --> weight them w/ small 1e-50 weights
  if(k>1){if(any(rowSums(keep)==0)){
    keep[rowSums(keep)==0,]=1
  }}

  ## Manual turn-off of estimating phi/covariates
  # if(init_parms){
  #   est_phi=rep(0,g)
  #   est_covar = rep(0,g)
  # }
  # all_temp_list = list()
  # all_theta_list = list()

  # lower_K=FALSE         # tracks when number of a mixture component --> 0

  cl_agreement = rep(NA,maxit_EM)
  par_X=rep(list(list()),g)         # store M_step results

  nits_IRLS=rep(0,g)
  nits_CDA=rep(0,g)

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
            clusterEvalQ(cl=clust,library(MASS))
            clusterEvalQ(cl=clust,library(stats))
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
        #disc_ids[j]=any(theta_list[[j]]!=0)
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

    # M step
    Mstart=as.numeric(Sys.time())

    # minibatching starts at iteration 5 (let EM stabilize first: turned off. see if this is stable)
    # if(a<5){
    #   mb_genes = 1:g
    # } else{
    if(!is.null(mb_size)){
      mb_genes = sample(1:g,mb_size,replace=F)
    } else{ mb_genes=1:g}
    # }

    if(ncores>1){
      # M_step parallelized across ncores
      if(Sys.info()[['sysname']]=="Windows"){
        ## use parLapply for Windows ##
        clust = makeCluster(ncores)
        clusterEvalQ(cl=clust,library(FSCseq))
        clusterExport(cl=clust,varlist=c("XX","y","p","a","k","wts","keep","offsets",
                                         "theta_list","coefs","phi","cl_phi","est_phi","est_covar",
                                         "lambda","alpha","IRLS_tol","CDA_tol","maxit_IRLS","maxit_CDA",
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
                      lambda, alpha, IRLS_tol, CDA_tol, maxit_IRLS, maxit_CDA,
                      Tau, disc_ids_list, par_X)
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
                             cl_phi=cl_phi,est_covar=est_covar[j],
                             lambda=lambda,alpha=alpha,
                             IRLS_tol=IRLS_tol,maxit_IRLS=maxit_IRLS,
                             CDA_tol=CDA_tol,maxit_CDA=maxit_CDA)
        if(est_phi[j]==1){
          ids = (c(t(keep))==1)
          mu = 2^(XX %*% par_X[[j]]$coefs_j + offsets)

          phi_g_temp=phi.ml(y=as.integer(rep(y[j,],k))[ids],
                            mu=mu[ids],
                            dfr=sum(ids)-1,
                            weights=c(t(wts))[ids],
                            limit=25,trace=F)
          par_X[[j]]$phi_j = rep(phi_g_temp,k)
        } else{par_X[[j]]$phi_j = phi[j,]}
      }
    }
    Mend=as.numeric(Sys.time())
    if(trace){cat(paste("M Step Time Elapsed:",Mend-Mstart,"seconds.\n"))}

    for(j in mb_genes){
      if(Tau<=1 & a>6){if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){next}}
      coefs[j,] <- par_X[[j]]$coefs_j
      theta_list[[j]] <- par_X[[j]]$theta_j
      # # correct for small computational/numerical inconsistencies: if theta = 0, set coefs to be exactly equal (fixed in Rcpp)
      # if( length(unique(coefs[j,1:k])) != length(unique(theta_list[[j]][1,])) ){
      #   fused_mean=rep(NA,k)    # calculate first. fix, then replace --> force all fused coefs to be equal
      #   for(c in 1:k){
      #     fused_mean[c] = mean(coefs[j,theta_list[[j]][,c]==0])
      #   }
      #   for(c in 1:k){
      #     coefs[j,c] = fused_mean[c]
      #   }
      # }
      disc_ids[j]=any(theta_list[[j]]!=0)
      temp_list[[j]] <- if(p>0){cbind(par_X[[j]]$temp_beta, par_X[[j]]$temp_gamma)}else{par_X[[j]]$temp_beta}
      if(cl_phi==1){
        phi[j,] <- par_X[[j]]$phi_j
      } else if(cl_phi==0){
        phi_g[j] <- (par_X[[j]]$phi_j)[1]
        phi[j,] = rep(phi_g[j],k)
      }
      nits_IRLS[j]=par_X[[j]]$nits_IRLS
      nits_CDA[j]=par_X[[j]]$nits_CDA

      # calculate LFCs (should just be 0 if k=1)
      LFCs[j,a] = max(coefs[j,1:k])-min(coefs[j,1:k])

      if(a>6){
        # set relative change in phi across 5 iterations
        if(cl_phi==1){
          diff_phi[a,j]=mean(abs(phi[j,]-phi_list[[a-5]][j,])/phi_list[[a-5]][j,])
        } else if(cl_phi==0){
          diff_phi[a,j]=abs(phi_g[j]-phi_list[[a-5]][j,1])/phi_list[[a-5]][j,1]
        }

        if(diff_phi[a,j]<0.01 & all(current_clusters==prev_clusters)){
          est_phi[j]=0
        } else{
          est_phi[j]=1
        }
      }
    }
    if(trace){
      cat(paste("#genes to continue est. phi next iter:",sum(est_phi),".\n"))
      cat(paste("Avg % diff in phi est (across 5 its) = ",mean(diff_phi[a,]),"\n"))
    }


    if(trace){cat(paste("Average # IRLS iterations:",mean(nits_IRLS[mb_genes]),"\n"))}
    if(trace){cat(paste("Average # CDA iterations:",mean(nits_CDA[mb_genes]),"\n"))}


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

    # nb log(f_k(y_i))
    l<-matrix(rep(0,times=k*n),nrow=k)
    if(covars){
      covar_coefs = matrix(coefs[,(k+1):(k+p)],ncol=p)
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
      if(cl_agreement[a]==1 & Tau <= 10){
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
    Estep_fit=E_step(wts,l,pi,CEM,Tau,PP_filt)
    wts=Estep_fit$wts; keep=Estep_fit$keep; Tau=Estep_fit$Tau; CEM=Estep_fit$CEM
    end_E = as.numeric(Sys.time())
    if(trace){cat(paste("E Step Time Elapsed:",end_E-start_E,"seconds\n"))}

    # if any cluster has 0 samples to estimate with, set keep to all samples --> weight them w/ small 1e-50 weights
    if(k>1){if(any(rowSums(keep)==0)){
      keep[rowSums(keep)==0,]=1
    }}

    # Diagnostics Tracking
    current_clusters = apply(wts,2,which.max)

    # if all PP's for a particular cluster is < PP_filt (threshold for keep).
    # shouldn't happen now: if any rowSums(keep) == 0, those rows are set to 1 at beginning of EM
    # if(any(rowSums(keep)==0)){
    #   lower_K=TRUE
    #   finalwts=wts
    #   pi[rowSums(keep)==0] = 1e-50      # this mixt. proportion --> 0
    #   warning(sprintf("No samples in a cluster %dth E step",a))
    #   break
    # }

    #print(current_clusters)
    if(trace){
      cat(paste("EM iter",a,"% of samples whose cls unchanged (from previous):",mean(current_clusters==prev_clusters),"\n"))
      if(!is.null(true_clusters)){cat(paste("ARI =",adjustedRandIndex(true_clusters,current_clusters),"\n"))}
      cat(paste("Cluster proportions:",pi,"\n"))

      # if true_disc gene list is input
      if(!is.null(true_disc)){

        # pick a disc gene
        if(sum(true_disc)==0){
          disc_gene=mb_genes[1]
          cat(sprintf("No discriminatory genes. Printing Gene%d instead\n",disc_gene))
        } else{ disc_gene = mb_genes[true_disc[mb_genes]][1] }      # select first of minibatch genes that is disc (would be random each time)

        # pick a nondisc gene
        if(sum(true_disc)==length(true_disc)){
          nondisc_gene=mb_genes[2]
          cat(sprintf("No nondiscriminatory genes. Printing Gene%d instead\n",nondisc_gene))
        } else{ nondisc_gene = mb_genes[!true_disc[mb_genes]][1] }   # select first of minibatch genes that is nondisc (would be random each time)

        # print number of IRLS iters, coefs, and phi for picked disc gene
        cat(paste("Disc Gene",disc_gene,": # of IRLS iterations used in M step:",nrow(temp_list[[disc_gene]][rowSums(temp_list[[disc_gene]])!=0,]),"\n"))
        cat(paste("coef:",coefs[disc_gene,],"\n"))
        cat(paste("phi:",phi[disc_gene,],"\n"))

        # print number of IRLS iters, coefs, and phi for picked nondisc gene
        cat(paste("Nondisc Gene",nondisc_gene,": # of IRLS iterations used in M step:",nrow(temp_list[[nondisc_gene]][rowSums(temp_list[[nondisc_gene]])!=0,]),"\n"))
        cat(paste("coef:",coefs[nondisc_gene,],"\n"))
        cat(paste("phi:",phi[nondisc_gene,],"\n"))

      } else{
        # if true_disc not specified, just track one gene (first one tincluded in minibatch for that M step iter)
        sel_gene = mb_genes[1]
        cat(paste(sprintf("Gene%d: # of IRLS iterations used in M step:",sel_gene),nrow(temp_list[[sel_gene]][rowSums(temp_list[[sel_gene]])!=0,]),"\n"))
        cat(paste("coef:",coefs[sel_gene,],"\n"))
        cat(paste("phi:",phi[sel_gene,],"\n"))
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

  # if(BIC_penalty %in% c("n","n2","en","en2")){eff_n = n}else if(BIC_penalty %in% c("ng","ng2","eng","eng2")){eff_n = n*g}

  # # add log(1+g*n_eff) for all log2 baseline estimates (betas)
  # BIC_penalty_term2n=0
  # BIC_penalty_term2ng=0
  # for(j in 1:g){
  #   unique_thetas = unique(theta_list[[j]][1,])  # also equal to the number of clusters for that particular gene.
  #   for(i in 1:length(unique_thetas)){
  #     cls_ids = which(theta_list[[j]][1,]==unique_thetas[i])      # returns all cluster indices corr to unique theta value
  #     n_cls_ids = sum(final_clusters %in% cls_ids)                # returns number of samples in that cluster
  #     BIC_penalty_term2ng = BIC_penalty_term2ng + log(1 + g*n_cls_ids)  # adds to penalty term log(1+g*(#samples involved in estimating that cl's beta))
  #     BIC_penalty_term2n = BIC_penalty_term2n + log(1 + n_cls_ids)  # adds to penalty term log(1+g*(#samples involved in estimating that cl's beta))
  #   }
  # }
  # # add in mixture proportions, gene disps, and cov effects
  # BIC_penalty_term2ng = BIC_penalty_term2ng + ((k-1) + (p+1)*g)*log(g*n)
  # BIC_penalty_term2n = BIC_penalty_term2n + ((k-1) + (p+1)*g)*log(n)

  ## for straight BIC approx log(n) or log(ng) (NA for the '2's)
  # BIC_penalty_term1 = if(BIC_penalty %in% c("n","en")){ log(n)*num_est_params }else if(BIC_penalty %in% c("ng","eng")){ log(n*g)*num_est_params } else{NA}

  # P = g*(k+p+1) + (k-1)
  #   kappa = log(P)/log(eff_n)
  #   gamma = 1 - 1/(2*kappa)
  # log_choose=function(P,j){
  #   # returns log(choose(P,j)), avoiding overflow issues
  #   if(j<P){
  #     return( sum(log(seq(j+1,P,1))) - sum(log(seq(1,P-j,1))) )
  #   } else if(j==P){
  #     return( 0 )
  #   }
  # }
  # eBIC_term = 2*gamma*log_choose(P,num_est_params)
  # BIC_penalty_term = if(BIC_penalty %in% c("n","ng","en","eng")){
  #       BIC_penalty_term1
  #     }else if(BIC_penalty %in% c("ng2","eng2")){
  #       BIC_penalty_term2ng
  #     }else if(BIC_penalty %in% c("n2","en2")){
  #       BIC_penalty_term2n
  #     }

  # BIC = if(BIC_penalty %in% c("n","n2","ng","ng2")){
  #   -2*log_L + BIC_penalty_term
  # } else if(BIC_penalty %in% c("en","en2","eng","eng2")){
  #   -2*log_L + BIC_penalty_term + eBIC_term
  # }
  # BIC_n = -2*log_L + log(n)*num_est_params
  # BIC_n2 = -2*log_L + BIC_penalty_term2n
  # BIC_ng = -2*log_L + log(n*g)*num_est_params
  # BIC_ng2 = -2*log_L + BIC_penalty_term2ng
  # eBIC_n = -2*log_L + log(n)*num_est_params + eBIC_term
  # eBIC_n2 = -2*log_L + BIC_penalty_term2n + eBIC_term
  # eBIC_ng = -2*log_L + log(n*g)*num_est_params + eBIC_term
  # eBIC_ng2 = -2*log_L + BIC_penalty_term2ng + eBIC_term

  eff_n=n*g
  BIC_penalty_term = log(eff_n)*num_est_params
  BIC = -2*log_L + BIC_penalty_term


  # if(lower_K){
  #   stop("something's wrong: this shouldn't be happening anymore")
  #   # print("K not optimal. clusters identified: cls")
  #   # print(unique(current_clusters))
  #   # k=length(unique(current_clusters))
  #   # BIC=Inf
  # }

  if(trace){
    # cat(paste("log(n) =",log(n),"\n"))
    # cat(paste("sum(log(1+n_eff)) =",BIC_penalty_term2n,"\n"))
    # cat(paste("sum(log(1+n_eff*g)) =",BIC_penalty_term2ng,"\n"))
    # cat(paste("BIC(n,n2,ng,ng2,en,en2,eng,eng2)=\n",BIC_n,"\n",BIC_n2,"\n",BIC_ng,"\n",BIC_ng2,"\n",eBIC_n,"\n",eBIC_n2,"\n",eBIC_ng,"\n",eBIC_ng2,"\n"))
    cat(paste("total # coefs estimated =",num_est_coefs,"\n"))
    cat(paste("total # params estimated =",num_est_params,"\n"))
    cat(paste("-2log(L) =",-2*log_L,"\n"))
    cat(paste("log(n*g) =",log(n*g),"\n"))
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
               # BIC_n=BIC_n, BIC_ng=BIC_ng, BIC_ng2=BIC_ng2,
               # eBIC_n=eBIC_n, eBIC_ng=eBIC_ng, eBIC_ng2=eBIC_ng2, eBIC_term=eBIC_term, gamma=gamma,
               discriminatory=!(nondiscriminatory),
               init_clusters=init_cls,#init_coefs=init_coefs,init_phi=init_phi,
               clusters=final_clusters,
               phi=phi,num_est_params=num_est_params,m=m,
               log_L=log_L,
               wts=wts,
               time_elap=time_elap,
               lambda=lambda,
               alpha=alpha,
               size_factors=size_factors,#norm_y=norm_y,
               DNC=DNC,n_its=a#,lower_K=lower_K
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


  if(length(size_factors_pred) != ncol(y_pred)){stop("length of prediction size factors must be the same as the number of columns in prediction data")}

  idx=fit$discriminatory
  if(sum(idx)==0){warning('No disc genes. Using all genes');idx=rep(TRUE,nrow(y_pred))}
  n=ncol(y_pred)
  y_pred=matrix(y_pred[idx,],ncol=n)      # subset to just disc genes found by FSCseq run
  # g=nrow(y_pred)
  g=sum(idx)
  pi=fit$pi
  k=length(pi)

  cl_phi=!is.null(dim(fit$phi))  # dimension of phi is null when gene-wise (vector)

  covars = !is.null(X)
  if(covars){
    p=ncol(X)
    init_coefs=matrix(fit$coefs[idx,],nrow=sum(idx))
  } else{
    cat("No covariates specified")
    p=0
    init_coefs=matrix(fit$coefs[idx,1:k],nrow=sum(idx))
  }

  # fit is the output object from the EM() function

  # # if EM fit lost mixture components (turned off now. can lose mixt component)
  # existing_cls = unique(fit$clusters)[order(unique(fit$clusters))]      # accounts for if lower K was selected
  # k=length(existing_cls)
  # # store coefs of just existing clusters (existing_cls, from FSCseq fit), of just disc genes (idx)
  # init_coefs=if(p>0){    # assumes covariate effects are last p columns of coefs
  #   matrix(fit$coefs[idx,c(existing_cls,(ncol(fit$coefs)-(p-1)):ncol(fit$coefs))],nrow=g,ncol=k+p)
  # } else{matrix(fit$coefs[idx,existing_cls],nrow=g,ncol=k)}
  # pi=fit$pi[existing_cls]

  init_phi=if(!cl_phi){fit$phi[idx]}else{matrix(fit$phi[idx,],nrow=sum(idx))}

  offsets=log2(size_factors_pred)


  # nb log(f_k(y_i))
  l<-matrix(0,nrow=k,ncol=n)
  if(covars){
    covar_coefs = matrix(init_coefs[,-(1:k)],ncol=p)
    cov_eff = X %*% t(covar_coefs)         # n x g matrix of covariate effects
  } else {cov_eff=matrix(0,nrow=n,ncol=g)}

  for(i in 1:n){
    for(c in 1:k){
      if(cl_phi){
        l[c,i]<-sum(dnbinom(y_pred[,i],size=1/init_phi[,c],mu=2^(init_coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))    # posterior log like, include size_factor of subj
      } else if(!cl_phi){
        l[c,i]<-sum(dnbinom(y_pred[,i],size=1/init_phi,mu=2^(init_coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))
      }
    }    # subtract out 0.1 that was added earlier
  }

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
