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
#' theta.ml() with small ridge penalty to stabilize estimates, and using more efficient
#' digamma and trigamma functions.
#'
#' @param y cts matrix
#' @param mu mean parameters
#' @param n number of samples
#' @param weights number of weights
#' @param t0 initial estimate of dispersion scale parameter
#' @param limit max number of iterations
#' @param eps tolerance
#' @param trace boolean whether to trace progress
#'
#' @return theta
#'
theta.ml2=function (y, mu, n = sum(weights), weights, t0=0, limit = 10, eps = .Machine$double.eps^0.25,
                    trace = FALSE)
{
  lambda=1E-25
  # score <- function(n, th, mu, y, w,lambda=1E-25) {sum(w * (digamma(th +
  #                                                                     y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +
  #                                                                                                                        th)/(mu + th))) + th*lambda}
  # info <- function(n, th, mu, y, w,lambda=1E-25) {sum(w * (-trigamma(th +
  #                                                                      y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +
  #                                                                                                                           th)^2)) + lambda}
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
  if(t0==0){t0 <- n/sum(weights * (y/mu - 1)^2)}
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
SCAD_soft_thresholding=function(diff_beta,lambda,alpha){
  a=3.7
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
    # CEM AE-step update on weights
    logdenom = apply((1/Tau)*(log(pi)+l),2,logsumexpc)
    for(c in 1:k){
      wts[c,]<-exp((1/Tau)*(log(pi[c])+l[c,])-logdenom)
    }
    if(Tau>1){
      Tau = 0.9*Tau
    } else{
     Tau=1       # after Tau hits 1 --> EM (off)
     CEM=F
    }
  }

  if(CEM){
    # C step
    draw_wts=wts                 # initialize
    for(i in 1:n){
      set.seed(i) # for reproducibility. stabilizes param ests in beginning iterations
      draw_wts[,i] = rmultinom(1,1,wts[,i])
    }
    seed_mult=1
    # have one sample per cluster --> stability in CEM
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
      # if no sample in one cluster still, have last sample's PP=1/k for all cls
      if(seed_mult>250){
        draw_wts[,n]=rep(1/k,k)
        break
      }
    }
    wts=draw_wts
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
phi.ml=function(y,mu,dfr,weights,t0,limit,trace){
  theta_g = NULL
  #try(theta_g <- theta_ml_g(y=y, mu=mu, wts=weights, t0=t0, limit=limit, trace=(trace)^2),silent=TRUE)
  try(theta_g <- theta.ml2(y=y, mu=mu, weights=weights, limit=limit, t0=t0, trace=trace),silent=TRUE)

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
               t0=0,
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
               t0=1/phi[j,1],
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
M_step_par = function(j){
  if(Tau<=1 & a>6){
    if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){
      res=par_X[[j]]
      return(res)
    }}
  res = M_step(X=XX, y_j=as.numeric(rep(y[j,],k)), p=p, j=j, a=a, k=k,
               all_wts=wts, keep=c(t(keep)), offset=rep(offsets,k),
               theta=theta_list[[j]],coefs_j=coefs[j,],phi_j=phi[j,],
               cl_phi=cl_phi,est_covar=est_covar[j],est_beta=est_beta[j],
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
                      t0=1/phi[j,1],
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
M_step_par2 = function(j, XX, y, p, a, k,
                       wts, keep, offsets,
                       theta_list, coefs, phi,
                       cl_phi, est_phi, est_covar, est_beta,
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
               cl_phi=cl_phi,est_covar=est_covar[j],est_beta=est_beta[j],
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
                      t0=1/phi[j,1],
                      limit=25,trace=F)
    res$phi_j = rep(phi_g_temp,k)
  } else{res$phi_j = phi[j,]}
  return(res)
}


#' Wrapper for main FSCseq function
#'
#' Run main CEM/EM clustering and feature selection algorithm.
#'
#' @param ncores integer, number of cores to utilize in parallel computing (default 1)
#' @param X (optional) design matrix of dimension n by p
#' @param y count matrix of dimension g by n
#' @param k integer, number of clusters
#' @param lambda numeric penalty parameter, lambda >= 0
#' @param alpha numeric penalty parameter, 0 <= alpha < 1
#' @param size_factors numeric vector of length n, factors to correct for subject-specific variation of sequencing depth
#' @param norm_y count matrix of dimension g by n, normalized for differences in sequencing depth
#' @param true_clusters (optional) integer vector of true groups, if available, for diagnostic tracking
#' @param true_disc (optional) logical vector of true discriminatory genes, if available, for diagnostic tracking
#' @param init_parms logical, TRUE: custom parameter initializations, FALSE (default): start from scratch
#' @param init_coefs matrix of dimension g by k, only if init_parms = TRUE
#' @param init_phi vector of dimension g (gene-specific dispersions), only if init_parms = TRUE
#' @param init_cls (optional) vector of length n, initial clustering.
#' @param init_wts (optional) matrix of dim k by n to denote initial clustering (allowing partial membership). If both init_cls and init_wts specified, init_wts will be ignored and init_cls used as initial clusters
#' @param n_rinits integer, number of additional random initializations to be searched (default 1)
#' @param maxit_inits integer, maximum number of iterations for each initialization search (default 100)
#' @param maxit_EM integer, maximum number of iterations for full CEM/EM run (default 100)
#' @param maxit_IRLS integer, maximum number of iterations for IRLS algorithm, in M step (default 50)
#' @param maxit_CDA integer, maximum number of iterations for CDA loop (default 50)
#' @param EM_tol numeric, tolerance of convergence for EM/CEM, default is 1E-6
#' @param IRLS_tol numeric, tolerance of convergence for IRLS, default is 1E-4
#' @param CDA_tol numeric, tolerance of convergence for IRLS, default is 1E-4
#' @param method string, either "EM" or "CEM" (default)
#' @param init_temp numeric, default for CEM: init_temp = nrow(y), i.e. number of genes. temp=1 for EM
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
#' @param trace.file (optional) string, file into which interim diagnostics will be printed
#' @param mb_size minibatch size: # of genes to include per M step iteration
#' @param PP_filt numeric between (0,1), threshold on PP for sample/cl to be included in M step estimation. Default is 1e-3
#'
#' @return list containing outputs from EM_run() function
#'
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
#'
#' @export
FSCseq<-function(ncores=1,X=NULL, y, k,
                 lambda=0,alpha=0,
                 size_factors,norm_y,
                 true_clusters=NULL, true_disc=NULL,
                 init_parms=FALSE,init_coefs=NULL,init_phi=NULL,init_cls=NULL,init_wts=NULL,
                 n_rinits=if(method=="EM"){20}else if(method=="CEM"){1},         # fewer searches for CEM to minimize computational cost
                 maxit_inits=if(method=="EM"){15}else if(method=="CEM"){100},  # for CEM, tolerates end temp (Tau) of 2 at end of initialization
                 maxit_EM=100,maxit_IRLS=50,maxit_CDA=50,EM_tol=1E-6,IRLS_tol=1E-4,CDA_tol=1E-4,
                 disp="gene",
                 method="CEM",init_temp=nrow(y),
                 trace=F,trace.file=NULL,
                 mb_size=NULL,PP_filt=1e-3){
  # disp = "gene". (disp="cluster" turned off)
  disp="gene"
  # set PP_filt = NULL or 0 < PP_filt < 1e-50: turn off PP_filt in M step (not recommended due to computation time)
  # y: raw counts
  # k: #clusters
  # size_factors: SF's derived from DESeq2
  # norm_y: counts normalized for sequencing depth by DESeq2
  # true_clusters: if applicable. For diagnostics tracking of ARI
  # true_disc: if applicable. For diagnostics tracking of disc/nondisc genes
  # init_parms: TRUE if initial coefficient estimates/dispersion estimates are input
  # init_coefs & init_phi: Initial estimates, if applicable
  # disp = c(gene, cluster), depending on whether dispersions are gene-specific (default) or cluster-specific (turned off)
  # init_cls: Initial clustering
  # n_rinits: Number of initial clusterings searched with maxit=15. More initializations = more chance to attain global max

  start_FSC = Sys.time()

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

  }


  # CEM is set to F for final run to convergence (after initialization).
  results=EM_run(ncores,X,y,k,lambda,alpha,size_factors,norm_y,true_clusters,true_disc,
                 init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi,disp=disp,
                 init_cls=init_cls,init_wts=init_wts,CEM=F,init_Tau=1,
                 maxit_EM=maxit_EM,maxit_IRLS=maxit_IRLS,maxit_CDA=maxit_CDA,
                 EM_tol=EM_tol,IRLS_tol=IRLS_tol,CDA_tol=CDA_tol,
                 trace=trace,mb_size=mb_size,PP_filt=PP_filt)

  end_FSC = Sys.time()
  results$total_time_elap = as.numeric(end_FSC-start_FSC,units="secs")

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
#' @param ncores integer, number of cores to utilize in parallel computing (default 1)
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
#' @param init_cls vector of length n, initial clustering.
#' @param init_wts matrix of dim k x n: denotes cluster memberships, but can have partial membership. init_wts or init_cls must be initialized
#' @param CEM logical, TRUE for CEM (default), FALSE for EM
#' @param init_Tau numeric, initial temperature for CEM. Default is g for CEM (set to 1 for EM)
#' @param maxit_EM integer, maximum number of iterations for full CEM/EM run (default 100)
#' @param maxit_IRLS integer, maximum number of iterations for IRLS loop, in M step (default 50)
#' @param maxit_CDA integer, maximum number of iterations for CDA loop (default is 50)
#' @param EM_tol numeric, tolerance of convergence for EM/CEM, default is 1E-6
#' @param IRLS_tol numeric, tolerance of convergence for IRLS, default is 1E-4
#' @param CDA_tol numeric, tolerance of convergence for CDA, default is 1E-4
#' @param disp string, either "gene" (default) or "cluster"
#' @param trace logical, TRUE: output diagnostic messages, FALSE (default): don't output
#' @param mb_size minibatch size: # of genes to include per M step iteration
#' @param PP_filt numeric between (0,1), threshold on PP for sample/cl to be included in M step estimation. Default is 1e-3
#'
#' @return FSCseq object with clustering results, posterior probabilities of cluster membership, and cluster-discriminatory status of each gene
#'
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
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
                   CEM=T,init_Tau=nrow(y),
                   maxit_EM=100, maxit_IRLS = 50,maxit_CDA=50,EM_tol = 1E-6,IRLS_tol = 1E-4, CDA_tol=1E-4, disp,trace=F,
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
  #coefs_list <- list()
  MSE_coefs = rep(0,maxit_EM)
  LFCs = matrix(0,nrow=g,ncol=maxit_EM)

  offsets=log2(size_factors)

  Q<-rep(0,times=maxit_EM)

  # if both init_wts and init_cls are specified, will default to init_wts for current clustering
  if(!is.null(init_wts)){
    wts=init_wts

    cls=apply(wts,2,which.max)
    current_clusters=cls
  } else if(!is.null(init_cls)){
    cls=init_cls
    current_clusters<-cls

    # Initialize weights
    wts<-matrix(0,nrow=k,ncol=n)
    for(c in 1:k){
      wts[c,]=(cls==c)^2
    }
  }

  # For use in CEM in E step #
  Tau = init_Tau

  phi_g = rep(0,times=g)
  DNC=0

  disc_ids_list = list()
  disc_ids=rep(T,g)

  diff_phi=matrix(0,nrow=maxit_EM,ncol=g)

  est_phi=rep(1,g)                          # 1 for true, 0 for false
  est_covar = if(covars){rep(1,g)}else{rep(0,g)}
  est_beta = rep(1,g)

  # if PP_filt is not set to NULL --> keep only obs/cl in M step > PP_filt threshold
  if(!is.null(PP_filt)){
    keep = (wts > PP_filt)^2
  }else{keep = matrix(1,nrow=nrow(wts),ncol=ncol(wts))}

  # if any cluster has 0 samples to estimate with, set keep to all samples --> weight them w/ small 1e-50 weights
  if(k>1){if(any(rowSums(keep)==0)){
    keep[rowSums(keep)==0,]=1
  }}

  ## Manual turn-off of estimating phi/covariates & dummy trackers
  # if(init_parms){
  #   est_phi=rep(0,g)
  #   est_covar = rep(0,g)
  # }
  # all_temp_list = list()
  # all_theta_list = list()
  # lower_K=FALSE         # tracks when number of a mixture components --> 0

  cl_agreement = rep(NA,maxit_EM)
  par_X=rep(list(list()),g)         # store M_step results

  nits_IRLS=rep(0,g)
  nits_CDA=rep(0,g)

  ########### M / E STEPS #########
  for(a in 1:maxit_EM){

    EMstart= Sys.time()

    prev_clusters = apply(wts,2,which.max)     # set previous clusters (or initial if a=1)
    if(a==1){         # Initializations for 1st EM iteration
      start=Sys.time()
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
          # parallelized Poisson glm + phi initialization
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
          # unparallelized loop Poisson glm + phi init
          for(j in 1:g){
            #j,y_j,XX,k,covars,p,offsets,wts,keep
            init_fit=glm.init(j,y[j,],XX,k,offsets,wts,keep)
            coefs[j,]=init_fit$coefs_j
            phi[j,]=rep(init_fit$phi_g,k)
            phi_g[j]=init_fit$phi_g
          }
        }
      }
      #pi=rowMeans(wts)
      for(j in 1:g){
        theta<-matrix(rep(0,times=k^2),nrow=k)
        for(c in 1:k){
          for(cc in 1:k){
            # run just once in R, for initialization
            if(cc==c){
              theta[cc,c]=0
            }else if(cc>c){
              theta[cc,c]<-SCAD_soft_thresholding(coefs[j,cc]-coefs[j,c],lambda,alpha)
            }else{theta[cc,c]=-theta[c,cc]}
          }
        }
        ### theta[lower.tri(theta)] = -t(theta)[lower.tri(theta)] # if upper tri matrix created --> fill in lower-tri mat

        theta_list[[j]] <- theta
        #disc_ids[j]=any(theta_list[[j]]!=0)
      }

      end=Sys.time()
      if(trace){
        cat(paste("Parameter Estimates Initialization Time Elapsed:",as.numeric(end-start,units="secs"),"seconds.\n"))
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

    # # update on pi_hat, and UB & LB on pi
    # for(c in 1:k){
    #   pi[c]=mean(wts[c,])
    #   if(pi[c]<1E-6){
    #     if(trace){warning(paste("cluster proportion", c, "close to 0"))}
    #     pi[c]=1E-6
    #   } # lowerbound for pi
    #   if(pi[c]>(1-1E-6)){
    #     if(trace){warning(paste("cluster proportion", c, "close to 1"))}
    #     pi[c]=(1-1E-6)
    #   } # upperbound for pi
    # }

    pi=rowMeans(wts)
    pi[pi<1e-6]=1e-6; pi[pi>(1-1e-6)]=1-1e-6  # for stability in log(pi) in Q function

    # M step
    Mstart=Sys.time()

    if(!is.null(mb_size)){
      mb_genes = sample(1:g,mb_size,replace=F)
    } else{ mb_genes=1:g}

    if(ncores>1){
      # M_step parallelized across ncores
      if(Sys.info()[['sysname']]=="Windows"){
        ## use parLapply for Windows ##
        clust = makeCluster(ncores)
        clusterEvalQ(cl=clust,library(FSCseq))
        clusterExport(cl=clust,varlist=c("XX","y","p","a","k","wts","keep","offsets",
                                         "theta_list","coefs","phi","cl_phi","est_phi","est_covar","est_beta",
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
                      cl_phi, est_phi, est_covar, est_beta,
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
                             cl_phi=cl_phi,est_covar=est_covar[j],est_beta=est_beta[j],
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
                            t0=1/phi[j,1],
                            limit=25,
                            trace=F)
          par_X[[j]]$phi_j = rep(phi_g_temp,k)
        } else{par_X[[j]]$phi_j = phi[j,]}
      }
    }
    Mend=Sys.time()
    if(trace){cat(paste("M Step Time Elapsed:",as.numeric(Mend-Mstart,units="secs"),"seconds.\n"))}

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
      nits_IRLS[j]=par_X[[j]]$nits_IRLS
      nits_CDA[j]=par_X[[j]]$nits_CDA

      # calculate LFCs (should just be 0 if k=1)
      LFCs[j,a] = max(coefs[j,1:k])-min(coefs[j,1:k])

      if(a>6){
        # set relative change in phi across 5 iterations
        if(cl_phi==1){
          diff_phi[a,j]=mean(abs(phi[j,]-phi_list[[a-5]][j,])/(phi_list[[a-5]][j,] + 0.001))
        } else if(cl_phi==0){
          diff_phi[a,j]=abs(phi_g[j]-phi_list[[a-5]][j,1])/(phi_list[[a-5]][j,1] + 0.001)
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

    ## dummy trackers: turned off
    # all_temp_list[[a]] = temp_list
    # all_theta_list[[a]] = theta_list
    # coefs_list[[a]] = coefs

    # Marker of all nondisc genes (T for disc, F for nondisc)
    if(trace){
      cat(paste("Disc genes:",sum(disc_ids),"of",g,"genes.\n"))
    }

    # save current objects in lists (to track)
    disc_ids_list[[a]] = disc_ids
    phi_list[[a]] <- phi

    if(a>1){
      #MSE_coefs[a]=mean(abs((coefs_list[[a]]-coefs_list[[a-1]])))
      MSE_coefs[a] = mean(abs(coefs-prev_coefs))
      if(trace){cat(paste("MSE_coef:",MSE_coefs[a],"\n"))}
      cl_agreement[a] = mean(current_clusters==prev_clusters)
    }

    prev_coefs = coefs    # store prev coefs for comparison in MSE_coefs next iteration

    # nb log(f_k(y_i))
    l<-matrix(rep(0,times=k*n),nrow=k,ncol=n)
    # if(covars){
    #   covar_coefs = matrix(coefs[,(k+1):(k+p)],ncol=p)
    #   cov_eff = X %*% t(covar_coefs)         # n x g matrix of covariate effects
    # } else {cov_eff=matrix(0,nrow=n,ncol=g)}
    # for(i in 1:n){
    #   for(c in 1:k){
    #     # All genes. nondisc genes should be cancelled out in E step calculation
    #     # if(cl_phi==1){
    #     #   l[c,i]<-sum(dnbinom(y[,i],size=1/phi[,c],mu=2^(coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))    # posterior log like, include size_factor of subj
    #     # } else if(cl_phi==0){
    #     #   l[c,i]<-sum(dnbinom(y[,i],size=1/phi_g,mu=2^(coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))
    #     # }
    #     ## phi[j,] = rep(phi_g[j],k) ##
    #     l[c,i]<-sum(dnbinom(y[,i],size=1/phi[,c],mu=2^(coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))
    #   }
    # }

    if(covars){
      covar_coefs = matrix(coefs[,(k+1):(k+p)],ncol=p)
      cov_eff = covar_coefs %*% t(X)         # g x n matrix of covariate effects
    } else {cov_eff=matrix(0,nrow=g,ncol=n)}
    offset_eff = matrix(offsets,nrow=g,ncol=n,byrow=T)
    for(c in 1:k){
      # l[c,] = colSums(
      #   dnbinom(y, size=1/matrix(phi[,c],nrow=g,ncol=n,byrow=F),
      #           mu=2^(matrix(coefs[,c],nrow=g,ncol=n,byrow=F) + cov_eff + offset_eff),log=TRUE)
      #   )
      l[c,] = colSums(
        dnbinom(y, size=1/phi[,c],
                mu=2^(coefs[,c] + cov_eff + offset_eff),log=TRUE)
      )
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
        }
      }
    }
    if(a>10){if(all(cl_agreement[(a-9):a]==1) & Tau==1){
      # if cls are not changing for a while
      finalwts<-wts
      break
    }}
    if(a==maxit_EM){
      finalwts<-wts
      if(trace){warning("Reached max iterations.")}
      DNC=1
      break
    }

    # E step
    if(trace){cat(paste("Tau:",Tau,"\n"))}
    start_E = Sys.time()
    Estep_fit=E_step(wts,l,pi,CEM,Tau,PP_filt)
    wts=Estep_fit$wts; keep=Estep_fit$keep; Tau=Estep_fit$Tau; CEM=Estep_fit$CEM
    end_E = Sys.time()
    if(trace){cat(paste("E Step Time Elapsed:",as.numeric(end_E-start_E,units="secs"),"seconds\n"))}

    # if any cluster has 0 samples to estimate with, set keep to all samples --> weight them w/ small 1e-50 weights
    if(k>1){if(any(rowSums(keep)==0)){
      keep[rowSums(keep)==0,]=1
    }}

    # Diagnostics Tracking
    current_clusters = apply(wts,2,which.max)

    #print(current_clusters)
    if(trace){
      cat(paste("EM iter",a,"% of samples whose cls unchanged (from previous):",mean(current_clusters==prev_clusters),"\n"))
      if(!is.null(true_clusters)){cat(paste("ARI =",adjustedRandIndex(true_clusters,current_clusters),"\n"))}
      cat(paste("Cluster proportions:",pi,"\n"))

      # if true_disc gene list is input
      if(!is.null(true_disc)){
        TPR=sum(true_disc & disc_ids)/sum(true_disc)
        FPR=sum(!true_disc & disc_ids)/sum(!true_disc)
        cat(paste("TPR:",TPR,", FPR:",FPR,"\n"))

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
      EMend = Sys.time()
      cat(paste("EM iter",a,"time elapsed:",as.numeric(EMend-EMstart,units="secs"),"seconds.\n"))
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
    # m<-rep(0,times=g)
    # for(j in 1:g){
    #   m[j]=length(unique(theta_list[[j]][1,]))
    #   if(m[j]==1){nondiscriminatory[j]=TRUE}
    # }
    m = sapply(theta_list, function(x){return(length(unique(x[1,])))})
    nondiscriminatory[m==1] = TRUE
  }

  num_est_coefs = sum(m)
  num_est_params =
    if(cl_phi==1){ # cl_phi turned off
      sum(m)+(k-1)+p*g + k*g                # p*g for covariates, sum(m) for coefs, k*g for cluster phis, k-1 for mixture proportions
    } else{ sum(m)+(k-1)+g+p*g }            # sum(m) for coef for each discriminatory clusters (cl_phi=1). sum(m) >= g

  log_L<-sum(apply(log(pi) + l, 2, logsumexpc))

  eff_n=n*g
  BIC_penalty_term = log(eff_n)*num_est_params
  BIC = -2*log_L + BIC_penalty_term

  if(trace){

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
  time_elap <- as.numeric(end_time-start_time,units="secs")

  if(cl_phi==0){
    phi = phi_g
  }

  # omit some output for memory reasons
  result<-list(k=k,
               pi=pi,
               coefs=coefs,
               Q=Q[1:a],
               BIC=BIC,
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
               DNC=DNC,LFCs=LFCs,n_its=a#,lower_K=lower_K
  )
  return(result)
}

#' Prediction via FSCseq
#'
#' Performs prediction of cluster membership via trained model
#'
#' @param X (optional) design matrix of dimension n by p
#' @param fit FSCseq results object
#' @param cts_train Training count data matrix of dimension g by n_train
#' @param cts_pred Prediction/test count data matrix of dimension g by n_pred
#' @param SF_train Vector of length n_train, size factors for training subjects. Can be accessed from the fit.
#' @param SF_pred Vector of length n_pred, size factors for prediction subjects
#' @param maxit Maximum number of iterations to run if prediction batches are estimated, default is 100.
#' @param eps Tolerance for relative change in Q function convergence criterion if prediction batches are estimated, default is 1e-4.
#'
#' @return list containing outputs
#' final_clusters: vector of length n of resulting clusters,
#' wts: k by n matrix of E step weights
#'
#' @author David K. Lim, \email{deelim@live.unc.edu}
#' @references \url{https://github.com/DavidKLim/FSCseq}
#'
#' @export
FSCseq_predict <- function(X=NULL, fit, cts_train=NULL, cts_pred,
                           SF_train=NULL, SF_pred, maxit=100, eps=1e-4){   # NEED TO CHANGE ORDER OF INPUT IN FSCseq_workflow
  # fit: Output of EM
  # cts_pred: Data to perform prediction on
  # SF_pred: SF's of new data
  # offsets: Additional offsets per sample can be incorporated

  covars = !is.null(X)
  idx=fit$discriminatory
  if(sum(idx)==0){warning('No disc genes. Using all genes'); idx=rep(TRUE,nrow(cts_pred))}

  if(covars){ cts = cbind(cts_train[idx,], cts_pred[idx,]) } else{ cts = cts_pred[idx,] }

  unique_cls_ids = unique(fit$clusters)[order(unique(fit$clusters))]   # unique clusters, in increasing order
  k=length(unique_cls_ids)

  train_clusters = fit$clusters
  for(i in 1:length(train_clusters)){
    train_clusters[i] = which(unique_cls_ids==fit$clusters[i])     # renumber into 1,2,3,...,K. if already in this order, won't change anything
  }
  n_train = length(train_clusters)
  n_pred = ncol(cts_pred)
  n = ncol(cts)
  g = sum(idx)
  print(paste("n:",n,", g:",g))
  pi_all = fit$pi  # including penalized out clusters here. pi is changed later
  coefs = fit$coefs

  # B = length(unique(batch_train))
  if(is.null(X)){
    p_train = ncol(fit$coefs) - length(pi_all)   # number of non-log2 cluster means estimated in training
  }else{p_train=ncol(X)}

  # coefs2 = matrix(0,nrow=nrow(fit$coefs),ncol=ncol(fit$coefs)+p_train)
  # coefs2[,1:ncol(fit$coefs)] = fit$coefs
  # fit$coefs = coefs2

  # apply(table(fit$clusters,cls),2,function(x){which(x==max(x))})      # returns
  # cls_match_ids

  if(length(SF_pred) != ncol(cts_pred)){stop("length of prediction size factors must be the same as the number of columns in prediction data")}

  cl_X = matrix(0,nrow=k*n,ncol=k)
  ident_k = diag(k)
  for(i in 1:k){ cl_X[((i-1)*n+1):(i*n),] = matrix(rep(ident_k[i,],n),ncol=k,nrow=n,byrow=T) }

  #cts_pred=matrix(cts_pred[idx,],ncol=ncol(cts_pred))      # subset to just disc genes found by FSCseq run
  #cts_train=matrix(cts_train[idx,],ncol=ncol(cts_train))

  # cl_phi=(!is.null(dim(fit$phi)))^2  # dimension of phi is null when gene-wise (vector)
  cl_phi=0

  if(covars){
    XX = do.call("rbind", replicate(k, X,simplify=FALSE))
    # print(dim(cl_X))
    # print(dim(XX))
    XX = cbind(cl_X,XX)
    p = ncol(XX) - k
    coefs = matrix(0,nrow=g,ncol=k+p)
    # coefs[,1:(k+p_train)] = fit$coefs[idx,1:(k+p_train)]
    betas=fit$coefs[idx,unique_cls_ids]; gammas = fit$coefs[idx,(length(pi_all)+1):ncol(fit$coefs)]
    coefs[,1:(k+p_train)] = cbind(betas,gammas)

    # initialize estimates of new batch effects as the mean of the other batch effects across batches for each gene
    # coefs[,(k+p_train+1):ncol(coefs)] <- if((length(pi_all)+1) == ncol(fit$coefs)){ fit$coefs[idx,(length(pi_all)+1):ncol(fit$coefs)]
    #                                                                                 } else{ rowMeans(fit$coefs[idx,(length(pi_all)+1):ncol(fit$coefs)]) }

    # coefs = cbind(coefs, if((length(pi_all)+1) == ncol(fit$coefs)){ fit$coefs[idx,(length(pi_all)+1):ncol(fit$coefs)]
    # } else{ rowMeans(fit$coefs[idx,(length(pi_all)+1):ncol(fit$coefs)]) })
    ## Initialize weights here (how do I initialize for new samples??):
    wts=matrix(0, nrow=k, ncol=n)

    init_cls_pred = sample(1:k,ncol(cts_pred),replace=T)

    # init_cls_pred = cls_pred ## TRUE INITIALIZATION
    #table(train_clusters,cls) # cls==2 --> train_clusters==1, cls==1 --> train_clusters==2
    # init_cls_pred[cls_pred==1]=2; init_cls_pred[cls_pred==2]=1    if using true initialization, need to worry about permutations..

    init_cls = c(train_clusters, init_cls_pred)

    for(c in 1:k){ wts[c,]=(init_cls==c)^2 }
    wts[,1:ncol(cts_train)] = fit$wts[unique_cls_ids,]

    #keep = wts > 0.001
    keep = matrix(1,nrow=k,ncol=n)  # just use all samples in prediction (not as many samples, so it shouldn't be too time consuming)

  } else{
    p=0
    coefs=matrix(fit$coefs[idx,unique_cls_ids],nrow=sum(idx))
    XX=cl_X
    wts=matrix(0, nrow=k, ncol=n)
    keep=matrix(1,nrow=k,ncol=n)
    pi =pi_all[unique_cls_ids]
  }

  # table(init_cls,c(cls,cls_pred))
  # adjustedRandIndex(init_cls,c(cls,cls_pred))
  # table(init_cls,apply(wts,2,which.max))
  # adjustedRandIndex(init_cls,apply(wts,2,which.max))

  # fit is the output object from the EM() function
  phi=if(cl_phi==0){matrix(fit$phi[idx], nrow=sum(idx), ncol=k)}else{matrix(fit$phi[idx,],nrow=sum(idx), ncol=k)}

  if(!is.null(SF_train)){SF_train = fit$size_factors}
  offsets=log2(c(SF_train,SF_pred))

  est_beta = rep(0,g)
  est_phi = rep(0,g)                          # 1 for true, 0 for false
  est_covar = if(covars){rep(1,g)}else{rep(0,g)}
  #lambda=fit$lambda; alpha=fit$alpha    # not estimating beta, so penalty parameter values don't matter
  lambda=0; alpha=0

  nits=ifelse(covars, maxit, 1)    # 1 iteration if no covariates (no prediction batch effects estimated), maxit otherwise
  # Tau=ifelse(covars,g^2,1); CEM=T
  Tau=1; CEM=F

  ###### INSERT LOOP HERE TO ESTIMATE NEW BATCH EFFECTS ######
  ###### IF(covars): M STEP WITH EST_BETA=0, EST_PHI=0, EST_COVAR=1 ######
  ######### NOTE: need to concatenate init_coefs (change this to "coefs) with new batch
  ###### ELSE: CALCULATE l[c,i] AND COMPUTE WTS, and don't iterate #####

  interm_cls = list()
  delta_cls = rep(NA,nits-1)
  Q = rep(NA,nits)
  for(a in 1:nits){
    if(covars){
      par_X=list(); temp_list=list(); nits_IRLS=rep(NA,g); nits_CDA=rep(NA,g)
      pi=rowMeans(wts)
      for(j in 1:g){
        theta<-matrix(rep(0,times=k^2),nrow=k)
        for(c in 1:k){
          for(cc in 1:k){
            # run just once in R, for initialization
            if(cc==c){
              theta[cc,c]=0
            }else if(cc>c){
              theta[cc,c]<-SCAD_soft_thresholding(coefs[j,cc]-coefs[j,c],lambda,alpha)
            }else{theta[cc,c]=-theta[c,cc]}
          }
        }
        # if(a==1){ ######## initialize batch effects #########
        #   # library(MASS)
        #   # X0 = (X-X[,1])[,-1]
        #   # glm_fit0 = glm.nb(cts[j,] ~ 1 + X0 + offset(offsets + coefs[j,init_cls]))
        #   # glm_fit = glm.nb(cts[j,] ~ 0 + X + offset(offsets + coefs[j,init_cls]))
        #   covariate_X = XX[,(k+1):(k+p)]
        #   # covariate_X0 = (covariate_X - covariate_X[,1])[,-1]
        #   # glm_fit0 = glm.nb(as.integer(rep(cts[j,],k)) ~ 1 + covariate_X0 + offset(rep(offsets,k) + rep(coefs[j,1:k],each=n)),weights=c(t(wts)))  # no intercept
        #   glm_fit = glm.nb(as.integer(rep(cts[j,],k)) ~ 0 + XX[,(k+1):(k+p)] + offset(rep(offsets,k) + rep(coefs[j,1:k],each=n)),weights=c(t(wts)))  # intercept
        #
        #   glm_fit = glm.nb(cts_train[idx,][j,] ~ 0 + as.factor(cls) + as.factor(batch_train) + offset(log2(SF_train)),
        #                    init.theta = 1/phi[j,1])
        #   log2(exp(glm_fit$coefficients)); 1/glm_fit$theta
        #
        #   glm_fit = glm.nb(cts[j,] ~ 0 + as.factor(init_cls) + as.factor(c(batch_train,batch_pred)) + offset(log2(c(SF_train, SF_pred))),
        #                    init.theta = 1/phi[j,1])
        #   log2(exp(glm_fit$coefficients)); 1/glm_fit$theta
        #
        #   coefs[j,(k+1):(k+p)] = glm_fit$coefficients
        # }


        ### Estimate all batch effects

        par_X[[j]] <- M_step(X=as.matrix(XX), y_j=as.numeric(rep(cts[j,],k)), p=p, j=j, a=a, k=k,
                               all_wts=wts, keep=c(t(keep)), offset=rep(offsets,k),
                               theta=theta,coefs_j=coefs[j,],phi_j=phi[j,],
                               cl_phi=cl_phi,est_covar=est_covar[j],est_beta=est_beta[j],
                               lambda=lambda,alpha=alpha,
                               IRLS_tol=1e-4,maxit_IRLS=50L,
                               CDA_tol=1e-4,maxit_CDA=50L)
        coefs[j,] <- par_X[[j]]$coefs_j


        if(est_phi[j]==1){
          ids = (c(t(keep))==1)
          mu = 2^(XX %*% par_X[[j]]$coefs_j + offsets)

          phi_g_temp=phi.ml(y=as.integer(rep(cts[j,],k))[ids],
                            mu=mu[ids],
                            dfr=sum(ids)-1,
                            weights=c(t(wts))[ids],
                            t0=1/phi[j,1],
                            limit=25,
                            trace=F)
          par_X[[j]]$phi_j = rep(phi_g_temp,k)
        } else{par_X[[j]]$phi_j = phi[j,]}

        #coefs[j,(ncol(coefs)-p_pred+1):ncol(coefs)] = par_X[[j]]$coefs_j[(ncol(coefs)-p_pred+1):ncol(coefs)]
        temp_list[[j]] <- if(p>0){cbind(par_X[[j]]$temp_beta, par_X[[j]]$temp_gamma)}else{par_X[[j]]$temp_beta}
        if(cl_phi==1){
          phi[j,] <- par_X[[j]]$phi_j
        } else if(cl_phi==0){
          phi[j,] <- (par_X[[j]]$phi_j)[1]
        }
        nits_IRLS[j]=par_X[[j]]$nits_IRLS
        nits_CDA[j]=par_X[[j]]$nits_CDA

      }
    }

      # nb log(f_k(y_i))
      l<-matrix(0,nrow=k,ncol=n)
      # if(covars){
      #   covar_coefs = matrix(coefs[,-(1:k)],ncol=p)
      #   cov_eff = X %*% t(covar_coefs)         # n x g matrix of covariate effects
      # } else {cov_eff=matrix(0,nrow=n,ncol=g)}
      # for(i in 1:n){
      #   for(c in 1:k){
      #     if(cl_phi){
      #       l[c,i]<-sum(dnbinom(cts[,i],size=1/phi[,c],mu=2^(coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))    # posterior log like, include size_factor of subj
      #     } else if(!cl_phi){
      #       l[c,i]<-sum(dnbinom(cts[,i],size=1/phi,mu=2^(coefs[,c] + cov_eff[i,] + offsets[i]),log=TRUE))
      #     }
      #   }    # subtract out 0.1 that was added earlier
      # }
      if(covars){
        covar_coefs = matrix(coefs[,(k+1):(k+p)],ncol=p)
        cov_eff = covar_coefs %*% t(X)         # g x n matrix of covariate effects
      } else {cov_eff=matrix(0,nrow=g,ncol=n)}
      offset_eff = matrix(offsets,nrow=g,ncol=n,byrow=T)
      # print("phi")
      # print(dim(phi))
      # print(head(phi))
      # print("coefs")
      # print(dim(coefs))
      # print(head(coefs))
      # print('cov_eff')
      # print(dim(cov_eff))
      # print(head(cov_eff))
      # print("offset_eff")
      # print(dim(offset_eff))
      # print(head(offset_eff))
      # print("l")
      # print(dim(l))
      # print("mu")
      # c=1
      # print(dim(2^(coefs[,c] + cov_eff + offset_eff)))
      # print(head(2^(coefs[,c] + cov_eff + offset_eff)))
      # print("cts")
      # print(head(cts))
      # print(dim(cts))
      for(c in 1:k){
        size = 1/phi[,c]
        mus = 2^(coefs[,c] + cov_eff + offset_eff)
        l[c,] = colSums(
          dnbinom(cts, size=size, mu=mus, log=TRUE)
        )
      }

      Estep_fit=E_step(wts,l,pi,CEM=CEM,Tau,1e-3)
      #keep=matrix(1,nrow=nrow(wts),ncol=ncol(wts))
      keep=Estep_fit$keep
      wts=Estep_fit$wts; Tau=Estep_fit$Tau

      if(covars){ wts[,1:ncol(cts_train)] = fit$wts[unique_cls_ids,] }

      interm_cls[[a]] = apply(wts,2,which.max)

      ### BREAK CRITERION ##
      #Cluster-label based
      #sum number of samples whose cluster labels change, and break if no change for 10 iterations
      # if(a>1){
      #   delta_cls[a-1] = sum((interm_cls[[a]] != interm_cls[[a-1]])^2)
      # }
      # if(a>10){
      #   if(all(delta_cls[(a-10):(a-1)] == 0)){
      #     break
      #   }
      # }

      #Q function based
      #if relative change of the Q function is within eps (default = 1e-4)
      Q[a]<- (log(pi)%*%rowSums(wts)) + sum(wts*l)
      if(a>5){
        if(abs((Q[a]-Q[a-5])/Q[a-5]) < eps){
          break
        }
      }
  }

  final_clusters = apply(wts,2,which.max)

  if(covars){
    pred_cls = final_clusters[(n_train+1):n]
    train_wts = wts[,1:n_train]
    pred_wts = wts[,(n_train+1):n]
  }else{
    pred_cls = final_clusters
    train_wts = fit$wts[unique_cls_ids,]
    pred_wts = wts
  }

  return(list(train_cls=train_clusters,
              pred_cls=pred_cls,
              train_wts=train_wts,
              pred_wts=pred_wts,
              coefs=coefs))
}
