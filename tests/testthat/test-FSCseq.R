test_that("simulateData works", {
  n=100
  k=3
  g=1000

  expect_equal(dim(simulateData(n,k,g,
               pi=rep(1/k,k),coefs=matrix(8,nrow=g,ncol=k),disp="gene",phi=rep(0.2,g),
               size_factors=NULL,batch_effects=0,batch=rep(1,n))$y),c(g,n))

  # wrong dimensions for coefs
  expect_error(simulateData(n,k,g,coefs=matrix(8,nrow=n,ncol=k),phi=rep(0.2,g),
                            size_factors=NULL,batch_effects=0,batch=rep(1,n)))

  # wrong dimension for phi with specified disp
  expect_error(simulateData(n,k,g,coefs=matrix(8,nrow=g,ncol=k),disp="cluster",phi=rep(0.2,g),
                            size_factors=NULL,batch_effects=0,batch=rep(1,n)))

  # wrong value for disp
  expect_error(simulateData(n,k,g,coefs=matrix(8,nrow=g,ncol=k),disp="?",phi=rep(0.2,g),
                            size_factors=NULL,batch_effects=0,batch=rep(1,n)))

  # two batches specified with one batch effect
  expect_error(simulateData(n,k,g,coefs=matrix(8,nrow=g,ncol=k),disp="gene",phi=rep(0.2,g),
                            size_factors=NULL,batch_effects=0,batch=c(rep(1,n/2),rep(2,n/2))))

  # can't have negative or 0 size factors
  expect_error(simulateData(n,k,g,coefs=matrix(8,nrow=g,ncol=k),disp="gene",phi=rep(0.2,g),
                            size_factors=rep(-0.1,n),batch_effects=0,batch=rep(1,n)))
  })

test_that("processData works",{
  n=100
  k=3
  g=1000

  sim.dat=simulateData(n,k,g,
      pi=rep(1/k,k),coefs=matrix(8,nrow=g,ncol=k),disp="gene",phi=rep(0.2,g),
      size_factors=NULL,batch_effects=0,batch=rep(1,n))

  # count data missing
  expect_error(processData(y=NA,geoMeans=NULL,estimateSFtype="ratio",
                           med_filt=TRUE,MAD_filt=TRUE,
                           med_thresh=100,MAD_quant_thresh=50))

  # undefined estimateSFtype
  expect_error(processData(y=sim.dat$y,geoMeans=NULL,estimateSFtype="?",
                           med_filt=TRUE,MAD_filt=TRUE,
                           med_thresh=100,MAD_quant_thresh=50))

  # MAD quantile outside of range 0 to 100
  expect_error(processData(y=sim.dat$y,geoMeans=NULL,estimateSFtype="ratio",
                           med_filt=TRUE,MAD_filt=TRUE,
                           med_thresh=100,MAD_quant_thresh=-2))

})

test_that("FSCseq works",{
  n=100
  k=3
  g=1000

  sim.dat=simulateData(n,k,g,
                       pi=rep(1/k,k),coefs=matrix(8,nrow=g,ncol=k),disp="gene",phi=rep(0.2,g),
                       size_factors=NULL,batch_effects=0,batch=rep(1,n))
  processed.dat=processData(y=sim.dat$y,geoMeans=NULL,estimateSFtype="ratio",
                           med_filt=TRUE,MAD_filt=TRUE,
                           med_thresh=100,MAD_quant_thresh=50)

  y=sim.dat$y
  norm_y=processed.dat$norm_y
  size_factors=processed.dat$size_factors

  expect_error(FSCseq(X=rep(1,99), y, k,
                      lambda=0,alpha=0,
                      size_factors=size_factors,
                      norm_y=norm_y,
                      offsets=rep(0,ncol(y)),
                      true_clusters=NULL, true_disc=NULL,
                      init_parms=FALSE,
                      init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                      init_phi=matrix(0,nrow=nrow(y),ncol=k),
                      init_cls=NULL,
                      maxit_EM=100,maxit_IRLS=50,EM_tol=1E-8,IRLS_tol=1E-6,
                      disp="gene",
                      method="EM",init_temp=nrow(y),trace=F))
  expect_error(FSCseq(X=NULL, y, k,
                      lambda=-1,alpha=0,
                      size_factors=size_factors,
                      norm_y=norm_y,
                      offsets=rep(0,ncol(y)),
                      true_clusters=NULL, true_disc=NULL,
                      init_parms=FALSE,
                      init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                      init_phi=matrix(0,nrow=nrow(y),ncol=k),
                      init_cls=NULL,
                      n_rinits=1,
                      maxit_EM=100,maxit_IRLS=50,EM_tol=1E-8,IRLS_tol=1E-6,
                      disp="gene",
                      method="EM",init_temp=nrow(y),trace=F))
  expect_error(FSCseq(X=NULL, y, k,
                      lambda=0,alpha=1,
                      size_factors=size_factors,
                      norm_y=norm_y,
                      offsets=rep(0,ncol(y)),
                      true_clusters=NULL, true_disc=NULL,
                      init_parms=FALSE,
                      init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                      init_phi=matrix(0,nrow=nrow(y),ncol=k),
                      init_cls=NULL,
                      n_rinits=1,
                      maxit_EM=100,maxit_IRLS=50,EM_tol=1E-8,IRLS_tol=1E-6,
                      disp="gene",
                      method="?",init_temp=nrow(y),trace=F))
  expect_error(FSCseq(X=NULL, y, k,
                      lambda=0,alpha=1,
                      size_factors=size_factors,
                      norm_y=norm_y,
                      offsets=rep(0,ncol(y)),
                      true_clusters=NULL, true_disc=NULL,
                      init_parms=FALSE,
                      init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                      init_phi=matrix(0,nrow=nrow(y),ncol=k),
                      init_cls=NULL,
                      n_rinits=1,
                      maxit_EM=100,maxit_IRLS=50,EM_tol=1E-8,IRLS_tol=1E-6,
                      disp="?",
                      method="CEM",init_temp=nrow(y),trace=F))
})

test_that("FSCseq_predict works",{
  n=20
  k=3
  g=100

  sim.dat=simulateData(n,k,g,
                       pi=rep(1/k,k),coefs=matrix(8,nrow=g,ncol=k),disp="gene",phi=rep(0.2,g),
                       size_factors=NULL,batch_effects=0,batch=rep(1,n))
  processed.dat=processData(y=sim.dat$y,geoMeans=NULL,estimateSFtype="ratio",
                            med_filt=TRUE,MAD_filt=TRUE,
                            med_thresh=100,MAD_quant_thresh=50)

  y=sim.dat$y
  norm_y=processed.dat$norm_y
  size_factors=processed.dat$size_factors

  fit=FSCseq(X=NULL, y, k,
             lambda=0,alpha=0,
             size_factors=size_factors,
             norm_y=norm_y,
             offsets=rep(0,ncol(y)),
             true_clusters=NULL, true_disc=NULL,
             init_parms=FALSE,
             init_coefs=matrix(0,nrow=nrow(y),ncol=k),
             init_phi=matrix(0,nrow=nrow(y),ncol=k),
             init_cls=NULL,n_rinits=1,
             maxit_EM=100,maxit_IRLS=50,EM_tol=1E-8,IRLS_tol=1E-6,
             disp="gene",
             method="EM",init_temp=nrow(y),trace=F)


  expect_error(FSCseq_predict(X=NULL,fit,y[,1:5],size_factors[1:6]))

})
