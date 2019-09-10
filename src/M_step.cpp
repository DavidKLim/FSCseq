//
//  M.cpp
//
//
//  Created by KYUNG T LIM on 1/19/18.
// [[Rcpp::depends(RcppArmadillo)]]

/*
depends(RcppEigen)
depends(RcppNumerical)
#include <RcppNumerical.h>
 using namespace Numer;
 typedef Eigen::Map<Eigen::VectorXd> MapVec;
*/

#include <string>
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>


int sign(double x) {
    return (x > 0) - (x < 0);
}

double SCAD_soft_thresh(double theta, double lambda, double alpha){ /* ST of SCAD penalty */
  double STval;
	double a=3.7;

    /*if(fabs(theta)<=2*lambda*alpha){*/
	if(fabs(theta)<=(alpha/(1-alpha))+lambda*alpha){
  		if(fabs(theta)<=alpha/(1-alpha)){
  			STval = 0;
  		} else{
  			STval = sign(theta)*(fabs(theta)-alpha/(1-alpha));
  		}
    } else if(fabs(theta)>(alpha/(1-alpha))+lambda*alpha && fabs(theta)<=a*lambda*alpha){
		double omega = ((a-1)*theta)/(a-1-1/(lambda*(1-alpha)));
		if(fabs(omega) - (a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha))) <= 0){
			STval=0;
		} else{
			STval = sign(omega)*(fabs(omega)-(a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha))));
		}
    } else{
		STval = theta;
	}
    return(STval);
}

double lasso_soft_thresh(double alpha, double lambda){
	double STval;
	if(fabs(alpha)-lambda<0){
        STval = 0;
    } else {
        STval = sign(alpha) * (fabs(alpha)-lambda);
    }
	return(STval);
}

/*Rcpp::List score_info(int N, double ph, arma::vec mu, arma::vec y, arma::vec wts){
    double lambda = 1e-25;
    double score1 = 0, info1 = 0;
    double inv_ph = 1/ph + lambda;  // try adding lambda to stabilize (when phi very large)?
    ph=ph/(1+ph*1e-25);
    double scorei, infoi;

    int n = y.size();


	for(int i=0; i<n; i++){

        scorei = wts(i) * -pow(inv_ph,2) * (R::digamma(inv_ph+y(i)) - R::digamma(inv_ph) + log(inv_ph) - log(inv_ph+mu(i)) + 1 - (y(i)+inv_ph)/(inv_ph+mu(i)));
        infoi = wts(i) * pow(inv_ph,4) * (R::trigamma(inv_ph) + 2/(inv_ph+mu(i)) - R::trigamma(inv_ph+y(i)) - ph - (y(i)+inv_ph)/pow(inv_ph+mu(i),2)) -
          2 * pow(inv_ph,3) * scorei;

        score1 += scorei;
        info1 += infoi;
	}



    double score = score1 + 2*lambda*ph;
    double info = info1 + 2*lambda;

    return Rcpp::List::create(score,info);
}

double phi_ml_g(arma::vec y, arma::vec mu, arma::vec wts, int limit, double p0, int trace){
    double eps = 0.0001220703; // from R
    double N = accu(wts);
    int n = y.size();

    if(n==1){
        return(p0=0);
    }
    if(accu(mu)<1e-13){
        return(p0=0);
    }

    if(p0==0){
	    for(int i=0; i<n; i++){
		    p0 += wts(i)*pow(y(i)/mu(i)-1,2);
	    }
      p0 = p0/N;
    }

    int it=0;
    double del=1;

    if(trace==1){
        Rprintf("phi_ml: iter %d 'phi = %f' \n",it,p0);
    }
    while(it < limit && fabs(del) > eps){
        it += 1;
        p0 = fabs(p0);
        Rcpp::List scoreinfo = score_info(N,p0,mu,y,wts);
        double score=scoreinfo[0], info=scoreinfo[1];
        del = score/info;
        p0 += del;
        if(trace==1){
            Rprintf("score: %f\n",score);
            Rprintf("info: %f\n",info);
            Rprintf("phi_ml: iter %d 'phi = %f' \n",it,p0);
        }
    }

    if(p0 < 0){
        p0 = 0;
        if(trace==1){
            Rprintf("estimate truncated at zero \n");
        }
    }
    if(it == limit && trace==1){
        Rprintf("iteration limit reached \n");
    }
	if(trace==1){Rprintf("phi_ml() converged after %d iterations\n",it);}

    return(p0);
}*/

/*
class phiBFGS: public MFuncGrad
{
private:
  const MapVec Y;
  const MapVec Mu;
  const MapVec Wts;
public:
  phiBFGS(const MapVec y_, const MapVec mu_, const MapVec wts_) : Y(y_), Mu(mu_), Wts(wts_) {}

  double f_grad(Constvec& phi, Refvec grad)
  {
    double phi0;
    if(phi[0]<=0){
      phi0=1E-6;
    } else{ phi0=phi[0];}
    double N = Y.size();
    double sumgrad=0; double sumNLL=0; double invphi=1/phi0;
    double scorei; double NLLi;
    for(int i=0; i<N; i++){
      scorei = Wts[i] * std::pow(invphi,2)*(R::digamma(Y[i]+invphi) - R::digamma(invphi) + std::log(invphi) - std::log(invphi + Mu[i]) + 1 - (Y[i]+invphi)/(invphi+Mu[i]));
      NLLi = -Wts[i] * R::dnbinom_mu(Y[i],invphi,Mu[i],true);
      sumgrad = sumgrad + scorei;
      sumNLL = sumNLL + NLLi;
    }
    grad[0]=sumgrad;
    const double NLL=sumNLL;
    return NLL;
  }
}; */




/* */

//' Perform M step of EM algorithm
//'
//' @param X design matrix
//' @export
// [[Rcpp::export]]
Rcpp::List M_step(arma::mat X, arma::vec y_j, int p, int j, int a, int k,
                  arma::mat all_wts, arma::vec keep, arma::vec offset,
                  arma::mat theta, arma::vec coefs_j, arma::vec phi_j,
                  int cl_phi, int est_covar,
                  double lambda, double alpha, double IRLS_tol, int maxit_IRLS,
                  std::string optim_method){

  int N = y_j.size();
  int n = N/k;    // N must be integer: N=n*k

  // wts: kxn. vectorise() goes column-wise. transpose to go row-wise
  arma::vec vec_wts = arma::vectorise(trans(all_wts));

  double phi_g = phi_j(0);
  arma::vec beta = coefs_j.subvec(0,k-1);
  arma::vec b = coefs_j;     // for use in simplifying X*b later

  arma::mat temp_beta(maxit_IRLS,k); temp_beta.zeros();
  int continue_beta = 1;

  arma::vec gamma(k+p);
  arma::mat temp_gamma(maxit_IRLS,k+p);   // initialize to have (k+p) cols (guaranteed to be >0), then reshape
  int continue_gamma = 0;

  if(p>0){
    gamma.resize(p); gamma = coefs_j.subvec(k,k+p-1);
    temp_gamma.reshape(maxit_IRLS,p); temp_gamma.zeros();
    continue_gamma = 1;
  } else{
    gamma.reshape(0,0);
    temp_gamma.reshape(0,0);
  }

	// arma::vec n_k(k);
	// for(int c=0; c<k; c++){
	//   n_k(c) = sum(all_wts.row(c));
	// }

	/* full_resid : y_tilde - sum_j(x_ij*beta_j) */
	arma::vec y_tilde(N), eta(N), mu(N), vec_W(N), resid(N), full_resid(N);
	y_tilde.zeros(); eta.zeros(); mu.zeros(); vec_W.zeros(); resid.zeros(); full_resid.zeros();

	/*//WLS estimate of covariates
	arma::vec MLE_beta(p+k);*/

	Rcpp::Timer timer;
	int output_timer = 0;
	int TIME_ITER = 0;

	/* PP filtering ids, gene */
	arma::uvec ids = find(keep==1);
	arma::uvec ids_c = find(keep==5);     /* just dummy initialization */


    /* IRLS */
  for(int i=0; i<maxit_IRLS; i++){

		/* Calculate eta and mu */
		eta = X*b + offset;
		for(int ii=0; ii<N; ii++){
			mu(ii)=pow(2,eta(ii));
		}

    /* Initiate temp matrix to track IRLS */
    for(int c=0; c<k; c++){
			temp_beta(i,c) = beta(c);
    }
    for(int c=0; c<p; c++){
      temp_gamma(i,c) = gamma(c);    // won't do anything if p=0
    }

		/* Calculate y_tilde and W matrix */
		for(int ii=0; ii<N; ii++){
			y_tilde(ii) = ( eta(ii)-offset(ii) ) + (y_j(ii)-mu(ii))/(mu(ii)*log(2));      /* link function: log_2(mu) = eta */
			double w_ii;
			if(cl_phi==0){
			  w_ii = sqrt(vec_wts(ii)*pow(mu(ii)*log(2),2)/(mu(ii)+mu(ii)*mu(ii)*phi_g));
			}else{
			  arma::vec all_phi = X.cols(0,k) * phi_j;      // phi for each cluster
			  w_ii = sqrt(vec_wts(ii)*pow(mu(ii)*log(2),2)/(mu(ii)+pow(mu(ii),2)*all_phi(ii)));
			}
			vec_W(ii) = w_ii;
		}

		/* Update theta matrix */
		for(int c=0; c<k; c++){
			for(int cc=0; cc<k; cc++){
				theta(c,cc) = SCAD_soft_thresh(beta(c)-beta(cc),lambda,alpha);
				/*theta(c,cc) = lasso_soft_thresh(beta(c)-beta(cc),lambda*alpha);*/
			}
		}
		if(i==TIME_ITER){
		  timer.step("est cl baselines start");
		}


    /* CDA */
    for(int c=0; c<k; c++){
      arma::rowvec c_track_vector(k); c_track_vector.zeros(); c_track_vector(c)=1;
			arma::uvec fused_ids = find(theta.row(c) == 0);        /* cl labels fused with cl c */
      arma::uvec other_fused_ids = find(theta.row(c) == 0 and c_track_vector == 0);

      //int fused_n_k = accu(n_k(fused_ids));             /* n_k = number of samples in fused cluster (if fused) */
      int num_fused_cls = fused_ids.n_elem;           /* tracks number of clusters that are fused with current cl c */

      // Have to test this in k>1 setting
      /*fused_ids.print("fused_ids");
      other_fused_ids.print("other_fused_ids");*/

			arma::uvec not_fused_ids = find(theta.row(c) != 0);
			arma::uvec current_fused_and_unfused_ids = find(theta.row(c) != 0 or c_track_vector == 1);

			arma::mat Xk_fused = X.cols(0,k-1);    // Xk_fused: (N rows, k-(num_fused_cls)+1 columns)

			Xk_fused.col(c) = sum(Xk_fused.cols(fused_ids),1);    // add columns of fused cl columns of X
			arma::uvec ids_c = find(Xk_fused.col(c) % keep == 1); // keep just samples with PP > 1e-3
			arma::vec Xc = Xk_fused.col(c);  // extract c'th column of fused X before shedding column (throws off indices)
			arma::vec Xc_ids_c = Xc.rows(ids_c);   // subset Xc to 'keep' rows

			arma::mat Xk_fusedc = Xk_fused.cols(not_fused_ids);   // Xk_fused: (N rows, k-(num_fused_cls) columns)
			arma::vec betac = beta.elem(not_fused_ids);           // subset beta to unfused betas for calculating resids later

			Xk_fused.shed_cols(other_fused_ids);        // remove all other fused cl in Xfused except current one
			arma::vec beta_fused = beta.elem(current_fused_and_unfused_ids);

			if(num_fused_cls>1){
			  int min_fused_id = fused_ids.index_min();       /* determine which of the fused clusters has the smallest label */
				/* skip beta and phi (cl-disp) calculations if c is not the first of the fused cls */
				if(c != min_fused_id){
					beta(c) = beta(min_fused_id);
					phi_j(c) = phi_j(min_fused_id);
					continue;
				}
			}
			arma::vec beta_fusedc = beta.elem(other_fused_ids);

			//initialize the size of matrices/vecs
			arma::mat X_fused(N,(k-num_fused_cls+1)+p), X_fusedc(N,k-num_fused_cls+p);
			arma::vec b_fused((k-num_fused_cls+1)+p), b_fusedc(k-num_fused_cls+p);
			if(p>0){
			  X_fused = join_rows(Xk_fused,X.cols(k,k+p-1)); X_fusedc = join_rows(Xk_fusedc,X.cols(k,k+p-1));
			  b_fused = join_cols(beta_fused,gamma); b_fusedc = join_cols(betac,gamma);
			} else{
			  X_fused = Xk_fused; X_fusedc = Xk_fusedc;
			  b_fused = beta_fused; b_fusedc = betac;
			}

			//b_fused.print("b_fused");   // p+k-(all other fused cls except current)
			//b_fusedc.print("b_fusedc"); // p+k-(#fused cls)
			//Rprintf("X_fused size:%d",X_fused.size());
			//Rprintf("X_fusedc size:%d",X_fusedc.size());

			arma::vec resid_c = y_tilde(ids_c) - X_fusedc.rows(ids_c)*b_fusedc;
			full_resid = y_tilde-X_fused*b_fused;

			//resid_c.print();
			//Rprintf("resid_c size:%d",resid_c.size());  // n*k
			//full_resid.print();
			//Rprintf("full_resid size:%d",full_resid.size());   // n*k

      /* Update beta */
      if(continue_beta==1){
				if(optim_method == "direct"){
				  beta(c) = ((1-alpha)*lambda*((accu(beta)-beta(c))+accu(theta.row(c))) + accu(vec_W(ids_c) % Xc_ids_c % resid_c)/(N) )  /
					((1-alpha)*lambda*(k-1) + accu(vec_W(ids_c) % Xc_ids_c)/(N) );
				} else if(optim_method == "GD"){
					beta(c) = beta(c) - 1 * (-accu(vec_W(ids_c) % Xc_ids_c % full_resid(ids_c))/N + lambda*(1-alpha)*(k*beta(c)-accu(beta)-accu(theta.row(c))));
				} else if(optim_method == "NR"){
					beta(c) = beta(c) - (-accu(vec_W(ids_c) % Xc_ids_c % full_resid(ids_c))/N + lambda*(1-alpha)*(k*beta(c)-accu(beta)-accu(theta.row(c))))  /
					(accu(vec_W(ids_c) % pow(Xc_ids_c,2))/N + lambda*(1-alpha)*(k-1));
				}

				//Rprintf("c=%d,beta(c)=%f",c,beta(c));

				if(beta(c) < (-100)){
					/* Rprintf("Cluster %d, gene %d truncated at -100",c+1,j); */
					beta(c) = -100;
				} else if(beta(c)>100){
					/* Rprintf("Cluster %d, gene %d truncated at +100",c+1,j); */
					beta(c) = 100;
				}
				/* if beta is NaN */
				if(beta(c) != beta(c)){
          beta(c) = log2(mean(y_j));
					Rprintf("Beta of cl%d is NaN. Substituting with mean of other cls\n",c+1);
				}

				b = join_cols(beta,gamma);
      }

    }

    if(i==TIME_ITER){
      timer.step("est cl baselines end");
    }

    /* Hard coded coordinate-wise estimation of gamma */
    if(est_covar==1 && continue_gamma==1){

      /* Hard coded coordinate-wise covar_beta: */
      for(int pp=0; pp<p; pp++){
        /* Xpp: X with p'th col removed, betapp: beta with p'th elt removed */
        arma::vec gammapp=gamma;

        arma::mat X_gammapp = X.cols(k,k+p-1);   // just covariates portion of design matrix
        X_gammapp.shed_col(pp);
        gammapp.shed_row(pp);                   // shed current covar effect (for resid calc)

        arma::mat Xpp = join_rows(X.cols(0,k-1),X_gammapp);
        arma::mat bpp = join_cols(beta,gammapp);      // re-concatenate X and beta without curr covar eff

        // resid = y_tilde - sum_(j!=p){x_ij*beta_j}. full_resid = y_tilde - sum_(j){x_ij*beta_j}
        resid = y_tilde - Xpp*bpp;
        full_resid = y_tilde - X*b;

        arma::vec Xcolpp=X.col(k+pp);    // column of X corr to pp'th covariate

        if(optim_method == "direct"){
          gamma(pp) = accu(vec_W(ids) % Xcolpp(ids) % resid(ids))/accu(vec_W(ids) % pow(Xcolpp(ids),2));
        } else if(optim_method == "GD"){
          gamma(pp) = gamma(pp) - 2 * (-accu(vec_W(ids) % Xcolpp(ids) % full_resid(ids))/n);
        } else if(optim_method == "NR"){
          gamma(pp) = gamma(pp) - (-accu(vec_W(ids) % Xcolpp(ids) % full_resid(ids)))/(accu(vec_W(ids) % pow(Xcolpp(ids),2)));
        }

        /*Rprintf("covar iter%d gamma=%f, top=%f, bottom=%f\n",i,gamma(pp),accu(subs_vec_W % subs_Xcolpp % subs_resid),accu(subs_vec_W % pow(subs_Xcolpp,2)));*/
        if(gamma(pp)>100){
          gamma(pp)=100;
        } else if(gamma(pp)<-100){
          gamma(pp)=-100;
        }   /* is this feasible to do? 2^100 is 1E30, 2^50 is 1.1E15 */

        /* if gamma is NaN */
        if(gamma(pp) != gamma(pp)){
          gamma(pp) = 1e-6;     // VERY small effect for gamma if no convergence --> not 0 for IRLS stopping condition.
          Rprintf("Gamma%d is NaN. Substituting with zero\n",pp+1);
        }

        b = join_cols(beta,gamma);
      }
    }
    if(i==TIME_ITER){
      timer.step("cov ests");
    }

    if(i>=1){
      double diff_beta=0;
			double diff_gamma=0;
      for(int cc=0; cc<k; cc++){
        diff_beta += fabs(temp_beta(i,cc)-temp_beta(i-1,cc))/fabs(temp_beta(i-1,cc)*k);
      }
			for(int cc=0; cc<p; cc++){
				if(p>0){
					diff_gamma += fabs(temp_gamma(i,cc)-temp_gamma(i-1,cc))/fabs(temp_gamma(i-1,cc)*p);
				}
			}

      if(diff_beta<IRLS_tol){
        continue_beta=0;
      }
			if(p>0){
				if(diff_gamma<IRLS_tol){
					continue_gamma=0;
				}
			}else{continue_gamma=0;}


      if(i==maxit_IRLS-1){
        break;
      }
      if(continue_beta==0 && continue_gamma==0){
        break;
      }

		  if(i==TIME_ITER){
			  timer.step("break conds");
		  }

		  if(i==TIME_ITER && output_timer==1){
			  Rcpp::NumericVector res(timer);
			  Rcpp::print(res);
		  }
    }
  }  // end of IRLS

    /* Estimate phi (TURNED OFF. done via theta.ml in R)*/
  /*if(est_phi==1 && cl_phi==0){
      // phi_ML NR estimation
    phi_g = phi_ml_g(y_j(ids),mu(ids),vec_wts(ids),10,phi_g,0);
    if(phi_g>100){
      phi_g=100;
    }
    for(int c=0; c<k; c++){
      phi_j(c) = phi_g;
    }
  } else if(cl_phi==1 && est_phi==1){
    for(int c=0; c<k; c++){
      phi_j(c) = phi_ml(y_j(ids_c),mu(ids_c),vec_wts(ids_c),10,phi_j(c),0);
    }
  }*/
  /*// phi BFGS estimation
   arma::vec temp_yy = y_j(ids);
   arma::vec temp_muu = mu(ids);
   arma::vec temp_wtss = vec_wts(ids);
   const MapVec yy = Rcpp::as<MapVec>(Rcpp::wrap(temp_yy));
   const MapVec muu = Rcpp::as<MapVec>(Rcpp::wrap(temp_muu));
   const MapVec wtss = Rcpp::as<MapVec>(Rcpp::wrap(temp_wtss));
   double fopt;
   phiBFGS nll(yy, muu, wtss);
   Eigen::VectorXd vxd_phi_g(1);
   vxd_phi_g[0]=phi_g;
   int status = optim_lbfgs(nll, vxd_phi_g, fopt);
   if(status < 0)
   Rcpp::stop("fail to converge");
   phi_g = vxd_phi_g[0];*/


    return Rcpp::List::create(Rcpp::Named("coefs_j")=b,
                              Rcpp::Named("theta_j")=theta,
                              Rcpp::Named("temp_beta")=temp_beta,
                              Rcpp::Named("temp_gamma")=temp_gamma); // removed Rcpp::Named("phi_j")=phi_j

}
