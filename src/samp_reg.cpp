#include <RcppArmadillo.h>
#include <Rmath.h>
#include "util.h"

// [[Rcpp::depends(RcppArmadillo)]]


using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List samp_reg(vec y, mat Xall,  ivec all_ids, Nullable<mat> Z2, vec u,   //data
                mat betaY, vec sig2, Nullable<mat> b2,         //Y params
                ivec Sy,       //cluster
                vec beta0, mat prec0, double beta_a0, double beta_b0,  //priors on Yparams
                Nullable<vec> sig2_b2, double a0_b, double b0_b  // variance for b
                ) {
 
  bool zexists = false;
  
  // these need to be declared outside the if statements
  mat Z, b;
  vec sig2_b;
  
  if ( Z2.isNotNull() ) {
    zexists = true;
    Z       = as<mat>(Z2);
    b       = as<mat>(b2);
    sig2_b  = as<vec>(sig2_b2);
  }
 
  int k = betaY.n_rows;
  int nbeta = betaY.n_cols;
  mat matX, matZ;
  vec tempu, tempy, resid;
  uvec k_ids, curr_ids;
  uvec xpos(nbeta); 
  
  int nZ;
  uvec zpos;
  mat new_sig;
  vec new_mu;
	
	if (zexists) {
		nZ = Z.n_cols;
		zpos.set_size(nZ);
  	for(unsigned int ii = 0; ii < nZ; ii++) zpos(ii) = ii;
 	}
	
  //updating regvar
  double an;
  mat bn;
  mat newprec;
  vec betan;
  
  for(unsigned int ii = 0; ii < nbeta; ii++) xpos(ii) = ii;
  
  for(int i = 0; i < k; i++) {
    
    //step 1: ids with Sy = k+1
    k_ids = find(Sy == i + 1);
    int nk = k_ids.n_elem;
    curr_ids = find(all_ids == k_ids(0) + 1);
    for(int ii = 1; ii < nk; ii++) curr_ids = join_cols(curr_ids, find(all_ids == k_ids(ii) + 1));

    //step 2: matX and matZ and y for those ids
    tempy = y(curr_ids);
    matX = Xall(curr_ids,xpos);
    if (zexists) matZ = Z(curr_ids,zpos);
    
    //step 3: u with those ids
    tempu = u(curr_ids);

    //step 4: resids Y - Z - u
    if (zexists) {
			resid = tempy - tempu - matZ * trans(b.row(i));
		} else {
			resid = tempy - tempu;
   	}


    //step 5: update sigreg
    an = beta_a0 + tempy.n_elem/2;
    newprec = matX.t()*matX + prec0;
    betan = newprec.i()*(prec0*beta0 + matX.t()*resid);
    bn = beta_b0 + 0.5 * (resid.t()*resid + beta0.t()*prec0*beta0 - betan.t()*newprec*betan);
    //Rprintf("an: %.2f\n",an);
    sig2(i) = 1/R::rgamma(an,1/bn(0,0));
    
    //step 6: update beta
    betaY.row(i) = trans(mvrnorm(betan,sig2(i) * newprec.i()));

    if (zexists) {
    	//step 7: update sig_b
			sig2_b(i) = 1/R::rgamma(a0_b + nZ/2, 1/(b0_b + 0.5*accu(pow(b.row(i),2))));
    
    	//step 8: resids Y - X - u
    	resid = tempy - matX * trans(betaY.row(i)) - tempu;

    	//step 9: update b
    	new_sig = matZ.t()*matZ / sig2(i) + eye(nZ, nZ) / sig2_b(i);
    	new_mu = new_sig.i()*matZ.t()*resid / sig2(i);
    	b.row(i) = trans(mvrnorm(new_mu,new_sig.i()));

    }

  }

  if (zexists) {
		return List::create(_["betaY"] = betaY,
    	                  _["sig2"] = sig2,
      	                _["b"] = b,
        	              _["sig2_b"] = sig2_b);
	} else {
		return List::create(_["betaY"] = betaY,
    	                  _["sig2"] = sig2);

	}

}

