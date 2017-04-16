#include <RcppArmadillo.h>
#include <Rmath.h>
#include "util.h"

// [[Rcpp::depends(RcppArmadillo)]]


using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
List samp_u(vec y, mat Xall, Nullable<mat> Z2,
            mat betaY, Nullable<mat> b2,
            vec ids, ivec Sy,
            double sig2_u, vec sig2,
            int n) {
  
  int nbeta = betaY.n_cols;
  uvec currID;
  uvec xpos(nbeta);
  vec tempy;
  mat matX;
  int nobs;
  double newsig, newmu;
  NumericVector u(n);//  u.set_size(n); u.zeros();


  for(unsigned int ii = 0; ii < nbeta; ii++) xpos(ii) = ii;

	bool zexists = Z2.isNotNull();
  int nZ;
  uvec zpos;
	mat Z, matZ, b;
	
	if (zexists) {
	  Z = as<mat>(Z2);
	  b = as<mat>(b2);
  	nZ = Z.n_cols;
		
		zpos.set_size(nZ);
  	for(unsigned int ii = 0; ii < nZ; ii++) zpos(ii) = ii;
	}

  for(int i = 0; i < n; i++) {
    currID = find(ids == i+1);
    nobs = currID.n_elem;
    matX = Xall(currID, xpos);
    if (zexists) {
			matZ = Z(currID, zpos);
    	tempy = y(currID) - matX * trans(betaY.row(Sy(i)-1)) - matZ * trans(b.row(Sy(i)-1));
		} else {
    	tempy = y(currID) - matX * trans(betaY.row(Sy(i)-1));
		}
    newsig = (nobs*sig2_u+sig2(Sy(i)-1))/(sig2_u*sig2(Sy(i)-1));
    newmu = accu(tempy)*sig2_u/(nobs*sig2_u+sig2(Sy(i)-1));

    u(i) = R::rnorm(newmu,sqrt(1/newsig));

  }
  
  return List::create(_["u"] = u);
  
}
