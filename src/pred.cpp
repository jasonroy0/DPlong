#include <RcppArmadillo.h>
#include <Rmath.h>
#include "util.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// Xonly      - matrix of X for observed y
// Xonly2     - matrix of X for unobserved y
// h0x        - matrix of values for predictions
// betaY      - current beta parameters
// sig2       - current sigma parameters
// breg       - current spline parameters
// xpipars, xmupars, xsigpars - current X parameters
// ptrt       - number of treatment vars
// p1         - number of discrete vars
// p2         - number of continuous vars
// alpha      - current value of alpha
// Sy         - vector of y clusters
// uniqueS    - unique combinations of clusters
// beta0, prec0, beta_a0, beta_b0, a0_b, b0_b - some priors
// timepoint  - time where predictions wifor subjects with y values
// timepoint2 - time where predictions for subjects without y values
// tZ         - spline of timepoint
// tZ2        - spline of timepoint2

// [[Rcpp::export]]
List pred(mat Xonly, Nullable<mat> Xonly2b, Nullable<mat> h0xb,
	  mat betaY, vec sig2, Nullable<mat> bregb,
	  Nullable<mat> xpipars2, Nullable<mat> xmupars2, Nullable<mat> xsigpars2,
	  int ptrt, int p1, int p2,
	  double alpha,
	  ivec Sy, ivec uniqueS,
	  vec beta0, mat prec0, double beta_a0, double beta_b0, double a0_b, double b0_b,
	  vec timepoint, Nullable<vec> timepoint2b, Nullable<mat> tZb, Nullable<mat> tZ2b) {
  
	bool spline_exists = bregb.isNotNull(); // are there spline values?
	bool newdata_exists = Xonly2b.isNotNull(); // are there predictions from new data?

	
	mat Xonly2, h0x, breg, tZ, tZ2, xpipars, xmupars, xsigpars;
	vec timepoint2;

	// cast nullable back to original types

  if (p1 + ptrt > 0) {
    xpipars = as<mat>(xpipars2);
  }
	
  if (p2 > 0) {
    xmupars  = as<mat>(xmupars2);
    xsigpars = as<mat>(xsigpars2);
  }

	if (spline_exists) {
	  breg = as<mat>(bregb);
	  tZ   = as<mat>(tZb);
	}
	
	if (newdata_exists) {
	  Xonly2      = as<mat>(Xonly2b);
	  h0x         = as<mat>(h0xb);
	  timepoint2  = as<vec>(timepoint2b);
	  
	  if (spline_exists) {
	    tZ2 = as<mat>(tZ2b);
	  }
	}
		
  // number of observations
  int nobs = Xonly.n_rows;
  int nobs2;
	
  // containers to hold predictions - 1 per person
  vec pred1, pred2;
  pred1.zeros(nobs);
  
  ivec pred2_clust;
	if (newdata_exists) {
  	nobs2 = Xonly2.n_rows;
		pred2.zeros(nobs2);
  
		//containers for predicted cluster
  	pred2_clust.zeros(nobs2);
	}

  // numY = number of unique y clusters
  // numBeta = number of beta parameters
  // numTotalCluster = number of total cluster combinations
  int numY, numBeta;
  //int numTotalCluster;


  numY = betaY.n_rows;
  numBeta = betaY.n_cols;
  //numTotalCluster = xmupars.n_rows;

  // vector with number of X clusters pertaining to each Y cluster
  ivec numXCluster; numXCluster.zeros(numY); 
  ivec cumXCluster; // cumulative number of clusters

  // fills in numXCluster -- check for accuracy
  uvec vdummy;
  int dummy;
  
  vec currX;
  
  //predictions for those with data
  for(int i = 0; i < nobs; i++) {
    currX = trans(Xonly.row(i));
    currX.insert_rows(0, 1);
    currX(0) = 1.0;
    
    if (!spline_exists) {
      currX.insert_rows(currX.size(), 1);
      currX(currX.size() - 1) = timepoint(i);
    }

		if (spline_exists) {
    	pred1(i) = R::rnorm( dot( betaY.row( Sy(i) - 1), currX) +
			 					 dot( breg.row(  Sy(i) - 1), tZ.row(i)   ),
			 					 sqrt(sig2( Sy(i) - 1 )));
		} else {
    	pred1(i) = R::rnorm( dot( betaY.row( Sy(i) - 1), currX),
			 					 sqrt(sig2( Sy(i) - 1 )));
		}

  }
  
	//
	//
	//
	//
	//
	// START PREDICTIONS FOR THOSE WITHOUT DATA
	//
	//
	//
	//
  
	
	if (newdata_exists) {	
		// probability in being in each Y cluster (or a new cluster)
		vec probs; probs.zeros(numY+1);

		// values for new cluster
		rowvec newbeta, newb; newbeta.zeros(numBeta);
		double newsig, newsigb;

		if (spline_exists) {
			newb.zeros(tZ.n_cols);
		}

		int nj, nlj, count, count2, chosenCluster;
		double h0, h02;
		
    //predictions for those without data
		for(int i = 0; i < nobs2; i++) {
			
			// Rprintf("Subject: %d\n",i);
			
			h0 = prod(h0x.row(i)); // product of row for subject
			
      //compute probabilities for each cluster + new cluster
			for(int j = 0; j < numY; j++) {

				vdummy = find( Sy == ( j + 1 ) ); // subjects in y cluster
				nj = vdummy.size();              // number of subjects in y cluster
        h02 = 1.0;
        
        
        if (ptrt + p1 > 0) {
          for(int ii = 0; ii < ptrt + p1; ii++) {
            h02 = h02 * R::dbinom(Xonly2(i,ii),1,xpipars(j,ii),0);
          }
        }
      
        if (p2 > 0) {
          for(int ii = 0; ii < p2; ii++) {
            h02 = h02 * R::dnorm(Xonly2(i,ptrt+p1+ii),xmupars(j,ii),xsigpars(j,ii),0);
          }
        }

        probs(j) =  ( nj / ( alpha + nobs ) ) * h02;
      }


		//Rprintf("DONE PART 1\n");
			probs(numY) = h0 * ( alpha / ( alpha + nobs ) );

			chosenCluster = rmultinomF(probs);
			pred2_clust(i) = chosenCluster;
			probs.zeros();
			
			currX = trans(Xonly2.row(i));
			currX.insert_rows(0, 1);
			currX(0) = 1.0;
			currX.insert_rows(currX.size(), 1);
			currX(currX.size() - 1) =  timepoint2(i);

			if(chosenCluster<=numY) 
				{
					if (spline_exists) {
						pred2(i) = R::rnorm( dot( betaY.row( chosenCluster - 1), currX) +
				 							 dot( breg.row(  chosenCluster - 1), tZ2.row(i)   ),
				 							 sqrt(sig2( chosenCluster - 1 )));
					} else {
						pred2(i) = R::rnorm( dot( betaY.row( chosenCluster - 1), currX),
				 							 sqrt(sig2( chosenCluster - 1 )));
					}
				} else {
					newsig = 1/R::rgamma(beta_a0, 1 / beta_b0);
					// newbeta = trans(mvrnorm(beta0,pow(newsig,2)*prec0));
					newbeta = trans(mvrnorm(beta0, newsig*prec0));
					if (spline_exists) {
						newsigb = 1/R::rgamma(a0_b, 1 / b0_b);
						for(int ii = 0; ii < tZ2.n_cols; ii++) {
							newb(ii) = R::rnorm(0.0, sqrt(newsigb));
						}
					
						pred2(i) = R::rnorm( dot( newbeta, currX) +
				 							 dot( newb   , tZ2.row(i)   ),
				 							 sqrt( newsig ));
					} else {

						pred2(i) = R::rnorm( dot( newbeta, currX),
				 							 sqrt( newsig ));
					}

				}

		}
		
		return List::create(_["pred1"] = pred1,
												_["pred2"] = pred2,
												_["pred2_clust"] = pred2_clust);

	} else {
  	return List::create(_["pred1"] = pred1);
	}
}

