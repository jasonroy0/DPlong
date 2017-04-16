#include <RcppArmadillo.h>
#include <Rmath.h>
#include "util.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List update_covs(mat matX, mat s, mat uniqueS,
                 Nullable<mat> xpipars2, Nullable<mat> xmupars2, Nullable<mat> xsigpars2,
                 int p1, int p2, int ptrt,
                 double a0, double b0,
                 double mu0, double nu0, double tau0, double c0) {
  
  bool bin_exists  = false;
  bool cont_exists = false;
  
  mat xpipars, xmupars, xsigpars;
  
  if ( xpipars2.isNotNull() ) {
    bin_exists = true;
    xpipars = as<mat>(xpipars2);
  }
  
  if ( xmupars2.isNotNull() ) {
    cont_exists = true;
    xmupars = as<mat>(xmupars2);
    xsigpars = as<mat>(xsigpars2);
  }
  
  uvec xpos(p1 + p2 + ptrt);
  for(unsigned int ii = 0; ii < p1 + p2 + ptrt; ii++) xpos(ii) = ii;
  
  mat tempX; // subsetted matrix of X values for appropriate cluster
  vec tempx; // vector from tempX of single column
  uvec vdummy;
  int  dummy, nx, sizex;
  double newdf, varx, numer, meanx, sumx, curtau, newvar, newmean;
  
  int count = 0;
  
  int k = max(uniqueS.col(0)) - 1; //number of unique Y clusters
  
  for(int j = 0; j < k; j++) {
    vdummy = find(uniqueS.col(0) == (j + 1));
    nx = vdummy.n_elem;
    //n.x <- length( unique( s[ s[ , 1 ] == j , 2] ) )  // number of x clusters within jth y cluster
    
    for(int l = 0; l < nx; l++) {
      vdummy = find(s.col(0) == j + 1 && s.col(1) == l + 1);
      tempX   = matX( vdummy, xpos ); 
       
      
      // binary covariates
      if ( bin_exists ) {
        for( int ii = 0; ii < ptrt + p1; ii++ ) {
          tempx   = tempX.col(ii);
          sumx    = sum(tempX.col(ii)); 
          sizex   = tempx.n_elem;
          xpipars( count, ii ) = R::rbeta( sumx + a0 , sizex - sumx + b0 );
        }	
      }
      
      // continuous covariates
      if ( cont_exists ) {
        for( int ii = 0; ii < p2; ii++ ) {
          tempx = tempX.col(p1 + ptrt + ii);
          sizex = tempx.n_elem;
          meanx = mean(tempx);
          newdf = nu0 + sizex;
          varx  = 0;
          if (sizex > 1) varx = var(tempx);
          numer = nu0 * tau0 + sizex * varx + (c0 * sizex / (c0 + sizex) ) * pow(meanx - mu0, 2);
          xsigpars(count, ii) = rinvchisq( newdf, numer / newdf );
          
          // posterior for mu. prior for mu|sigma^2 mean 0 prior sample size 2
          curtau = xsigpars(count, ii);
          newvar = 1 / ( c0 / curtau + sizex / curtau);
          newmean = ( mu0 * c0 / curtau + meanx * sizex / curtau ) * newvar;
          xmupars(count, ii) = R::rnorm( newmean, sqrt(newvar) );
        }
      }
      count++;
    }
  } // end covariate parameter update
  
  if ( bin_exists && !cont_exists) {
    return List::create(_["xpipars"] = xpipars);
  } else if ( !bin_exists && cont_exists) {
    return List::create(_["xmupars"] = xmupars,
                        _["xsigpars"] = xsigpars);
  } else {
    return List::create(_["xpipars"] = xpipars,
                      _["xmupars"] = xmupars,
                      _["xsigpars"] = xsigpars);
  }
  
}

