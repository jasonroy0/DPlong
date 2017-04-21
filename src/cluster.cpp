#include <RcppArmadillo.h>
#include <Rmath.h>
#include "util.h"

// todo: clean function -- write smaller functions for each part

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List cluster(vec y, mat Xonly, mat Xall, vec ntp, vec ids, Nullable<mat> Z2, vec u,   //data
						 mat betaY, vec sig2, Nullable<mat> b2,         //Y params
						 Nullable<mat> xpipars2, Nullable<mat> xmupars2, Nullable<mat> xsigpars2, //X params
						 int ptrt, int p1, int p2,     //# vars
						 double alpha,  //alpha
						 ivec Sy, ivec uniqueS,       //cluster
						 vec beta0, mat prec0, double beta_a0, double beta_b0,  //priors on Yparams
						 Nullable<vec> sig2_b2, double a0_b, double b0_b,  // variance for b
						 double a0, double b0, // priors on discrete X vars
						 double mu0, int nu0, double tau0, double c0, //priors on cont X vars
						 double alp_a0, double alp_b0, //priors on alpha
						 int m) {
 	
  mat xpipars, xmupars, xsigpars;

  if (p1 + ptrt > 0) {
    xpipars = as<mat>(xpipars2);
  }

  if (p2 > 0) {
    xmupars  = as<mat>(xmupars2);
    xsigpars = as<mat>(xsigpars2);
  }

  bool spl = Z2.isNotNull();

    
  mat Z, b;
  vec sig2_b;
  int nknots;

  if (spl) { 
    Z      = as<mat>(Z2);
    b      = as<mat>(b2);
    sig2_b = as<vec>(sig2_b2);

		nknots = b.n_cols;
  }

  int nobs = Xonly.n_rows;
  int nbeta = betaY.n_cols;
		
  int numY;
		
  int dummy1, dummy2;
  uvec vdummy1, vdummy2;
		
  int newCluster;
		
  double likeregy, prodx, prodx2, likeregb;
		
  // containers for auxiliary parameter values
  mat xpipars_aux;
  mat xmupars_aux;
  mat xsigpars_aux;
		
  mat betaY_aux;
  vec sig2_aux;
  mat b_aux;
  vec sig2_b_aux;
  mat Y_xpipars;
  mat Y_xmupars;
  mat Y_xsigpars;
  
  uvec currID;
  uvec xpos(nbeta);  
  for(unsigned int ii = 0; ii < nbeta; ii++) xpos(ii) = ii;
  uvec xpos2;
  if (spl) {
    xpos2.set_size(nknots); 
    for(unsigned int ii = 0; ii < nknots; ii++) xpos2(ii) = ii;
  }
  vec yvals;
  mat xvals;
  mat zvals;
		
	 
  for(int j = 1; j < nobs; j++) {
    //Rprintf("\n\nSubject: %d\n",j);
    currID = find(ids == j+1);
    yvals = y(currID);
    xvals = Xall(currID,xpos);
    if (spl) {
      zvals = Z(currID,xpos2);
    }
    
    // vector of people in same cluster as person j
    vdummy1 = find(Sy == Sy(j)); 
    dummy1 = vdummy1.size(); //# in the cluster
    
    //if lone person in cluster
    if(dummy1==1) { //if lone person in X-Y cluster
      
    
      betaY.shed_row(Sy(j)-1);
      sig2.shed_row(Sy(j)-1);
      if (spl) {
        b.shed_row(Sy(j)-1);
        sig2_b.shed_row(Sy(j)-1);
      }
      
      // delete X coef
      //should find row in uniqueS that corresponds to person i
      vdummy1 = find(uniqueS==Sy(j));
      
      if (p1 + ptrt > 0) {
        xpipars.shed_row(vdummy1(0)); 
      }
      
      if (p2 > 0) {
        xmupars.shed_row(vdummy1(0));
        xsigpars.shed_row(vdummy1(0)); 
      }
      
      //relabel Y cluster (if needed)
      if(dummy2==1) {
        for(int k = 0; k < Sy.size(); k++) {
          if(Sy(k) > Sy(j))
            Sy(k) = Sy(k) - 1;
        }
					
        for(int k = 0; k < uniqueS.n_elem; k++) {
          if(uniqueS(k) > Sy(j))
            uniqueS(k) = uniqueS(k) - 1;
        }
      }
				
      uniqueS.shed_row(vdummy1(0)); //get rid of row
    }  //end lone person in cluster
			
    
    Sy.shed_row(j);
			
			
    numY = betaY.n_rows;

    int totalposs = numY + m; 
    vec probs(totalposs);
    int count=0;
			
    int njwoi, nljwoi; //counts for # in appropriate Y and X cluster,
    //excluding the ith person

    double spline_part = 0.0;

    //  Rprintf("Calculate probabilities for existing clusters\n");
    for(int k = 0; k < numY; k++) {
      //FILL IN PROBS FOR EXISTING CLUSTERS
      //get count of number of X clusters within kth Y cluster
      
      //get number of subjects within kth cluster
      vdummy1 = find(Sy == (k+1));
      njwoi = vdummy1.size();
      
      likeregy = 1;
      
      for(int kk = 0; kk < ntp(j); kk++) {
        if (spl) { 
          spline_part = dot(zvals.row(kk), b.row(k));
        }
        likeregy = likeregy * R::dnorm(yvals(kk),dot(xvals.row(kk),betaY.row(k)) + spline_part + u(j), sqrt(sig2(k)), 0);
      }
      
      likeregb = 1;
      for(int kk = 0; kk < nknots; kk++) {
        //likeregb = likeregb * R::dnorm(b(k,kk),0,sig2_b(k),0);
      }
     
      prodx = 1;
      prodx2 = 1;
      
      
      //likelihood for binary covariates
      if (p1 + ptrt > 0) {
        for(int q = 0; q < ptrt+p1; q++) {
          prodx = prodx*R::dbinom(Xonly(j,q), 1, xpipars(k, q), 0);
        }
      }

      //likelihood for continuous covariates
      if (p2 > 0) {
        for(int q = 0; q < p2; q++) {
          prodx2 = prodx2*R::dnorm(Xonly(j,ptrt+p1+q), xmupars(k,q), sqrt(xsigpars(k,q)), 0 );
          //Rprintf("prodx2: %.8f\n", prodx2);
        }
      }
        
      probs(k) = ( njwoi / (nobs - 1 + alpha) ) * likeregy * prodx * prodx2 * likeregb;

      
    } // ends probs for existing clusters
			
			
    
    
    
    // COMPLETELY NEW CLUSTERS
    
    //set auxiliary parameters
    betaY_aux.set_size(m,nbeta); betaY_aux.zeros();
    sig2_aux.set_size(m); sig2_aux.zeros();
    if (spl) {
      b_aux.set_size(m,nknots); b_aux.zeros();
      sig2_b_aux.set_size(m); sig2_b_aux.zeros();
    }

    if (p1 + ptrt > 0) {
      Y_xpipars.set_size(m,p1+ptrt); Y_xpipars.zeros();
    }

    if (p2 > 0) {
      Y_xmupars.set_size(m,p2); Y_xmupars.zeros();
      Y_xsigpars.set_size(m,p2); Y_xsigpars.zeros();
    }

    for(int w = 0; w < m; w++) {
      
      sig2_aux(w) = 1/R::rgamma(beta_a0, 1 / beta_b0);
      //betaY_aux.row(w) = trans(mvrnorm(beta0,sig2_aux(w)*prec0));
      betaY_aux.row(w) = trans(mvrnorm(beta0,sig2_aux(w)*prec0.i()));
      
      if (spl) {
        sig2_b_aux(w) = 1/R::rgamma(a0_b,1 / b0_b);
        
        for(int ww = 0; ww < nknots; ww++) {
          b_aux(w,ww) = R::rnorm(0,sqrt(sig2_b_aux(w)));
        }
      }
    
      if (p1 + ptrt > 0) {
        for(int ww = 0; ww < ptrt + p1; ww++) {
          Y_xpipars(w,ww) = R::rbeta(a0,b0);
        }
      }
      
      if (p2 > 0) {
        for(int ww = 0; ww < p2; ww++) {
          Y_xsigpars(w,ww) = rinvchisq(nu0,tau0);
          Y_xmupars(w,ww) = R::rnorm(mu0,xsigpars_aux(w,ww)/sqrt(c0));
        }
      }
    }

    //Rprintf("Calculate probabilities for new clusters\n");
    for(int k = 0; k < m; k++) {
      
      likeregy = 1;
      
      for(int kk = 0; kk < ntp(j); kk++) {
        if (spl) { 
          spline_part = dot(zvals.row(kk), b_aux.row(k));
        }
        likeregy = likeregy * R::dnorm(yvals(kk), dot(xvals.row(kk), betaY_aux.row(k)) + spline_part + u(j), sqrt(sig2_aux(k)), 0);
      }
      
      likeregb = 1;
      for(int kk = 0; kk < nknots; kk++) {
        //likeregb = likeregb * R::dnorm(b_aux(k,kk),0,sig2_b_aux(k),0);
      }

      prodx = 1;
      prodx2 = 1;
      
      //likelihood for binary covariates
      if (p1 + ptrt > 0) {
        for(int q = 0; q < ptrt+p1; q++) {
          prodx = prodx*R::dbinom(Xonly(j,q), 1, Y_xpipars(k, q), 0);
        }
      }
      
      //likelihood for continuous covariates
      if (p2 > 0) {
        for(int q = 0; q < p2; q++) {
          prodx2 = prodx2*R::dnorm(Xonly(j,ptrt+p1+q), Y_xmupars(k,q), sqrt(Y_xsigpars(k,q)), 0 );
        }
      }
      
      probs(numY + k) = (alpha / m) / (nobs - 1 + alpha) * prodx * prodx2 * likeregy * likeregb;
    
    }
    
    // Rcout << "probs: " << probs.t() << std::endl;
    
    
    //USE MULTINOMIAL DISTRIBUTION TO CHOOSE NEW CLUSTER
    newCluster = rmultinomF(probs);
    //Rprintf("Chosen cluster: %d\n",newCluster);
    
    probs.zeros();
    
    
    
    if(newCluster <= numY) 
    {
      //Rprintf("Chose existing cluster\n");
      Sy.insert_rows(j,1);
      Sy(j) = uniqueS(newCluster-1);
    }
    else //new Y and X cluster
    {
      //Rprintf("Chose new cluster\n");
      
      Sy.insert_rows(j,1);
      Sy(j) = Sy.max()+1;
      //Rprintf("new Y cluster number: %d\n",Sy(j));
      
      uniqueS.insert_rows(numY,1);
      uniqueS(numY) = Sy(j);
      betaY.insert_rows(numY,1);
      betaY.row(numY) = betaY_aux.row(newCluster- numY -1);
      if (spl) {
        b.insert_rows(numY, 1);
        b.row(numY) = b_aux.row(newCluster - numY - 1);
        sig2_b.insert_rows(numY, 1);
        sig2_b.row(numY) = sig2_b_aux.row(newCluster - numY - 1);
      }
      sig2.insert_rows(numY, 1);
      sig2.row(numY) = sig2_aux.row(newCluster - numY - 1);
      if (p1 + ptrt > 0) {
        xpipars.insert_rows(numY,1);
        xpipars.row(numY) = Y_xpipars.row(newCluster - numY - 1);
      }
      if (p2 > 0) {
        xmupars.insert_rows(numY,1);
        xmupars.row(numY) = Y_xmupars.row(newCluster - numY - 1);
        xsigpars.insert_rows(numY,1);
        xsigpars.row(numY) = Y_xsigpars.row(newCluster- numY-1);
      }
    }
    
    //Rprintf("Y cluster: %d\n",Sy(j));
    
  } // end cluster loop (j)
  
  if (spl) {
    return List::create(_["uniqueS"] = uniqueS,
                         _["Sy"] = Sy,
                         _["xpipars"] = xpipars,
                         _["xmupars"] = xmupars,
                         _["xsigpars"] = xsigpars,
                         _["betaY"] = betaY,
                         _["sig2"] = sig2,
                         _["b"] = b,
                         _["sig2_b"] = sig2_b);
   } else {
    return List::create(_["uniqueS"] = uniqueS,
                         _["Sy"] = Sy,
                         _["xpipars"] = xpipars,
                         _["xmupars"] = xmupars,
                         _["xsigpars"] = xsigpars,
                         _["betaY"] = betaY,
                         _["sig2"] = sig2);
   }
}  
