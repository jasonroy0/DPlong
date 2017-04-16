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

	// first part if splines are included
	if ( Z2.isNotNull() ) {
    
    mat Z, b;
	  vec sig2_b;

	  
	  Z      = as<mat>(Z2);
	  b      = as<mat>(b2);
	  sig2_b = as<vec>(sig2_b2);
    
		int nobs = Xonly.n_rows;
		int nbeta = betaY.n_cols;
		int nknots = b.n_cols;
		
		int numY;
		
		int dummy1, dummy2;
		uvec vdummy1, vdummy2;
		
		int newCluster;
		
		double likeregy, prodx, prodx2, likeregb;
		
		bool onlyX;
		bool onlyY;
		
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
		//uvec xpos(1+p1+p2+ptrt+2);
		//for(unsigned int ii = 0; ii < 1+p1+p2+ptrt+2; ii++) xpos(ii) = ii;
		uvec xpos(nbeta);  //
		for(unsigned int ii = 0; ii < nbeta; ii++) xpos(ii) = ii;
		uvec xpos2(nknots); 
		for(unsigned int ii = 0; ii < nknots; ii++) xpos2(ii) = ii;
		vec yvals;
		mat xvals;
		mat zvals;
		
	 
		for(int j = 1; j < nobs; j++) {
			//Rprintf("\n\nSubject: %d\n",j);
			currID = find(ids == j+1);
			yvals = y(currID);
			xvals = Xall(currID,xpos);
			zvals = Z(currID,xpos2);
		  onlyY = false;
			
			vdummy1 = find(Sy == Sy(j)); 
			dummy1 = vdummy1.size(); //# in the cluster
			
			if(dummy1==1) { //if lone person in X-Y cluster
				onlyX = true;
				
				
			  betaY.shed_row(Sy(j)-1);
				sig2.shed_row(Sy(j)-1);
				b.shed_row(Sy(j)-1);
				sig2_b.shed_row(Sy(j)-1);
        
				if (p1 + ptrt > 0) {
          xpipars.shed_row(Sy(j)-1); // CHECK --- NEEDED TO DO THIS BECAUSE dummy1 is uvec
				}

        if (p2 > 0) {
          xmupars.shed_row(Sy(j)-1);
				  xsigpars.shed_row(Sy(j)-1); 
			  }


        // RELABEL CLUSTERS
        for(int k = 0; k < Sy.size(); k++) {
          if(Sy(k) > Sy(j))
            Sy(k) = Sy(k) - 1;
        }
          
        for(int k = 0; k < uniqueS.n_rows; k++) {
          if(uniqueS(k) > Sy(j))
            uniqueS(k) = uniqueS(k) - 1;
        }
      }
				
			
			// Delete rox in Sy and Sx;
			Sy.shed_row(j); Sx.shed_row(j);
			
			
			numY = betaY.n_rows;

			//Rprintf("Total number of clusters: %d\n",numTotalCluster);
			//Rprintf("Number of Y clusters: %d\n",numY);
			
      //existing clusters plus auxiliary parameters;
			int totalposs = numY + m; 
			
      //Rprintf("Number of possibilities: %d\n",totalposs);
			vec probs(totalposs);
			int count=0;
			
			int njwoi; // counts for # in cluster, exclusing ith person

			//  Rprintf("Calculate probabilities for existing clusters\n");
			
      for(int k = 0; k < numY; k++) {
				
				//get number of subjects within kth cluster
				vdummy1 = find(Sy == (k+1));
				njwoi = vdummy1.size();
				
				likeregy = 1;
				
				for(int kk = 0; kk < ntp(j); kk++) {
					likeregy = likeregy * R::dnorm(yvals(kk),dot(xvals.row(kk),betaY.row(k))+dot(zvals.row(kk),b.row(k))+u(j),sqrt(sig2(k)),0);
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
          }
        }	
        
        probs(k) = (njwoi/(nobs+alpha))*likeregy*prodx*prodx2*likeregb; // FIX THIS
				
			} // ends probs for existing clusters
			
			
			
			
			
			// COMPLETELY NEW CLUSTERS
			
			//set auxiliary parameters
			betaY_aux.set_size(m,nbeta); betaY_aux.zeros();
			sig2_aux.set_size(m); sig2_aux.zeros();
			b_aux.set_size(m,nknots); b_aux.zeros();
			sig2_b_aux.set_size(m); sig2_b_aux.zeros();
			
			if (p1 + ptrt > 0) {
			  Y_xpipars.set_size(m,p1+ptrt); Y_xpipars.zeros();
			}

      if (p2 > 0) {
        Y_xmupars.set_size(m,p2); Y_xmupars.zeros();
        Y_xsigpars.set_size(m,p2); Y_xsigpars.zeros();
      }

			for(int w = 0; w < m; w++) {
				
				sig2_aux(w) = 1/R::rgamma(beta_a0,beta_b0);
				//betaY_aux.row(w) = trans(mvrnorm(beta0,sig2_aux(w)*prec0));
				betaY_aux.row(w) = trans(mvrnorm(beta0,sig2_aux(w)*prec0.i()));
				
				sig2_b_aux(w) = 1/R::rgamma(a0_b,b0_b);
				
				for(int ww = 0; ww < nknots; ww++) {
					b_aux(w,ww) = R::rnorm(0,sqrt(sig2_b_aux(w)));
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
					likeregy = likeregy * R::dnorm(yvals(kk),dot(xvals.row(kk),betaY_aux.row(k))+dot(zvals.row(kk),b_aux.row(k))+u(j),sqrt(sig2_aux(k)),0);
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
				
				probs(numY + k) = (alpha/m)*prodx*prodx2*likeregy*likeregb;
			}
			
			// Rcout << "probs: " << probs.t() << std::endl;
			
			
			//USE MULTINOMIAL DISTRIBUTION TO CHOOSE NEW CLUSTER
			newCluster = rmultinomF(probs);
			//Rprintf("Chosen cluster: %d\n",newCluster);
			
			probs.zeros();
			
			
			
			if(newCluster<=numY) // EXISTING CLUSTER
			{
				Sy.insert_rows(j,1);
				Sy(j) = uniqueS(newCluster-1);
			}
			else                 // NEW CLUSTER
			{
				Sy.insert_rows(j,1);
				Sy(j) = Sy.max()+1;
				
				uniqueS.insert_rows(numY,1);
				uniqueS(numTotalCluster) = Sy(j);
				betaY.insert_rows(numY,1);
				betaY.row(numY) = betaY_aux.row(newCluster-numY-1);
				b.insert_rows(numY,1);
				b.row(numY) = b_aux.row(newCluster-numY-1);
				sig2_b.insert_rows(numY,1);
				sig2_b.row(numY) = sig2_b_aux.row(newCluster-numY-1);
				sig2.insert_rows(numY,1);
				sig2.row(numY) = sig2_aux.row(newCluster-numY-1);
				if (p1 + ptrt > 0) {
          xpipars.insert_rows(numY,1);
          xpipars.row(numTotalCluster) = Y_xpipars.row(newCluster-numY-1);
				}
        if (p2 > 0) {
          xmupars.insert_rows(numY,1);
          xmupars.row(numY) = Y_xmupars.row(newCluster-numY-1);
          xsigpars.insert_rows(numY,1);
          xsigpars.row(numY) = Y_xsigpars.row(newCluster-numY-1);
			  }
      }
			
		} // end cluster loop (j)
		
		
		//Rprintf("Loop succesfully completed.\n");
		return List::create(_["uniqueS"] = uniqueS,
												 _["Sy"] = Sy,
												 _["Sx"] = Sx,
												 _["xpipars"] = xpipars,
												 _["xmupars"] = xmupars,
												 _["xsigpars"] = xsigpars,
												 _["betaY"] = betaY,
												 _["sig2"] = sig2,
												 _["b"] = b,
												 _["sig2_b"] = sig2_b);
	} else {  
  
	// first part if splines are NOT included
  
		
		int nobs = Xonly.n_rows;
		int nbeta = betaY.n_cols;
		
		int numY;
		int numTotalCluster;
		ivec numXCluster;
		
		int dummy1, dummy2;
		uvec vdummy1, vdummy2;
		
		int newCluster;
		
		double likeregy, prodx, prodx2;
		
		bool onlyX;
		bool onlyY;
		
		// containers for auxiliary parameter values
		mat xpipars_aux;
		mat xmupars_aux;
		mat xsigpars_aux;
		
		mat betaY_aux;
		vec sig2_aux;
		mat Y_xpipars;
		mat Y_xmupars;
		mat Y_xsigpars;
		
		uvec currID;
		//uvec xpos(1+p1+p2+ptrt+2);
		//for(unsigned int ii = 0; ii < 1+p1+p2+ptrt+2; ii++) xpos(ii) = ii;
		uvec xpos(nbeta);  //
		for(unsigned int ii = 0; ii < nbeta; ii++) xpos(ii) = ii;
		vec yvals;
		mat xvals;
		
	 
		for(int j = 1; j < nobs; j++) {
			//Rprintf("\n\nSubject: %d\n",j);
			currID = find(ids == j+1);
			yvals = y(currID);
			xvals = Xall(currID,xpos);
			onlyX = false; onlyY = false;
			
			// vector of people in same cluster as person j
			vdummy1 = find(Sy == Sy(j) && Sx == Sx(j)); 
			dummy1 = vdummy1.size(); //# in the cluster
			
				//if lone person in cluster
			if(dummy1==1) { //if lone person in X-Y cluster
				//DELETE ASSOCIATED COEFFICIENTS IN Y AND X CLUSTER
				onlyX = true;
				
				vdummy2 = find(Sy == Sy(j)); //check if only person in Y cluster too
				dummy2 = vdummy2.size();
				
				// delete Y coef if only one in Y cluster
				if(dummy2==1) {
					onlyY = true;
					betaY.shed_row(Sy(j)-1);
					sig2.shed_row(Sy(j)-1);
				}
				
				// delete X coef
				//should find row in uniqueS that corresponds to person i
				vdummy1 = find(uniqueS.col(0)==Sy(j) && uniqueS.col(1)==Sx(j));
				
				if (p1 + ptrt > 0) {
          xpipars.shed_row(vdummy1(0)); // CHECK --- NEEDED TO DO THIS BECAUSE dummy1 is uvec
				}
        
        if (p2 > 0) {
          xmupars.shed_row(vdummy1(0));
				  xsigpars.shed_row(vdummy1(0)); 
				}
				
				//relabel X cluster
				for(int k = 0; k < Sx.size(); k++) {
					if(Sy(k) == Sy(j) && Sx(k) > Sx(j))
						Sx(k) = Sx(k) - 1;
				}
				
				for(int k = 0; k < uniqueS.n_rows; k++) {
					if(uniqueS(k,0) == Sy(j) && uniqueS(k,1) > Sx(j))
						uniqueS(k,1) = uniqueS(k,1) - 1;
				}
				
				//relabel Y cluster (if needed)
				if(dummy2==1) {
					for(int k = 0; k < Sy.size(); k++) {
						if(Sy(k) > Sy(j))
							Sy(k) = Sy(k) - 1;
					}
					
					for(int k = 0; k < uniqueS.n_rows; k++) {
						if(uniqueS(k,0) > Sy(j))
							uniqueS(k,0) = uniqueS(k,0) - 1;
					}
				}
				
				uniqueS.shed_row(vdummy1(0)); //get rid of row
			}  //end lone person in cluster
			
			// Delete rox in Sy and Sx;
			Sy.shed_row(j); Sx.shed_row(j);
			
			
			numY = betaY.n_rows;
			if (p2 > 0) {
        numTotalCluster = xmupars.n_rows;
      } else {
        numTotalCluster = xpipars.n_rows;
      }
      
			//Rprintf("Total number of clusters: %d\n",numTotalCluster);
			//Rprintf("Number of Y clusters: %d\n",numY);
			//existing clusters plus auxiliary parameters;
			int totalposs = numY*m+numTotalCluster+m; 
			//Rprintf("Number of possibilities: %d\n",totalposs);
			vec probs(totalposs);
			int count=0;
			
			int njwoi, nljwoi; //counts for # in appropriate Y and X cluster,
			//excluding the ith person
			int numXj; //# of X clusters within each Y cluster

			//  Rprintf("Calculate probabilities for existing clusters\n");
			for(int k = 0; k < numY; k++) {
				//FILL IN PROBS FOR EXISTING CLUSTERS
				//get count of number of X clusters within kth Y cluster
				vdummy1 = find(uniqueS.col(0) == (k+1));
				numXj = vdummy1.size();
				
				//get number of subjects within kth cluster
				vdummy1 = find(Sy == (k+1));
				njwoi = vdummy1.size();
				
				likeregy = 1;
				
				for(int kk = 0; kk < ntp(j); kk++) {
					likeregy = likeregy * R::dnorm(yvals(kk),dot(xvals.row(kk),betaY.row(k))+u(j),sqrt(sig2(k)),0);
				}
				
			 
				for(int l = 0; l < numXj; l++) {
					prodx = 1;
					prodx2 = 1;
					
					vdummy1 = find(Sy == (k+1) && Sx == (l+1));
					nljwoi = vdummy1.size();
					
					//likelihood for binary covariates
          if (p1 + ptrt > 0) {
            for(int q = 0; q < ptrt+p1; q++) {
              prodx = prodx*R::dbinom(Xonly(j,q), 1, xpipars(count, q), 0);
            }
          }
					
					//likelihood for continuous covariates
          if (p2 > 0) {
            for(int q = 0; q < p2; q++) {
              prodx2 = prodx2*R::dnorm(Xonly(j,ptrt+p1+q), xmupars(count,q), sqrt(xsigpars(count,q)), 0 );
            }
					}

					probs(count) = ((njwoi*nljwoi)/(njwoi+alphapsi))*likeregy*prodx*prodx2;
					if(probs(count)<0) {
						Rprintf("NEGATIVE PROBABILITY ALERT\n");
						Rprintf("prodx = %.8f\n",prodx);
						Rprintf("vals: %.1f, %.1f, %.1f\n",Xonly(j,0),Xonly(j,1),Xonly(j,2));
						Rprintf("prodx2 = %.8f\n",prodx2);
						Rprintf("vals: %.2f, %.2f\n",Xonly(j,ptrt+p1),Xonly(j,ptrt+p1+1));
						Rprintf("likeregy = %.8f\n",likeregy);
					}
					//Rprintf("count + prob: %.2f, %d\n",count,probs(count));
					count++;
				}
				
			} // ends probs for existing clusters
			
			
			
			
			// EXISTING Y CLUSTER BUT NEW X CLUSTERS
			//set auxiliary parameters
			if (p1 + ptrt > 0) {
        xpipars_aux.set_size(m*numY,p1+ptrt); xpipars_aux.zeros();
      }
      if (p2 > 0) {
			  xmupars_aux.set_size(m*numY,p2); xmupars_aux.zeros();
			  xsigpars_aux.set_size(m*numY,p2); xsigpars_aux.zeros();
			}

			for(int w = 0; w < m*numY; w++) {
			
        if (p1 + ptrt > 0) {
          for(int ww = 0; ww < ptrt+p1; ww++) {
            xpipars_aux(w,ww) = R::rbeta(a0,b0);
          }
        }
			
        if (p2 > 0) {
          for(int ww = 0; ww < p2; ww++) {
            xsigpars_aux(w,ww) = rinvchisq(nu0,tau0);
            xmupars_aux(w,ww) = R::rnorm(mu0,xsigpars_aux(w,ww)/sqrt(c0));
          }
        }
			}
			
			//begin probs for existing Y clusters, but new X clusters
			int count2 = 0;
			//Rprintf("Calculate probabilities for new X clusters\n");
			for(int k = 0; k < numY; k++) {
				//FILL IN PROBS FOR NEW X CLUSTERS IN EXISTING Y CLUSTERS
				vdummy1 = find(Sy == (k+1));
				njwoi = vdummy1.size();
				
				likeregy = 1;
				
				for(int kk = 0; kk < ntp(j); kk++) {
					likeregy = likeregy * R::dnorm(yvals(kk),dot(xvals.row(kk),betaY.row(k))+u(j),sqrt(sig2(k)),0);
				} //end for kk

				
				for(int kk = 0; kk < m; kk++) {
					prodx=1;prodx2=1;
					//likelihood for binary covariates
          if (p1 + ptrt > 0) {
            for(int q = 0; q < ptrt+p1; q++) {
              prodx = prodx*R::dbinom(Xonly(j,q), 1, xpipars_aux(count2, q), 0);
            } //end for q
          }
					
					//likelihood for continuous covariates
          if (p2 > 0) {
            for(int q = 0; q < p2; q++) {
              prodx2 = prodx2*R::dnorm(Xonly(j,ptrt+p1+q), xmupars_aux(count2,q), sqrt(xsigpars_aux(count2,q)), 0 );
            } //end for q
					}
					
					probs(count) = ((njwoi*alphapsi/m)/(njwoi+alphapsi))*likeregy*prodx*prodx2;
					if(probs(count)<0) {
						Rprintf("NEGATIVE PROBABILITY ALERT\n");
						Rprintf("prodx = %.8f\n",prodx);
						Rprintf("vals: %.1f, %.1f, %.1f\n",Xonly(j,0),Xonly(j,1),Xonly(j,2));
						Rprintf("prodx2 = %.8f\n",prodx2);
						Rprintf("vals: %.2f, %.2f\n",Xonly(j,ptrt+p1),Xonly(j,ptrt+p1+1));
						Rprintf("likeregy = %.8f\n",likeregy);
					}
					//Rprintf("count + prob: %.2f, %d\n",count,probs(count));
					count++;count2++;
					
				}
				
			}
			
			
			
			
			
			
			
			// COMPLETELY NEW CLUSTERS
			
			//set auxiliary parameters
			betaY_aux.set_size(m,nbeta); betaY_aux.zeros();
			sig2_aux.set_size(m); sig2_aux.zeros();
			
		  if (p1 + ptrt > 0) {	
			  Y_xpipars.set_size(m,p1+ptrt); Y_xpipars.zeros();
			}
      if (p2 > 0) {
        Y_xmupars.set_size(m,p2); Y_xmupars.zeros();
			  Y_xsigpars.set_size(m,p2); Y_xsigpars.zeros();
      }

			for(int w = 0; w < m; w++) {
				
				sig2_aux(w) = 1/R::rgamma(beta_a0,beta_b0);
				//betaY_aux.row(w) = trans(mvrnorm(beta0,sig2_aux(w)*prec0));
				betaY_aux.row(w) = trans(mvrnorm(beta0,sig2_aux(w)*prec0.i()));
				
				
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
					likeregy = likeregy * R::dnorm(yvals(kk),dot(xvals.row(kk),betaY_aux.row(k))+u(j),sqrt(sig2_aux(k)),0);
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
				
				probs(count) = (alphatheta/m)*prodx*prodx2*likeregy;
				//Rprintf("count + prob: %.2f, %d\n",count,probs(count));
				if(probs(count)<0) {
					Rprintf("NEGATIVE PROBABILITY ALERT\n");
					Rprintf("prodx = %.8f\n",prodx);
					Rprintf("vals: %.1f, %.1f, %.1f\n",Xonly(j,0),Xonly(j,1),Xonly(j,2));
					Rprintf("prodx2 = %.8f\n",prodx2);
					Rprintf("vals: %.2f, %.2f\n",Xonly(j,ptrt+p1),Xonly(j,ptrt+p1+1));
					Rprintf("likeregy = %.8f\n",likeregy);
				}
				count++;
			}
			
			// Rcout << "probs: " << probs.t() << std::endl;
			
			
			//USE MULTINOMIAL DISTRIBUTION TO CHOOSE NEW CLUSTER
			newCluster = rmultinomF(probs);
			//Rprintf("Chosen cluster: %d\n",newCluster);
			
			probs.zeros();
			
			
			
			if(newCluster<=numTotalCluster) 
			{
				//Rprintf("Chose existing cluster\n");
				Sy.insert_rows(j,1);
				Sy(j) = uniqueS(newCluster-1,0);
				Sx.insert_rows(j,1);
				Sx(j) = uniqueS(newCluster-1,1);
			}
			else if(newCluster>numTotalCluster+m*numY)  //new Y and X cluster
			{
				//Rprintf("Chose new cluster\n");
				Sx.insert_rows(j,1);
				Sx(j) = 1;
				
				Sy.insert_rows(j,1);
				Sy(j) = Sy.max()+1;
				//Rprintf("new Y cluster number: %d\n",Sy(j));
				
				uniqueS.insert_rows(numTotalCluster,1);
				uniqueS(numTotalCluster,0) = Sy(j);
				uniqueS(numTotalCluster,1) = Sx(j);
				betaY.insert_rows(numY,1);
				betaY.row(numY) = betaY_aux.row(newCluster-numTotalCluster-m*numY-1);
				sig2.insert_rows(numY,1);
				sig2.row(numY) = sig2_aux.row(newCluster-numTotalCluster-m*numY-1);
				if (p1 + ptrt > 0) {
          xpipars.insert_rows(numTotalCluster,1);
				  xpipars.row(numTotalCluster) = Y_xpipars.row(newCluster-numTotalCluster-m*numY-1);
				}
        if (p2 > 0) {
          xmupars.insert_rows(numTotalCluster,1);
          xmupars.row(numTotalCluster) = Y_xmupars.row(newCluster-numTotalCluster-m*numY-1);
          xsigpars.insert_rows(numTotalCluster,1);
          xsigpars.row(numTotalCluster) = Y_xsigpars.row(newCluster-numTotalCluster-m*numY-1);
        }
			}
			else  //new X cluster within Y cluster
			{
				//Rprintf("Chose new X cluster\n");
				Sy.insert_rows(j,1);
				Sy(j) = ceil((1.0*newCluster-numTotalCluster)/m);
				vdummy1 = find(uniqueS.col(0) == Sy(j));
				numXj = vdummy1.size();
				Sx.insert_rows(j,1);
				Sx(j) = numXj + 1;
				vdummy1 = find(uniqueS.col(0) <= Sy(j));
				numXj = vdummy1.size();
				uniqueS.insert_rows(numXj,1);
				uniqueS(numXj,0) = Sy(j);
				uniqueS(numXj,1) = Sx(j);
				if (p1 + ptrt > 0) {
          xpipars.insert_rows(numXj,1);
				  xpipars.row(numXj) = xpipars_aux.row(newCluster-numTotalCluster-1);
				}
        if (p2 > 0) {
          xmupars.insert_rows(numXj,1);
				  xmupars.row(numXj) = xmupars_aux.row(newCluster-numTotalCluster-1);
				  xsigpars.insert_rows(numXj,1);
				  xsigpars.row(numXj) = xsigpars_aux.row(newCluster-numTotalCluster-1);
        }
			}
			
			//Rprintf("Y cluster: %d\n",Sy(j));
			//Rprintf("X cluster: %d\n",Sx(j));
			
		} // end cluster loop (j)
		
		
		//Rprintf("Loop succesfully completed.\n");
		return List::create(_["uniqueS"] = uniqueS,
												 _["Sy"] = Sy,
												 _["Sx"] = Sx,
												 _["xpipars"] = xpipars,
												 _["xmupars"] = xmupars,
												 _["xsigpars"] = xsigpars,
												 _["betaY"] = betaY,
												 _["sig2"] = sig2);
		

	}


}  // end function

