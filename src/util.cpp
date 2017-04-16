#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

int rmultinomF(vec const& p) {
  vec csp = cumsum(p/sum(p));
  double rnd = runif(1)[0];
  int res = 0;
  int psize = p.size();

  for(int i = 0; i < psize; i++) {
    if(rnd>csp(i)) res = res+1;
  }

  return(res+1); //think need to return res not res+1 -- edit:: maybe not
}

// taken from: http://gallery.rcpp.org/articles/simulate-multivariate-normal/
vec mvrnorm(vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(1, ncols);
  return mu + trans(Y * chol(sigma));
}

double rinvchisq(double df, double scale) {
  double res = as<double>(rchisq(1,df));
  return( (df*scale)/res );
}

//adapted from dmvnorm function from mvtnorm package in R
double dmvn(vec x, vec mu, mat sig, int logt) {
  mat dec = chol(sig);
  vec tmp = solve(dec,x-mu); //check
  double rss = sum(tmp%tmp);
  
  double logretval = -sum( log (dec.diag() ) ) - 0.5 * sig.n_cols * log (2 * PI) - 0.5 * rss;
  if(logt==0)
    logretval = exp(logretval);
  return logretval;
}

/*
double newalphatheta(int numYclus, double prioralp, double alp_a0, 
		     double alp_b0, int nobs) {
  double eta = R::rbeta(prioralp+1,nobs);
  double pieta = (numYclus/(nobs*(1-log(eta))))/(1+numYclus/(nobs*(1-log(eta))));
  int whichmix = R::rbinom(1,pieta);
  double newalp = whichmix*R::rgamma(alp_a0+numYclus,alp_b0-log(eta)) +
    (1-whichmix)*R::rgamma(alp_a0+numYclus-1,alp_b0-log(eta));
  return newalp;
}
*/

/*

newalppsi<-function(alpha){
  sortuniquey<-sort(unique(s[,1]))
  ss<-numeric(length(sortuniquey))
  for(j in 1:length(sortuniquey)){
    ss[j]<-sum(s[,1]==sortuniquey[j])
  }
  likecur<-dgamma(alpha,alpa0,alpb0)*(alpha^(nrow(unique(s))-k))*prod((alpha+ss)*beta(alpha+1,ss))
  proposed<-rgamma(1,2,1)
    likeprop<-dgamma(proposed,alpa0,alpb0)*(proposed^(nrow(unique(s))-k))*prod((proposed+ss)*beta(proposed+1,ss))
  rat<-(likeprop)/(likecur)
  if(runif(1,0,1)<rat){alpha<-proposed}
  return(alpha)
}
*/
/*
double newalphapsi(int numYclus, int numTotalCluster, double alpha, double alp_a0, double alp_b0, mat uniqueS) {
  ivec ss(numYclus);
  vec ss2(numYclus);
  vec ss3(numYclus);
  for(int ii = 0; ii < numYclus; ii++) {
    uvec dum_ind = find(uniqueS.col(0) == (ii+1));
    int numdum = dum_ind.size();
    ss(ii) = numdum;
    ss2(ii) = (alpha+ss(ii))*R::beta(alpha+1,ss(ii));
  }
  
  double likecur = R::dgamma(alpha,alp_a0,alp_b0,0)*(alpha^(numTotalCluster-numYclus))*prod(ss2);
  double proposed = R::rgamma(2,1); //vars?
  for(int ii = 0; ii < numYclus; ii++) {
    ss3(ii) = (proposed+ss(ii))*R::beta(proposed+1,ss(ii));
  }
  double likeprop = R::dgamma(proposed,alp_a0,alp_b0,0)*(proposed^(numTotalCluster-numYclus))*prod(ss3);
  double rat = likeprop/likecur;
  return alpha;
    
}
*/
