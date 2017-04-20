#ifndef __UTIL_H__
#define __UTIL_H__

#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

int rmultinomF(vec const&p); //random multinomial
vec mvrnorm(vec mu, mat sigma); //multivariate normal draw
double rinvchisq(double df, double scale); //random inverse chi sq
double dmvn(vec x, vec mu, mat sig, int logt); //evaluate MVN density
// double newalphatheta(int numYclus, double prioralp, double alp_a0, double alp_b0,int nobs);
// double newalphapsi(int numYclus, int numTotalCluster, double alpha, double alp_a0, double alp_b0);

#endif
