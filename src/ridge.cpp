#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat woodbury_rreg(arma::mat X, arma::mat Y, double lambda){
  arma::mat inverse = arma::inv(arma::eye(X.n_rows, X.n_rows) + (1/lambda) * X * X.t());
  return (1/lambda)*X.t()*Y - (1/(lambda*lambda))*X.t()*inverse*X*X.t()*Y;
}

// [[Rcpp::export]]
arma::mat svd_rreg(arma::mat X, arma::mat Y, double lambda){
  arma::mat U, D, V;
  arma::vec s;
  arma::svd(U,s,V,X);

  D = arma::diagmat(s);

  //int bound = D.n_rows;
  V=V.head_cols(X.n_rows);

  arma::mat R = U*D;

  arma::mat s_inverse = arma::inv(R.t()*R+lambda*arma::eye(X.n_rows, X.n_rows));
  return V*s_inverse*R.t()*Y;
}
