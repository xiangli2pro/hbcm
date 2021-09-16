#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_qc(int n, int p, int centers, 
                      NumericVector ppi, 
                      NumericVector hsigma, NumericVector hlambda, 
                      NumericMatrix alpha_mu, NumericMatrix X, 
                      List alpha_mu_inter){
  
  NumericMatrix qc(centers, p);
  
  for(int j = 0; j < p; ++j){
    for(int k = 0; k < centers; ++k){
      
      double p1 = log(ppi[k]);
      double p2 = -n / 2 * log(std::pow(hsigma[j], 2));
      double p3 = 0.0;
      
      for(int i = 0; i < n; ++i){
        
        NumericMatrix alpha_mu_inter_i = alpha_mu_inter[i];
        
        p3 += (std::pow(X(i,j), 2) - 2 * hlambda[j] * X(i,j) * alpha_mu(i,k) +
          std::pow(hlambda[j], 2) * alpha_mu_inter_i(k, k));
      }
      
      p3 = -1 / (std::pow(hsigma[j], 2)) * p3 * 1 / 2;
      qc(k,j) = p1 + p2 + p3;
    }
  }
  
  return qc;
}

