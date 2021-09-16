#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_qalpha_mu_inter(int n, int p, int centers, 
                          NumericMatrix alpha_mu, 
                          NumericMatrix alpha_cov){
  
  List alpha_mu_inter(n);
  
  for(int i = 0; i < n; ++i){
    
    NumericMatrix mat(centers, centers);
    
    for(int j = 0; j < centers; ++j){
      for(int k = 0; k < centers; ++k){
        
        mat(j,k) = alpha_cov(j, k) + alpha_mu(i, j) * alpha_mu(i, k);
      }
    }
    
    alpha_mu_inter[i] = mat;
  }
  
  return alpha_mu_inter;
}