#include <Rcpp.h>
#include <vector>
#include <unordered_map>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericMatrix cpp_beta_mat(IntegerVector population, double beta, NumericMatrix distance, 
                           double nu = 0, double mu = 0, double epsilon = 0, double iota = 0, 
                           double connection_strength = 1, bool mobility = false) {
  
  int n_patch = population.size();
  NumericMatrix beta_matrix(n_patch, n_patch);
  
  // Compute the beta_matrix
  NumericVector density_dependent(n_patch);
  
  if (iota) { 
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        if (i != j) {
          density_dependent[i] += pow(population[j], nu) / pow(distance(i,j), epsilon);
        }
      }
    } 
  }
  
  
  for (int i = 0; i < n_patch; ++i) {
    for (int j = 0; j < n_patch; ++j) {
      if (j != i) {
        if (mobility) {
          beta_matrix(i,j) = beta * connection_strength * pow(population[i], nu) * pow(population[j], mu) * pow(distance(i,j), epsilon) / pow(density_dependent[i], iota); 
        } else {
          beta_matrix(i,j) = beta * connection_strength * pow(population[i], nu) * pow(population[j], mu) / (pow(distance(i,j), epsilon) * pow(density_dependent[i], iota)); 
        }
      } else {
        beta_matrix(i,j) = beta;
      }
    }
  }
  

  return beta_matrix;
  
}

//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericMatrix cpp_gamma_mat(IntegerVector population, NumericMatrix distance, 
                           double nu = 0, double mu = 0, double epsilon = 0, double iota = 0, 
                           double connection_strength = 1, bool mobility = false) {
  
  int n_patch = population.size();
  NumericMatrix beta_matrix(n_patch, n_patch);
  
  // Compute the beta_matrix
  NumericVector density_dependent(n_patch);
  
  if (iota) { 
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        if (i != j) {
          density_dependent[i] += pow(population[j], nu) / pow(distance(i,j), epsilon);
        }
      }
    } 
  }
  
  
  for (int i = 0; i < n_patch; ++i) {
    for (int j = 0; j < n_patch; ++j) {
      if (j != i) {
        if (mobility) {
          beta_matrix(i,j) = connection_strength * pow(population[i], nu) * pow(population[j], mu) * pow(distance(i,j), epsilon) / pow(density_dependent[i], iota);
        } else {
          beta_matrix(i,j) = connection_strength * pow(population[i], nu) * pow(population[j], mu) / (pow(distance(i,j), epsilon) * pow(density_dependent[i], iota));
        }
      }
    }
  }
  
  
  return beta_matrix;
  
}


//////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericMatrix cpp_beta_adj(IntegerVector population, NumericVector beta, NumericMatrix distance, 
                           double nu = 0, double mu = 0, double epsilon = 0, double iota = 0,
                           double connection_strength = 1, bool mobility = false) {
  
  int n_patch = population.size();
  NumericMatrix beta_matrix(n_patch, n_patch);
  NumericVector density_dependent(n_patch);
  
 if (iota) {
   for (int i = 0; i < n_patch; ++i) {
     for (int j = 0; j < n_patch; ++j) {
       if (i != j && distance(i,j) != 0) {
         density_dependent[i] += pow(population[j], nu) / pow(distance(i,j), epsilon);
       }
     }
   }
 }


 for (int i = 0; i < n_patch; ++i) {
   for (int j = 0;j < n_patch; ++j) {
     if (j != i && distance(i,j) != 0) {
       if (mobility) {
         beta_matrix(i,j) = beta[j] * connection_strength * pow(population[i], nu) * pow(population[j], mu) * pow(distance(i,j), epsilon) / pow(density_dependent[i], iota);
       } else {
         beta_matrix(i,j) = beta[j] * connection_strength * pow(population[i], nu) * pow(population[j], mu) / (pow(distance(i,j), epsilon) * pow(density_dependent[i], iota)); 
       }
     } else if (j == i){
       beta_matrix(i,j) = beta[i];
     }
   }
 }


 return beta_matrix;
    
}
  
