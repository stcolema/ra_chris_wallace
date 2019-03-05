# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Compares how similar two points are with regards to their clustering across 
// all iterations.
// Works for unsupervised methods (i.e. allows label flipping)
double point_similarity(arma::uword point, 
                        arma::uword comparison_point,
                        arma::umat cluster_record,
                        arma::uword num_iter) {
  double out = 0.0;
  
  arma::uvec ind_1(num_iter);
  arma::uvec ind_2(num_iter);
  
  ind_1 = arma::trans(cluster_record.row(point));
  ind_2 = arma::trans(cluster_record.row(comparison_point));
  
  arma::umat out_1(num_iter, 1);
  
  // Compare vector of allocations element-wise
  out_1.col(0) = (ind_1 == ind_2);
  
  // Similarity is the sum of the above divided by the number of entries
  // Convert the sum to a double as otherwise is integer divison and does not 
  // work
  out = (double)arma::sum(out_1.col(0)) / (double)num_iter;
  
  return out;
}

// Constructs a similarity matrix comparing all points clustering across the 
// iterations
// [[Rcpp::export]]
arma::mat similarity_mat(arma::umat cluster_record){
  arma::uword sample_size = cluster_record.n_rows;
  arma::uword num_iter = cluster_record.n_cols;
  arma::mat out(sample_size, sample_size);
  out.zeros();
  
  // if not doing diagonal, restrict to sample size - 1
  for (arma::uword point = 0; point < sample_size; point++){ 
    for (arma::uword comparison_point = point; 
         comparison_point < sample_size;
         comparison_point++){
      out(point, comparison_point) = point_similarity(point, 
          comparison_point,
          cluster_record,
          num_iter);
      out(comparison_point, point) = out(point, comparison_point);
    }
  }
  return out;
}
