#include "catch.hpp"
#include "main.h"
//#include "keyword_parameters.hpp"
#include "armadillo"

TEST_CASE("Comparing eigenvalues with armadillo library")
{ 
  arma::mat A = arma::zeros <arma::mat> (5, 5);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);

  //Finding eigenvalues with armadillo
  arma::cx_vec eig_val;
  arma::cx_mat eig_vec;
  arma::eig_gen(eig_val, eig_vec, A);
  
  //Finding eigenvalues with Jacobi algorithm
  Jacobi_Algorithm(A);
  arma::vec calculated_eig_vals = A.diag();
  REQUIRE(norm(sort(A.diag()) - sort(real(eig_val))) <= eps);
}
