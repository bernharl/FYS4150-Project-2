#include "catch.hpp"
#include "main.h"
#include "armadillo"



TEST_CASE("Comparing eigenvalues with armadillo library")
{ 
  n = 5;
  eps = 1e-8;

  
  rhoN = 5;
  h = (rhoN - rho0) / ((double) n + 1);
  d = 2. / (h * h);
  a = -1. / (h * h);

  arma::mat A = arma::zeros <arma::mat> (n, n);
  arma::mat E = arma::eye <arma::mat> (n, n);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);

  //Finding eigenvalues with armadillo
  arma::cx_vec eig_val;
  arma::cx_mat eig_vec;
  arma::eig_gen(eig_val, eig_vec, A);
  
  //Finding eigenvalues with Jacobi algorithm
  Jacobi_Algorithm(A, E, n, h, eps);
  arma::vec calculated_eig_vals = A.diag();
  REQUIRE(arma::norm(arma::sort(A.diag()) - arma::sort(arma::real(eig_val))) == Approx(0).epsilon(eps));
}