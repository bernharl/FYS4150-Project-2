#include "catch.hpp"
#include "jacobi_algorithm.h"
#include "armadillo"



TEST_CASE("Comparing eigenvalues with armadillo library")
{ 
  int n = 5;
  double eps = 1e-8;

  double rho0 = 0.0;
  double rhoN = 5.0;
  double h = (rhoN - rho0) / ((double) n + 1);
  double d = 2. / (h * h);
  double a = -1. / (h * h);

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

TEST_CASE("Verifies that max_offdiag method choosesthe correct maximum value and checks indices")
{
  
  double t = -1e3; //Number with larger fabs than 0
  int N = 5;
  int k, l;
  arma::mat T = arma::zeros <arma::mat> (N, N); //Making matrix of zeros
  T(2, 1)     = t; // Setting one element to high value
  double max_val = max_offdiag(T, k, l, N);
  REQUIRE(k == 2);
  REQUIRE(l == 1);
  REQUIRE(max_val == fabs(t));
}

TEST_CASE("Verifies that the offdiagonal elements is zero after calling Jacobi_Algorithm method.")
{
  int n = 5;
  double eps = 1e-8;

  double rho0 = 0.0;
  double rhoN = 5.0;
  double h = (rhoN - rho0) / ((double) n + 1);
  double d = 2. / (h * h);
  double a = -1. / (h * h);

  arma::mat A = arma::zeros <arma::mat> (n, n);
  arma::mat E = arma::eye <arma::mat> (n, n);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);
  Jacobi_Algorithm(A, E, n, h, eps);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i != j){
        REQUIRE(fabs(A(i, j) * A(i, j)) <= eps);
      }
    }
  }
}
TEST_CASE("Verifies that the columns in A are orthogonal after calling the Jacobi_Algorithm method.")
{
  int n = 3;
  double eps = 1e-8;

  double rho0 = 0.0;
  double rhoN = 5.0;
  double h = (rhoN - rho0) / ((double) n + 1);
  double d = 2. / (h * h);
  double a = -1. / (h * h);
  arma::mat A = arma::zeros <arma::mat> (n, n);
  arma::mat E = arma::eye <arma::mat> (n, n);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);
  Jacobi_Algorithm(A, E, n, h, eps);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i != j){
        REQUIRE(arma::dot(E.col(i), E.col(j)) <= eps);
      }
    }
  }
}