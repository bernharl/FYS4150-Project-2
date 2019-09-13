#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "main.h"
//#include "keyword_parameters.h"
using namespace std;


double max_offdiag(const arma::mat &A, int &k, int &l, int n)
{
  double max_val = 0.0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i != j && fabs(A(i, j)) > max_val){
        max_val = fabs(A(i, j));
        k = i;
        l = j;
      }
    }
  }
  return max_val;
}

void Jacobi_Algorithm(arma::mat &A)
{

  arma::mat B;
  int k, l;

  double max_val = max_offdiag(A, k, l, n);
  double a_ll, a_kk, a_ik, a_il, a_kl;
  double t_val, tau, c, s; 

  while(max_val * max_val > eps)
  {
    a_kl = A(k, l);
    a_ll = A(l, l);
    a_kk = A(k, k);

    tau = (a_ll - a_kk) / (2. * a_kl);
    if(tau > 0)
    {
      t_val = 1.0 / (tau + sqrt(1. + tau * tau));
    }
    else
    {
      t_val = -1.0 / (-tau + sqrt(1. + tau * tau));
    }
    c = 1. / sqrt(1. + t_val * t_val);
    s = t_val * c; 

    A(k, k) = a_kk * c * c - 2.0 * a_kl * c * s + a_ll * s * s;
    A(l, l) = a_ll * c * c + 2.0 * a_kl * c * s + a_kk * s * s;
    A(l, k) = A(k, l) = 0;

    for (int i = 0; i < n; i++)
    {
      a_ik = A(i, k);
      a_il = A(i, l);
      if (i != k && i != l)
      {
        A(i, k) = A(k, i) = a_ik * c - a_il * s;
        A(i, l) = A(l, i) = a_il * c + a_ik * s;
      }

    }
    max_val = max_offdiag(A, k, l, n);
  }
    //cout << A << endl;
}


void TEST_OFFMAX()
{ 
  double t = -1e3; //Number with larger fabs than 0
  int n = 5;
  int k, l;
  arma::mat T = arma::zeros <arma::mat> (5, 5); //Making matrix of zeros
  T(2, 1)     = t; // Setting one element to high value
  double max_val = max_offdiag(T, k, l, n);
  assert(max_val == fabs(t) && k == 2 && l == 1);
}

void TEST_OFFDIAG_IS_ZERO()
{
  arma::mat A = arma::zeros <arma::mat> (5, 5);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);
  Jacobi_Algorithm(A);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i != j){
        assert(fabs(A(i, j) * A(i, j)) <= eps);
      }
    }
  }

}

int main()
{ 

  TEST_OFFMAX();
  TEST_OFFDIAG_IS_ZERO();
  arma::mat A = arma::zeros <arma::mat> (n, n);
  A.diag(0).fill(d);
  A.diag(1).fill(a);
  A.diag(-1).fill(a);
  return 0;
}

