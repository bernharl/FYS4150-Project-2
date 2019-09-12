#include <armadillo>
#include <iostream>

#include "keyword_parameters.h"
using namespace std;


inline double tau(double a_ll, double a_kk, double a_kl)
{
  return (a_ll - a_kk) / (2. * a_kl);
}


inline double t(double tau)
{
  double t = 1.0;
  if(tau > 0)
  {
    t /= (tau + sqrt(1. + tau * tau));
  }
  else
  {
    t /= - (-tau + sqrt(1. + tau * tau));
  }
  return t;
}

inline double c(double t)
{
  return 1. / sqrt(1. + t * t);
}

inline double s(double c, double t)
{
  return t * c;
}

double max_offdiag(const arma::mat &A, int &k, int &l, int n)
{ 
  double max_val = 0.0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (fabs(A(i, j)) > max_val && i != j){
        max_val = fabs(A(i, j));
        k = i;
        l = j;
      }
    }
  }
  cout << "jeff" << max_val << endl;
  return max_val;
}

int main()
{
  //arma::mat A = arma::mat(n, n, arma::fill::zeros);
  arma::mat zero_diag_matrix = arma::ones <arma::mat> (n, n);
  zero_diag_matrix.diag(0).fill(0);
  arma::mat A = arma::zeros <arma::mat> (n, n);
  arma::mat B = arma::zeros <arma::mat> (n, n);
  A.diag(0).fill(d);
  A.diag(1).fill(a);
  A.diag(-1).fill(a);
  int k, l;

  double max_val = max_offdiag(A, k, l, n);
  //cout << max_val << endl;
  double a_ll, a_kk, a_ik, a_il, a_kl;//, a_lk;
  cout << max_val << " " << k << " " << l << endl;

  a_kl = A(k, l);
  //a_lk = A(l, k);
  //cout << A << endl;
  while(max_val * max_val > eps)
  { 
    
    
    double tau_val = tau(a_ll, a_kk, a_kl);
    double t_val = t(tau_val);
    double c_val = c(t_val);
    double s_val = s(c_val, t_val);
    
    B = A;
    a_ll = A(l, l);
    a_kk = A(k, k);
  
    for (int i = 0; i < n; i++)
    { 
      a_ik = A(i, k);
      a_il = A(i, l);
      if (i != k && i != l)
      {
        B(i, i) = A(i, i);
        B(i, k) = a_ik * c_val - a_il * s_val;
        B(i, l) = a_il * c_val + a_ik * s_val;
      }
      B(k, k) = a_kk * c_val * c_val - 2 * a_kl * c_val * s_val + a_ll * s_val * s_val;
      B(l, l) = a_ll * c_val * c_val + 2 * a_kl * c_val * s_val + a_kk * s_val * s_val;
      B(l, k) = B(k, l) = (a_kk - a_ll) * c_val * s_val + a_kl * (c_val * c_val - s_val * s_val);
      A = B;
    }
    //cout << A << endl;
    max_val = max_offdiag(A, k, l, n);  
  }
  cout << A << endl;
  return 0;
}
