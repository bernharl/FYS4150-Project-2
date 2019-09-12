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
  double t = -tau;
  if(tau > 0)
  {
    t += sqrt(1. + tau * tau);
  }
  else
  {
    t -= sqrt(1. + tau * tau);
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

int main()
{
  //arma::mat A = arma::mat(n, n, arma::fill::zeros);
  arma::mat zero_diag_matrix = arma::ones <arma::mat> (n, n);
  zero_diag_matrix.diag(0).fill(0);
  arma::mat A = arma::zeros <arma::mat> (n, n);
  A.diag(0).fill(d);
  A.diag(1).fill(a);
  A.diag(-1).fill(a);

  while((A.*zero_diag_matrix)*(A.*zero_diag_matrix)).max() > eps)
  {

  }
  return 0;
}
