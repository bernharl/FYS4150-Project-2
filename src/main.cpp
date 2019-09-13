#include <armadillo>
#include <iostream>

#include "keyword_parameters.h"
using namespace std;


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
  return max_val;
}


int main()
{
 
  arma::mat A = arma::zeros <arma::mat> (n, n);
  arma::mat B;
  A.diag(0).fill(d);
  A.diag(1).fill(a);
  A.diag(-1).fill(a);
  arma::cx_vec eig_val;
  arma::cx_mat eig_vec;
  arma::eig_gen(eig_val, eig_vec, A);
  cout << eig_val << endl;
  int k, l;

  double max_val = max_offdiag(A, k, l, n);
  cout << k << " " << l << endl;
  double a_ll, a_kk, a_ik, a_il, a_kl;//, a_lk;
  //cout << max_val << " " << k << " " << l << endl;

  double t_val, tau_val, c_val, s_val; 
  while(max_val * max_val > eps)
  {
    a_kl = A(k, l);
    a_ll = A(l, l);
    a_kk = A(k, k);

    tau_val = (a_ll - a_kk) / (2. * a_kl); // tau(a_ll, a_kk, a_kl);
    if(tau_val > 0)
    {
      t_val = 1.0 / (tau_val + sqrt(1. + tau_val * tau_val));
    }
    else
    {
      t_val = -1.0 / (-tau_val + sqrt(1. + tau_val * tau_val));
    }
    c_val = 1. / sqrt(1. + t_val * t_val); //c(t_val);
    s_val = t_val * c_val; // s(c_val, t_val);

    // A = A;
    //cout << "hei " << a_kk << " " << a_ll << endl;
    //cout << "hei2 " << a_kk << " " << a << endl;
    A(k, k) = a_kk * c_val * c_val - 2.0 * a_kl * c_val * s_val + a_ll * s_val * s_val;
    A(l, l) = a_ll * c_val * c_val + 2.0 * a_kl * c_val * s_val + a_kk * s_val * s_val;
    A(l, k) = A(k, l) = 0;

    for (int i = 0; i < n; i++)
    {
      a_ik = A(i, k);
      a_il = A(i, l);
      if (i != k && i != l)
      {
        A(i, k) = A(k, i) = a_ik * c_val - a_il * s_val;
        A(i, l) = A(l, i) = a_il * c_val + a_ik * s_val;
      }

    }

    // A = A;
    //cout << A << endl;
    max_val = max_offdiag(A, k, l, n);
    cout << "Max value " << max_val << endl;
  }
    cout << A << endl;
  return 0;
}

