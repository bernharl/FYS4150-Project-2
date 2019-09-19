#include <armadillo>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <iomanip>
#include <ctime>
#include <cmath>

#include "keyword_parameters.h"
using namespace std;


double max_offdiag(const arma::mat &A, int &k, int &l, int n)
/*
Finds the maximum value in an nxn matrix A
excluding the diagonal.

Parameters
------------
A: arma::mat
  An armadillo object nxn matrix
n: int
 The dimensions of A
k,l: int
  The indices of the maximum value
 */
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

void Jacobi_Algorithm(arma::mat &A, arma::mat &E, int n, double h)
/*
Implements the jacobi algorithm for finding the
eigenvalues on an nxn matrix

Parameters
------------
A: arma::mat
  An armadillo object nxn matrix
E: arma::mat
  An armadillo object nxn matrix.
  Serves as an orthogonal basis of n dimensional space
n: int
 The dimensions of A
h: double
  The step size
 */
{

  int k, l;
  double max_val = max_offdiag(A, k, l, n);
  double a_ll, a_kk, a_ik, a_il, a_kl, e_ik, e_il;
  double t_val, tau, c, s; 
  double iterator = 0;

  while(max_val * max_val > eps && iterator <= tot_iterations)
  { 
    iterator ++;
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
      e_ik = E(i, k);
      e_il = E(i, l);
      E(i, k) = c * e_ik - s * e_il;
      E(i, l) = c * e_il + s * e_ik;
    }
    max_val = max_offdiag(A, k, l, n);
  }
}

void Harmonic_Potential(arma::vec &V, double rho0, double h, int n)
/*
Computes the harmonic oscillator potential a particle

Parameters
-----------
V: arma::vec
  An armadillo vector containing the potential
rho0: double
 The minumum value for the potential
h: double
 The step length
n: int
 The length of the vector V
 */
{ 
  for (int i = 0; i < n; i++)
  {
    V(i) = (rho0 + (i + 1) * h) * (rho0 + (i + 1) * h);
  }
}

void coulomb_potential(arma::vec &V, double rho0, double h, int n, double omega_r)
/*
Computes the coulomb potential for a particle

Parameters
-----------
V: arma::vec
  An armadillo vector containing the potential
rho0: double
  The minumum value for the potential
h: double
  The step length
n: int 
  The length of the vector V
omega_r: double
 The oscillator frequency reflecting
         the strength of the potential
 */
{ 
  double rho;
  for (int i = 0; i < n; i++)
  { 
    rho = rho0 + (i + 1) * h;
    V(i) = omega_r * omega_r * rho * rho + 1 / rho;
  }
}




int main()
{ 
  /* 
  TEST_OFFMAX();
  TEST_JACOBI_ALGORITHM();
  TEST_OFFDIAG_IS_ZERO();
  TEST_ORTHOGONALITY();

  arma::mat A = arma::zeros <arma::mat> (n, n);
  arma::mat E = arma::eye <arma::mat> (n, n);
  A.diag(0).fill(d);
  A.diag(1).fill(a);
  A.diag(-1).fill(a);

  Jacobi_Algorithm(A, E, n);
  
  double h;
  double rhoN;
  clock_t t_start = clock(); // Initializing timer
  arma::vec V = arma::zeros <arma::vec> (n);
  Harmonic_Potential(V, rho0, h, n);
  arma::mat P = arma::zeros <arma::mat> (n, n);
  arma::mat E = arma::eye <arma::mat> (n, n);
  P.diag(0) += d + V;
  P.diag(1).fill(a);
  P.diag(-1).fill(a);

  clock_t t_start = clock(); // Initializing timer
  Jacobi_Algorithm(P, E, n);
  clock_t t_end = clock(); // End timer
  double CPU_time = (t_end - t_start) / CLOCKS_PER_SEC; // Calculating CPU time [ms]
  arma:: vec diags = arma::sort(P.diag(0));
  cout << diags(0) << " " << diags(1) << " " 
       << diags(2) << " " << diags(3) << endl;
  cout << "Run time: " << CPU_time << " s " << endl;

  */

  double h;
  double rhoN;
  double d;
  double a;
  int length = 5;
  arma:: vec rho_max = arma::linspace<arma::vec>(1 ,20, length);

  double vars[5][length];

  clock_t t_start = clock(); // Initializing timer

  double progress = 0;
  //#pragma omp parallel for
  for (int i=0; i<length; i++) 
  {
    arma::mat P = arma::zeros <arma::mat> (n, n);
    rhoN = rho_max(i);
    h = (rhoN - rho0) / ((double) n + 1);
    d = 2. / (h * h);
    a = -1. / (h * h);

    arma::vec V = arma::zeros <arma::vec> (n);
    Harmonic_Potential(V, rho0, h, n);
    arma::mat E = arma::eye <arma::mat> (n, n);
    P.diag(0) += d + V;
    P.diag(1).fill(a);
    P.diag(-1).fill(a);
    Jacobi_Algorithm(P, E, n, h);
    arma:: vec diags = arma::sort(P.diag(0));

    vars[0][i] = rhoN;
    vars[1][i] = diags(0);
    vars[2][i] = diags(1);
    vars[3][i] = diags(2);
    vars[4][i] = diags(3);
    progress += 1;
    cout << "iterations left: "<< progress << "/" << length << endl;
  }
  clock_t t_end = clock(); // End timer
  double CPU_time = (t_end - t_start) / CLOCKS_PER_SEC;
  cout << "Run time: " << CPU_time << " s " << endl;

  ofstream outfile;
  outfile.open("eigendata.dat");
  outfile << setw(20) << "rhoN"
          << setw(20) << "l1"
          << setw(20) << "l2"
          << setw(20) << "l3"
          << setw(20) << "l4" << endl;


  for (int i=0; i<length; i++) 
  {
    outfile << setprecision(7) << setw(20) << vars[0][i]
            << setprecision(7) << setw(20) << vars[1][i]
            << setprecision(7) << setw(20) << vars[2][i]
            << setprecision(7) << setw(20) << vars[3][i]
            << setprecision(7) << setw(20) << vars[4][i] << endl;

  }
  outfile.close();

  return 0;
}
double rhoN = 5.4;
double h = (rhoN - rho0) / ((double) n + 1);
void TEST_JACOBI_ALGORITHM()
/*
Compares eigenvalues computed by the Jacobi_Algorithm method
with armadillo eig_gen method
 */
{ 
  int N = 5;
  arma::mat A = arma::zeros <arma::mat> (N, N);
  arma::mat E = arma::eye <arma::mat> (N, N);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);

  //Finding eigenvalues with armadillo
  arma::cx_vec eig_val;
  arma::cx_mat eig_vec;
  arma::eig_gen(eig_val, eig_vec, A);
  
  //Finding eigenvalues with Jacobi algorithm
  Jacobi_Algorithm(A, E, N, h);
  arma::vec calculated_eig_vals = A.diag();
  assert(arma::norm(arma::sort(A.diag()) - arma::sort(arma::real(eig_val))) <= eps);
}

void TEST_OFFMAX()
/*
Verifies that max_offdiag method chooses
the correct maximum value and checks indices
 */
{ 
  double t = -1e3; //Number with larger fabs than 0
  int N = 5;
  int k, l;
  arma::mat T = arma::zeros <arma::mat> (N, N); //Making matrix of zeros
  T(2, 1)     = t; // Setting one element to high value
  double max_val = max_offdiag(T, k, l, N);
  assert(k == 2 && l == 1 && max_val == fabs(t) );
}

void TEST_OFFDIAG_IS_ZERO()
/*
Verifies that the offdiagonal elements is zero after 
calling Jacobi_Algorithm method.
 */
{ 
  int N = 5;
  arma::mat A = arma::zeros <arma::mat> (N, N);
  arma::mat E = arma::eye <arma::mat> (N, N);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);
  Jacobi_Algorithm(A, E, N, h);
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      if (i != j){
        assert(fabs(A(i, j) * A(i, j)) <= eps);
      }
    }
  }
}

void TEST_ORTHOGONALITY()
/*
Verifies that the columns in A are orthogonal 
after calling the Jacobi_Algorithm method.
 */
{
  int N = 3;
  arma::mat A = arma::zeros <arma::mat> (N, N);
  arma::mat E = arma::eye <arma::mat> (N, N);
  A.diag(0).fill(2);
  A.diag(1).fill(-1);
  A.diag(-1).fill(-1);
  Jacobi_Algorithm(A, E, N, h);
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      if (i != j){
        assert(arma::dot(E.col(i), E.col(j)) <= eps);
      }
    }
  }
}