#include <armadillo>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <stdexcept>

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

void Jacobi_Algorithm(arma::mat &A, arma::mat &E, int n, double h, double eps)
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
eps: double
  Tolerance for treating off-diagonal elements in A as zero
 */
{

  int k, l;
  double max_val = max_offdiag(A, k, l, n);
  double a_ll, a_kk, a_ik, a_il, a_kl, e_ik, e_il;
  double t_val, tau, c, s; 
  double iterator = 0;
  int tot_iterations = 3 * n * n;

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

void Harmonic_Potential(arma::vec &V, double rho0, int n, double h)
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


//hallabrur

int main(int argc, char* argv[])
{ 
  int n;
  double eps;
  double rhoN;

  if( argc != 4 ){
    throw invalid_argument("Program takes two variables only through command line. Please insert dimension n and toleranse epsilon and rhoN.");
  }
  else{
    n = atoi(argv[1]);
    eps = atof(argv[2]);
    rhoN = atof(argv[3]);
  }

  double h;
  double d;
  double a;
  double rho0 = 0.0;


  h = (rhoN - rho0) / ((double) n + 1);
  d = 2. / (h * h);
  a = -1. / (h * h);

  arma::mat A = arma::zeros <arma::mat> (n, n);
  arma::mat E = arma::eye <arma::mat> (n, n);
  A.diag(0).fill(d);
  A.diag(1).fill(a);
  A.diag(-1).fill(a);

  Jacobi_Algorithm(A, E, n, h, eps);
  cout << A.diag(0) << endl;

 
  arma::vec V = arma::zeros <arma::vec> (n);
  Harmonic_Potential(V, rho0, n, h);
  arma::mat P = arma::zeros <arma::mat> (n, n);
  P.diag(0) += d + V;
  P.diag(1).fill(a);
  P.diag(-1).fill(a);

  
  Jacobi_Algorithm(P, E, n, h, eps);
  arma:: vec diags = arma::sort(P.diag(0));
  cout << diags(0) << " " << diags(1) << " " 
       << diags(2) << " " << diags(3) << endl;



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
    Harmonic_Potential(V, rho0, n, h);
    arma::mat E = arma::eye <arma::mat> (n, n);
    P.diag(0) += d + V;
    P.diag(1).fill(a);
    P.diag(-1).fill(a);
    Jacobi_Algorithm(P, E, n, h, eps);

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
