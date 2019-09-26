#include <armadillo>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <stdexcept>

#include "jacobi_algorithm.h"

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
  double a_ll, a_kk, a_ik, a_il, a_kl, e_ik, e_il;  // Offdiagonal elements of A and elements of basis matrix E
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