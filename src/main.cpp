#include <armadillo>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <iomanip>
#include <math.h> // tigonometric functions
#include <ctime>  // Run times  
#include <cmath>
#include <stdexcept>

#include "jacobi_algorithm.cpp"

#define PI 3.14159265

using namespace std;

inline double analytical_beam(int j, int N, double d , double a)
/* Compute analytical eigenvalues of buckling beam problem.

Parameters
------------
j: int
  j-th eigenvalue of buckling beam equation.
N: int
  Grid size.
d: double
  Diagonal elements.
a: double
  Upper and lower diagonal elements
*/
{
  return d + 2 * a * cos((double) j * PI / ((double) N + 1));
}

 
int main(int argc, char* argv[])
{ 
  int n;        // Grid size
  double eps;   // Tolerance
  double rhoN;  // Upper limit of rho

  if( argc != 4 ) {
    throw invalid_argument("Program takes three variables only through command line. Please insert dimension n, toleranse epsilon and rhoN.");
  }
  else {
    n = atoi(argv[1]);
    eps = atof(argv[2]);
    rhoN = atof(argv[3]);
  }
  double omega_r;   // Harmonic oscillator potential  
  double h;         // Step size
  double d;         // Diagonal elements
  double a;         // Upper and lower diagonal elements
  double rho0 = 0.0;  // Lower limit in rho
  int iterations;
  
  // Buckling beam problem:

  // Writing CPU time numerical and analytical eigenvalues to files
  ofstream outfile1;  
  ofstream outfile2;

  outfile1.open("bbeam_num.dat");
  outfile2.open("bbeam_analytical.dat");

  outfile1 << setw(20) << "N"
           << setw(20) << "lambda1" 
           << setw(20) << "lambda2"
           << setw(20) << "lambda 3" 
           << setw(20) << "CPUtime" 
           << setw(20) << "# Iterations" << endl;

  outfile2 << setw(20) << "N"
           << setw(20) << "lambda1" 
           << setw(20) << "lambda2"
           << setw(20) << "lambda3" << endl;

  for (int i = 50; i <= 400; i += 5 )
  {
    iterations = 0;
    h = (rhoN - rho0) / ((double) i + 1);
    d = 2. / (h * h);
    a = -1. / (h * h);
    arma::vec analy_eig = arma::zeros <arma::vec> (3);
    for (int j = 1; j <= 3; j++)
    {
      analy_eig(j-1) = analytical_beam(j, i, d, a);
    }

    arma::mat A = arma::zeros <arma::mat> (i, i);
    arma::mat E = arma::eye <arma::mat> (i, i);
    A.diag(0).fill(d);
    A.diag(1).fill(a);
    A.diag(-1).fill(a);

    clock_t t_start = clock();
    Jacobi_Algorithm(A, E, i, h, eps, iterations);
    clock_t t_end = clock();

    double CPU_time = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
    arma::vec num_eig = arma::sort(A.diag(0));

    outfile1 << setprecision(12) << setw(20) << i
             << setprecision(12) << setw(20) << num_eig(0)
             << setprecision(12) << setw(20) << num_eig(1)
             << setprecision(12) << setw(20) << num_eig(2)
             << setprecision(12) << setw(20) << CPU_time 
             << setprecision(12) << setw(20) << iterations << endl;

    outfile2 << setprecision(12) << setw(20) << i
             << setprecision(12) << setw(20) << analy_eig(0)
             << setprecision(12) << setw(20) << analy_eig(1)
             << setprecision(12) << setw(20) << analy_eig(2) << endl;

  }
  outfile1.close();
  outfile2.close();

  // One electron in a harmonic oscillator potential:

  h = (rhoN - rho0) / ((double) n + 1);
  d = 2. / (h * h);
  a = -1. / (h * h);

  arma::vec V = arma::zeros <arma::vec> (n);  
  Harmonic_Potential(V, rho0, n, h);
  arma::mat P = arma::zeros <arma::mat> (n, n);
  int length = 100;  
  arma:: vec rho_max = arma::linspace<arma::vec>(1, 10, length);
  double vars[5][length];
  double progress = 0;

  // Iterating over several rhoN for constant grid size n,
  // to find best rhoN:

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
    P.diag(0) += d + V; // Adding potential along diagonal
    P.diag(1).fill(a);
    P.diag(-1).fill(a);
    Jacobi_Algorithm(P, E, n, h, eps);

    arma:: vec diags = arma::sort(P.diag(0)); // Sorted eigenvalues

    vars[0][i] = rhoN;
    vars[1][i] = diags(0);
    vars[2][i] = diags(1);
    vars[3][i] = diags(2);
    vars[4][i] = diags(3);
    progress += 1;
    cout << "iterations left: "<< progress << "/" << length << endl;
  }
  // Writing first four eigenvalues to file
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
  
  // Coulomb interaction between two electrons
  // in a harmonic oscillator potential: 

  double omega_lin[6] = {0.01, 0.05, 0.25, 0.5, 1, 5};
  arma::mat P_col = arma::zeros <arma::mat> (n, n);
  arma::vec V_col = arma::zeros <arma::vec> (n);
  arma::mat E_col = arma::eye <arma::mat> (n, n);
  arma::vec energy_diag = arma::zeros <arma::vec> (n);
  arma::mat ground_states = arma::zeros <arma::mat> (n, 6);
  arma::mat lambdas = arma::zeros <arma::mat> (n, 6);
  arma::uword min_pos;

  // Calculating eigenvalue of ground state for different omega_r:

  for (int i = 0; i < 6; i++)
  {
    coulomb_potential(V_col, rho0, h, n, omega_lin[i]);
    P_col.diag(0) += d + V_col;   // Adding potential along diagonal
    P_col.diag(1).fill(a);
    P_col.diag(-1).fill(a);
    Jacobi_Algorithm(P_col, E_col, n, h, eps);
    energy_diag = arma::sort(P_col.diag(0));  // Sorted energy eigenvalues
    min_pos = P_col.diag(0).index_min();      // Extracting index of ground state
    ground_states.col(i) = E_col.col(min_pos);
    lambdas.col(i) = P_col.diag(0);
    cout << omega_lin[i] << " " << energy_diag(0) << endl;
    P_col.fill(0);
    E_col.fill(0);
    E_col.diag(0).fill(1);  // Resetting basis matrix
  }
  // Saving ground states and ground state energies in files
  ground_states.save("ground_states.txt", arma::arma_ascii);
  lambdas.save("lambdas.txt", arma::arma_ascii);
  // Write omega_r and rhoN to file
  ofstream outfile;
  outfile.open("omegas.txt");
  outfile << rhoN << endl;
  for (int i = 0; i < 6; i++)
  {
    outfile << omega_lin[i] << " ";
  }
  outfile << endl;
  outfile.close();
  
  return 0;
}
