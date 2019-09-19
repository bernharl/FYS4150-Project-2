#include <armadillo>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <stdexcept>

#include "jacobi_algorithm.cpp"
using namespace std;


 
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
  //cout << A.diag(0) << endl;

 
  arma::vec V = arma::zeros <arma::vec> (n);
  Harmonic_Potential(V, rho0, n, h);
  arma::mat P = arma::zeros <arma::mat> (n, n);
  P.diag(0) += d + V;
  P.diag(1).fill(a);
  P.diag(-1).fill(a);

  
  Jacobi_Algorithm(P, E, n, h, eps);
  arma:: vec diags = arma::sort(P.diag(0));
  //cout << diags(0) << " " << diags(1) << " " 
  //     << diags(2) << " " << diags(3) << endl;



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
