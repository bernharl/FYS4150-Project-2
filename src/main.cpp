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
    throw invalid_argument("Program takes three variables only through command line. Please insert dimension n, toleranse epsilon and rhoN.");
  }
  else{
    n = atoi(argv[1]);
    eps = atof(argv[2]);
    rhoN = atof(argv[3]);
  }
  double omega_r;
  double h;
  double d;
  double a;
  double rho0 = 0.0;


  h = (rhoN - rho0) / ((double) n + 1);
  d = 2. / (h * h);
  a = -1. / (h * h);
  /*
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
  */


  int length = 100;
  arma:: vec rho_max = arma::linspace<arma::vec>(1 ,10, length);

  double vars[5][length];

  clock_t t_start = clock(); // Initializing timer
  /*
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
  */
  length = 10;
  omega_r = 1./4;
  //arma:: vec omega_lin = arma::linspace<arma::vec>(0.01 , 5, );
  double omega_lin[6] = {0.01, 0.05, 0.25, 0.5, 1, 5};
  arma::mat P_col = arma::zeros <arma::mat> (n, n);
  arma::vec V_col = arma::zeros <arma::vec> (n);
  arma::mat E_col = arma::eye <arma::mat> (n, n);
  arma::vec energy_diag = arma::zeros <arma::vec> (n);
  arma::mat ground_states = arma::zeros <arma::mat> (n, 6);

  for (int i = 0; i < 6; i++)
  {
    coulomb_potential(V_col, rho0, h, n, omega_lin[i]);
    P_col.diag(0) += d + V_col;
    P_col.diag(1).fill(a);
    P_col.diag(-1).fill(a);
    Jacobi_Algorithm(P_col, E_col, n, h, eps);
    energy_diag = arma::sort(P_col.diag(0));
    ground_states.col(i) = E_col.col(0);
    cout << omega_lin[i] << " " << energy_diag(0) << endl;
    P_col.fill(0);
    E_col.fill(0);
    E_col.diag(0).fill(1);
  }
  ground_states.save("ground_states.txt", arma::arma_ascii);
  return 0;
}
