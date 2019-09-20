#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <assert.h>


double max_offdiag(const arma::mat &A, int &k, int &l, int n);
void Jacobi_Algorithm(arma::mat &A, arma::mat &E, int n, double h, double eps);
void coulomb_potential(arma::vec &V, double rho0, double h, int n, double omega_r);
void Harmonic_Potential(arma::vec &V, double rho0, int n, double h);
 
#endif 