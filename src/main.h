#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <assert.h>

int n;
double h;
double d;
double a;
double eps;
double rhoN;
double rho0 = 0.0;

double max_offdiag(const arma::mat &A, int &k, int &l, int n);
void Jacobi_Algorithm(arma::mat &A, arma::mat &E, int n, double h, double eps);

 
#endif 