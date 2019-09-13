#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <assert.h>

int n = 5;
double h = 1. / (double) n;
double eps = 1e-8;
double d = 2. / (h * h)* 0 + 2.;
double a = -1. / (h * h) * 0 +1.;


double max_offdiag(const arma::mat &A, int &k, int &l, int n);
void Jacobi_Algorithm(arma::mat &A);


#endif