double rhoN = 5.4;
double h = (rhoN - rho0) / ((double) n + 1);
double eps = 1e-8;
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
