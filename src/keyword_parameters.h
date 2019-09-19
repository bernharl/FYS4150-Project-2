int n = 400;
double rho0 = 0;
double rhoN = 5.4;
double h = (rhoN - rho0) / ((double) n + 1);
double eps = 1e-8;
double d = 2. / (h * h);
double a = -1. / (h * h);
int tot_iterations = 3 * n * n;

