#include "matrix.h"

mtx uptsyst(mtx matrix, mtx y);
mtx bottsyst(mtx matrix, mtx y);
mtx LU(mtx A, mtx B);
mtx cholesky(mtx A, mtx B);

mtx regression(mtx points, int n);
void draweq(mtx coefs);
long double func(mtx coefs, long double x, int degree);
long double r2(mtx coefs, mtx samples, int degree);
long double residvariance(mtx coefs, mtx samples, int degree);
