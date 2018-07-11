#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void draw(long double **matrix, int nrows, int ncols);
long double **load(char *filename, int nrows, int ncols);
void savematrix(char *filename, long double **matrix ,int nrows, int ncols);
long double **mtxsum(long double **matrix1, long double **matrix2, int nrows, int ncols);
long double **mtxmult(long double **matrix, long double q, int nrows, int ncols);
long double **mtxprod(long double **matrix1, long double **matrix2, int n, int nrows, int ncols);
long double** uptsyst(long double **matrix, long double **y, int nrows);
long double** bottsyst(long double **matrix, long double **y, int nrows);
long double** transpose(long double **matrix, int nrows, int ncols);
long double** matrixcopy(long double **A, int nrows, int ncols);
long double** LU(long double **A, long double **B, int n);
long double **cholesky(long double **A, long double **B, int n);
long double **regression(long double **points, int npoints, int n);
void draweq(long double **coefs, int n);
long double func(long double **coefs, long double x, int degree);
long double r2(long double **coefs, long double **samples, int degree, int n);
long double residvariance(long double **coefs, long double **samples, int degree, int n);
long double **cutmatrix(long double **matrix, int srowcut, int nrows, int scolumncut, int ncols);
long double **randmatrix(int nrows, int ncols);