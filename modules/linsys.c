#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

//solves an upper triangular system returning a column matrix
mtx uptsyst(mtx matrix, mtx y){
 int ncols=matrix.nrows;
 int nrows=ncols;
 mtx ans;
 ans.nrows=nrows;
 ans.ncols=1;
 long double **solut=(long double **)malloc(ncols*sizeof(long double *));
 int i,j;
 for(i=nrows-1; i>=0; i--){
  solut[i]=(long double *)malloc(sizeof(long double));
  solut[i][0]=y.data[i][0];
  for(j=i+1; j<ncols; j++){
   solut[i][0]=solut[i][0]-matrix.data[i][j]*solut[j][0];
  }
  solut[i][0]=solut[i][0]/matrix.data[i][i];
 }
 ans.data=solut;
 return ans;
 }
//solves a bottom triangular system returning a column matrix
mtx bottsyst(mtx matrix, mtx y){
 int ncols=matrix.nrows;
 int nrows=ncols;
 mtx ans;
 ans.nrows=nrows;
 ans.ncols=1;
 long double **solut=(long double **)malloc(ncols*sizeof(long double *));
 int i,j;
 for(i=0; i<nrows; i++){
  solut[i]=(long double *)malloc(sizeof(long double));
  solut[i][0]=y.data[i][0];
  for(j=0; j<i; j++){
   solut[i][0]=solut[i][0]-matrix.data[i][j]*solut[j][0];
  }
  solut[i][0]=solut[i][0]/matrix.data[i][i];
 }
 ans.data=solut;
 return ans;
 }

//solves a linear system by the LU decomposition
mtx LU(mtx A, mtx B){
 mtx mL, mU,X,Y;
 int n=A.nrows;
 long double **L=(long double **)malloc(n*sizeof(long double *));
 long double **U=(long double **)malloc(n*sizeof(long double *));
 long double max, change, r, m;
 int i, j, k, l, pivo, chpivo;
 mL.nrows=n;
 mU.nrows=n;
 mL.ncols=n;
 mU.ncols=n;
 for(i=0; i<n; i++){
  L[i]=(long double *)malloc(n*sizeof(long double *));
  U[i]=(long double *)malloc(n*sizeof(long double *));
 }
 U=matrixcopy(A).data;
 for(i=0; i<n; i++){
  pivo=i;
  max=fabsl(U[i][i]);
  for(j=i+1; j<n; j++){
   if(fabsl(U[j][i])>max){
    max=fabsl(U[j][i]);
    pivo=j;
   }
  }
  if(pivo!=i){
   for(k=0;k<n;k++){
    change=U[i][k];
    U[i][k]=U[pivo][k];
    U[pivo][k]=change;

    change=L[i][k];
    L[i][k]=L[pivo][k];
    L[pivo][k]=change;
   }
   change=B.data[i][0];
   B.data[i][0]=B.data[pivo][0];
   B.data[pivo][0]=change;
  }
  if(U[i][i]!=0.0){
   r=1.0/(U[i][i]);
   L[i][i]=1.0;
   for(k=i+1; k<n; k++){
    m=U[k][i]*r;
    L[k][i]=m;
    U[k][i]=0.0;
    for(l=i+1; l<n; l++){
     U[k][l]=U[k][l]-m*U[i][l];
    }
   }
  }
 }
 mL.data=L;
 mU.data=U;
 Y=bottsyst(mL,B);
 X=uptsyst(mU,Y);
 return X;
}
//solves a simetric linear system using cholesky decomposition
//it will not verify if the matrix A is simetric
mtx cholesky(mtx A, mtx B){
 int i,j,k,n;
 n=A.nrows;
 long double **L=(long double **)malloc(n*sizeof(long double *));
 long double **Lt=(long double **)malloc(n*sizeof(long double *));
 mtx X, Y, mL, mLt;
 mL.nrows=n;
 mL.ncols=n;
 for(i=0; i<n; i++){
  L[i]=(long double *)malloc(n*sizeof(long double));
  Lt[i]=(long double *)malloc(n*sizeof(long double));
  for(j=0; j<=i; j++){
   L[i][j]=A.data[i][j];
   for(k=0; k<j; k++){
    L[i][j]=L[i][j]-L[j][k]*L[i][k];
   }
   if(i==j){
    L[i][j]=sqrtl(L[i][j]);
   }
   else{
    L[i][j]=L[i][j]/L[j][j];
   }
  }
 }
 mL.data=L;
 mLt=transpose(mL);
 Y=bottsyst(mL,B);
 X=uptsyst(mLt,Y);
 return X;
}


//find a polynomial regression equation, returning the coefficients
//in a column matrix
mtx regression(mtx points, int n){
 mtx matrix, Y, B;
 long double sumX, sumY, xn;
 int i,j,k,npoints;
 npoints=points.nrows;
 matrix=nullmatrix(n+1,n+1);
 Y=nullmatrix(n+1,1);
 B=nullmatrix(n+1,1);
 for(i=0; i<=n; i++){
  sumY=0.0;
  for(j=0; j<=i; j++){
   sumX=0.0;
   for(k=0; k<npoints; k++){
    xn=powl(points.data[k][0], i+j);
    if(j==0){
     sumY=sumY+xn*points.data[k][1];
    }
    sumX=sumX+xn;
   }
   Y.data[i][0]=sumY;
   matrix.data[i][j]=sumX;
  }
 }
 B=cholesky(matrix, Y);
 return B;
}
//draw the equation on the output, receiving the coefficients in a
//column matrix
void draweq(mtx coeffs){
 int n=coeffs.nrows-1;
 for(int j=0; j<=n ; j++){
  if(j>0 && coeffs.data[j][0]>0){
   printf("+");
  }
  printf("%.17Lf*x**%d", coeffs.data[j][0],j);
 }
 printf("\n");
}
//return a value of a polynomial function, receiving the coefficients
//in a column matrix
long double func(mtx coeffs, long double x, int degree){
 long double y=coeffs.data[0][0];
 for(int i=1; i<=degree; i++){
  y=y+coeffs.data[i][0]*powl(x,i);
 }
 return y;
}
//return the value of r squared
long double r2(mtx coeffs, mtx samples, int degree){
 long double r2, d=0.0, y1=0.0, y2=0.0, y, x;
 int n=samples.nrows;
 for(int i=0; i<n; i++){
  x=samples.data[i][0];
  y=samples.data[i][1];
  y1=y1+y;
  y2=y2+y*y;
  d=d+powl((y-func(coeffs, x, degree)),2);
 }
 r2=1.0-d/(y2-y1*y1/n);
 return r2;
}
//return the value of the residual variance
long double residvariance(mtx coeffs, mtx samples, int degree){
 long double residvar, d=0.0, y, x;
 int n=samples.nrows;
 for(int i=0; i<n; i++){
  x=samples.data[i][0];
  y=samples.data[i][1];
  d=d+powl((y-func(coeffs, x, degree)),2);
 }
 residvar=d/(n-degree-1);
 return residvar;
}
