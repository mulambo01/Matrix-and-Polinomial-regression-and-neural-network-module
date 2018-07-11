#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//draw an nrows x ncols matrix on the console
void draw(long double **matrix, int nrows, int ncols){
 int j;
 for(int i=0; i<nrows ; i++){
  for(j=0; j<ncols ; j++){
   printf("%Lf\t", matrix[i][j]);
  }
  printf("\n");
 }
}
//load a matrix stored in a text file
//split columns with tab or space and rows with breaklines \n
long double **load(char *filename, int nrows, int ncols){
 long double **matrix=(long double **)malloc(nrows * sizeof(long double *));
 int j;
 FILE *file;
 file=fopen(filename, "r");
 for(int i=0; i<nrows; i++){
  matrix[i]=(long double *)malloc(ncols * sizeof(long double));
  for(j=0; j<ncols; j++){
   fscanf(file, "%Lf", &matrix[i][j]);
  }
 }
 return matrix;
}
void savematrix(char *filename, long double **matrix ,int nrows, int ncols){
 int j;
 FILE *file;
 file=fopen(filename, "w");
 for(int i=0; i<nrows; i++){
  for(j=0; j<ncols; j++){
   fprintf(file, "%.20Le", matrix[i][j]);
   if(j!=ncols-1)fprintf(file, " ");
  }
  fprintf(file, "\n");
 }
}

//matrix sum, it will not verify the dimensions of matrix
long double **mtxsum(long double **matrix1, long double **matrix2, int nrows, int ncols){
 long double **answ=(long double **)malloc(nrows*sizeof(long double *));
 int j;
 for(int i=0; i<nrows; i++){
  answ[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0;j<ncols;j++){
   answ[i][j]=matrix1[i][j]+matrix2[i][j];
  }
 }
 return answ;
}
//matrix multiplication by a real
long double **mtxmult(long double **matrix, long double q, int nrows, int ncols){
 long double **answ=(long double **)malloc(nrows*sizeof(long double *));
 int j;
 for(int i=0; i<nrows; i++){
  answ[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0;j<ncols;j++){
   answ[i][j]=q*matrix[i][j];
  }
 }
 return answ;
}
//matrix multiplication, it will not verify the dimension of matrix
long double **mtxprod(long double **matrix1, long double **matrix2, int n, int nrows, int ncols){
 long double **answ=(long double **)malloc(nrows*sizeof(long double *));
 int j,k;
 for(int i=0; i<nrows; i++){
  answ[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0;j<ncols;j++){
   answ[i][j]=0.0;
   for(k=0; k<n; k++){
    answ[i][j]=answ[i][j]+matrix1[i][k]*matrix2[k][j];
   }
  }
 }
 return answ;
}
//solves an upper triangular system returning a column matrix
long double** uptsyst(long double **matrix, long double **y, int nrows){
 int ncols=nrows;
 long double **solut=(long double **)malloc(ncols*sizeof(long double *));
 int i,j;
 for(i=nrows-1; i>=0; i--){
  solut[i]=(long double *)malloc(sizeof(long double));
  solut[i][0]=y[i][0];
  for(j=i+1; j<ncols; j++){
   solut[i][0]=solut[i][0]-matrix[i][j]*solut[j][0];
  }
  solut[i][0]=solut[i][0]/matrix[i][i];
 }
 return solut;
 }
//solves a bottom triangular system returning a column matrix
long double** bottsyst(long double **matrix, long double **y, int nrows){
 int ncols=nrows;
 long double **solut=(long double **)malloc(ncols*sizeof(long double *));
 int i,j;
 for(i=0; i<nrows; i++){
  solut[i]=(long double *)malloc(sizeof(long double));
  solut[i][0]=y[i][0];
  for(j=0; j<i; j++){
   solut[i][0]=solut[i][0]-matrix[i][j]*solut[j][0];
  }
  solut[i][0]=solut[i][0]/matrix[i][i];
 }
 return solut;
 }
//return the transpose of a matrix
long double** transpose(long double **matrix, int nrows, int ncols){
 long double **trans=(long double **)malloc(ncols*sizeof(long double *));
 int i,j;
 for(i=0; i<ncols; i++){
  trans[i]=(long double *)malloc(nrows*sizeof(long double));
  for(j=0; j<nrows; j++){
   trans[i][j]=matrix[j][i];
  }
 }
 return trans;
}
//copy a matrix to another address
long double** matrixcopy(long double **A, int nrows, int ncols){
 long double **A2=(long double **)malloc(nrows*sizeof(long double *));
 int j;
 for(int i=0; i<nrows; i++){
  A2[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0; j<ncols; j++){
   A2[i][j]=A[i][j];
  }
 }
 return A2;
}
//solves a linear system by the LU decomposition
long double** LU(long double **A, long double **B, int n){
 long double **X=(long double **)malloc(n*sizeof(long double *));
 long double **Y=(long double **)malloc(n*sizeof(long double *));
 long double **L=(long double **)malloc(n*sizeof(long double *));
 long double **U=(long double **)malloc(n*sizeof(long double *));
 long double max, change, r, m;
 int i, j, k, l, pivo, chpivo;
 for(i=0; i<n; i++){
  X[i]=(long double *)malloc(sizeof(long double));
  Y[i]=(long double *)malloc(sizeof(long double));
  L[i]=(long double *)malloc(n*sizeof(long double *));
  U[i]=(long double *)malloc(n*sizeof(long double *));
 }
 U=matrixcopy(A,n,n);
 for(i=0; i<n; i++){
  pivo=i;
  max=abs(U[i][i]);
  for(j=i+1; j<n; j++){
   if(abs(U[j][i])>max){
    max=abs(U[j][i]);
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
   change=B[i][0];
   B[i][0]=B[pivo][0];
   B[pivo][0]=change;
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
 Y=bottsyst(L,B,n);
 X=uptsyst(U,Y,n);
 return X;
}
//solves a simetric linear system using cholesky decomposition
//it will not verify if the matrix A is simetric
long double **cholesky(long double **A, long double **B, int n){
 int i,j,k;
 long double **L=(long double **)malloc(n*sizeof(long double *));
 long double **Lt=(long double **)malloc(n*sizeof(long double *));
 long double **X=(long double **)malloc(n*sizeof(long double *));
 long double **Y=(long double **)malloc(n*sizeof(long double *));
 for(i=0; i<n; i++){
  X[i]=(long double *)malloc(sizeof(long double));
  Y[i]=(long double *)malloc(sizeof(long double));
  L[i]=(long double *)malloc(n*sizeof(long double));
  Lt[i]=(long double *)malloc(n*sizeof(long double));
  for(j=0; j<=i; j++){
   L[i][j]=A[i][j];
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
 Lt=transpose(L,n,n);
 Y=bottsyst(L,B,n);
 X=uptsyst(Lt,Y,n);
 return X;
}
//find a polynomial regression equation, returning the coefficients
//in a column matrix
long double **regression(long double **points, int npoints, int n){
 long double **matrix=(long double **)malloc((n+1)*sizeof(long double *));
 long double **Y=(long double **)malloc((n+1)*sizeof(long double *));
 long double **B=(long double **)malloc((n+1)*sizeof(long double *));
 long double sumX, sumY, xn;
 int i,j,k;
 for(i=0; i<=n; i++){
  matrix[i]=(long double *)malloc((n+1)*sizeof(long double));
  Y[i]=(long double *)malloc(sizeof(long double));
  B[i]=(long double *)malloc(sizeof(long double));
  sumY=0.0;
  for(j=0; j<=i; j++){
   sumX=0.0;
   for(k=0; k<npoints; k++){
    xn=powl(points[k][0], i+j);
    if(j==0){
     sumY=sumY+xn*points[k][1];
    }
    sumX=sumX+xn;
   }
   Y[i][0]=sumY;
   matrix[i][j]=sumX;
  }
 }
 B=cholesky(matrix, Y, n+1);
 return B;
}
//draw the equation on the output, receiving the coefficients in a
//column matrix
void draweq(long double **coeffs, int n){
 for(int j=0; j<=n ; j++){
  if(j>0 && coeffs[j][0]>0){
   printf("+");
  }
  printf("%.17Lf*x**%d", coeffs[j][0],j);
 }
 printf("\n");
}
//return a value of a polynomial function, receiving the coefficients
//in a column matrix
long double func(long double **coeffs, long double x, int degree){
 long double y=coeffs[0][0];
 for(int i=1; i<=degree; i++){
  y=y+coeffs[i][0]*powl(x,i);
 } 
 return y;
}
//return the value of r squared
long double r2(long double **coeffs, long double **samples, int degree, int n){
 long double r2, d=0.0, y1=0.0, y2=0.0, y, x;
 for(int i=0; i<n; i++){
  x=samples[i][0];
  y=samples[i][1];
  y1=y1+y;
  y2=y2+y*y;
  d=d+powl((y-func(coeffs, x, degree)),2);
 }
 r2=1.0-d/(y2-y1*y1/n);
 return r2;
}
//return the value of the residual variance
long double residvariance(long double **coeffs, long double **samples, int degree, int n){
 long double stddev, d=0.0, y, x;
 for(int i=0; i<n; i++){
  x=samples[i][0];
  y=samples[i][1];
  d=d+powl((y-func(coeffs, x, degree)),2);
 }
 stddev=d/(n-degree-1);
 return stddev;
}

long double **cutmatrix(long double **matrix, int srowcut, int nrows, int scolumncut, int ncols){
 long double **result=(long double **)malloc(nrows*sizeof(long double *));
 int j, line, column;
 for(int i=0; i<nrows; i++){
  result[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0; j<ncols; j++){
   line=srowcut+i;
   column=scolumncut+j;
   result[i][j]=matrix[line][column];
  }
 }
 return result;
}

long double **randmatrix(int nrows, int ncols){
 int j;
 long double **matrix=(long double **)malloc(nrows*sizeof(long double *));
 for(int i=0; i<nrows; i++){
  matrix[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0; j<ncols; j++){
   matrix[i][j]=(rand()%1000000)/1000000.0;
  }
 }
 return matrix;
}