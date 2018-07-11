#include <stdio.h>
#include <stdlib.h>

typedef struct
{
 int nrows;
 int ncols;
 long double **data;
}mtx;


mtx nullmatrix(int nrows, int ncols){
 mtx matrix;
 matrix.nrows=nrows;
 matrix.ncols=ncols;
 matrix.data=(long double **)malloc(nrows*sizeof(long double *));
 for(int i=0; i<nrows; i++){
  matrix.data[i]=(long double *)calloc(ncols, sizeof(long double));
 }
 return matrix;
}

mtx randmatrix(int nrows, int ncols){
 int j;
 mtx ans;
 ans.nrows=nrows;
 ans.ncols=ncols;

 long double **matrix=(long double **)malloc(nrows*sizeof(long double *));
 for(int i=0; i<nrows; i++){
  matrix[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0; j<ncols; j++){
   matrix[i][j]=(rand()%1000000)/1000000.0;
  }
 }

 ans.data=matrix;
 return ans;
}

//draw an nrows x ncols matrix on the console
void draw(mtx matrix){
 int j;
 for(int i=0; i<matrix.nrows ; i++){
  for(j=0; j<matrix.ncols ; j++){
   printf("%Lf\t", matrix.data[i][j]);
  }
  printf("\n");
 }
}
//load a matrix stored in a text file
//split columns with tab or space and rows with breaklines \n
mtx load(char *filename, int nrows, int ncols){
 mtx ans;
 ans.nrows=nrows;
 ans.ncols=ncols;
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
 fclose(file);
 ans.data=matrix;
 return ans;
}
void savematrix(char *filename, mtx matrix){
 int j;
 FILE *file;
 file=fopen(filename, "w");
 for(int i=0; i<matrix.nrows; i++){
  for(j=0; j<matrix.ncols; j++){
   fprintf(file, "%.20Le", matrix.data[i][j]);
   if(j!=matrix.ncols-1)fprintf(file, " ");
  }
  fprintf(file, "\n");
 }
 fclose(file);
}

//matrix sum, it will not verify the dimensions of matrix
mtx mtxsum(mtx matrix1, mtx matrix2){
 mtx sum;
 sum.nrows=matrix1.nrows;
 sum.ncols=matrix1.ncols;
 long double **answ=(long double **)malloc(matrix1.nrows*sizeof(long double *));
 int j;
 for(int i=0; i<matrix1.nrows; i++){
  answ[i]=(long double *)malloc(matrix1.ncols*sizeof(long double));
  for(j=0;j<matrix1.ncols;j++){
   answ[i][j]=matrix1.data[i][j]+matrix2.data[i][j];
  }
 }
 sum.data=answ;
 return sum;
}
//matrix multiplication by a real
mtx mtxmult(mtx matrix, long double q){
 mtx mult;
 mult.nrows=matrix.nrows;
 mult.ncols=matrix.ncols;
 long double **answ=(long double **)malloc(matrix.nrows*sizeof(long double *));
 int j;
 for(int i=0; i<matrix.nrows; i++){
  answ[i]=(long double *)malloc(matrix.ncols*sizeof(long double));
  for(j=0;j<matrix.ncols;j++){
   answ[i][j]=q*matrix.data[i][j];
  }
 }
 mult.data=answ;
 return mult;
}
//matrix multiplication, it will not verify the dimension of matrix
mtx mtxprod(mtx matrix1, mtx matrix2){
 mtx prod;
 prod.nrows=matrix1.nrows;
 prod.ncols=matrix2.ncols;
 long double **answ=(long double **)malloc(prod.nrows*sizeof(long double *));
 int j,k;
 for(int i=0; i<prod.nrows; i++){
  answ[i]=(long double *)malloc(prod.ncols*sizeof(long double));
  for(j=0;j<prod.ncols;j++){
   answ[i][j]=0.0;
   for(k=0; k<matrix1.ncols; k++){
    answ[i][j]=answ[i][j]+matrix1.data[i][k]*matrix2.data[k][j];
   }
  }
 }
 prod.data=answ;
 return prod;
}
//return the transpose of a matrix
mtx transpose(mtx matrix){
 mtx answ;
 answ.ncols=matrix.nrows;
 answ.nrows=matrix.ncols;
 long double **trans=(long double **)malloc(answ.nrows*sizeof(long double *));
 int i,j;
 for(i=0; i<matrix.ncols; i++){
  trans[i]=(long double *)malloc(answ.ncols*sizeof(long double));
  for(j=0; j<matrix.nrows; j++){
   trans[i][j]=matrix.data[j][i];
  }
 }
 answ.data=trans;
 return answ;
}
//copy a matrix to another address
mtx matrixcopy(mtx A){
 mtx ans;
 int j,nrows, ncols;
 ans.nrows=A.nrows;
 ans.ncols=A.ncols;
 nrows=ans.nrows;
 ncols=ans.ncols;
 long double **A2=(long double **)malloc(nrows*sizeof(long double *));
 for(int i=0; i<nrows; i++){
  A2[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0; j<ncols; j++){
   A2[i][j]=A.data[i][j];
  }
 }
 ans.data=A2;
 return ans;
}

mtx cutmatrix(mtx matrix, int srowcut, int nrows, int scolumncut, int ncols){
 mtx result=nullmatrix(nrows, ncols);
 int j, line, column;
 for(int i=0; i<nrows; i++){
  for(j=0; j<ncols; j++){
   line=srowcut+i;
   column=scolumncut+j;
   result.data[i][j]=matrix.data[line][column];
  }
 }
 return result;
}

long double vectprod(mtx vec1, mtx vec2){
 long double result=0.0;
 int j;
 for(int i=0; i<vec1.nrows; i++){
  for(j=0; j<vec1.ncols; j++){
   result=result+vec1.data[i][j]*vec2.data[i][j];
  }
 }
 return result;
}
