#include <stdio.h>
#include <stdlib.h>

//struct of an object called matrix
typedef struct
{
 int nrows;
 int ncols;
 long double **data;
}mtx;

//return a matrix with equals elements of value
mtx crystalmatrix(int nrows, int ncols, long double value){
 mtx matrix;
 int j;
 matrix.nrows=nrows;
 matrix.ncols=ncols;
 matrix.data=(long double **)malloc(nrows*sizeof(long double *));
 for(int i=0; i<nrows; i++){
  matrix.data[i]=(long double *)malloc(ncols*sizeof(long double));
  for(j=0; j<ncols; j++){
   matrix.data[i][j]=value;
  }
 }
 return matrix; 
}
//return a null matrix
mtx nullmatrix(int nrows, int ncols){
 return crystalmatrix(nrows, ncols, 0.0);
}
//return a matrix of elements with random values
mtx randmatrix(int nrows, int ncols, unsigned int seed){
 int j;
 mtx ans;
 ans.nrows=nrows;
 ans.ncols=ncols;
 srand(seed);

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

//draw a matrix on console
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
mtx mtxload(char *filename, int nrows, int ncols){
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
//save a matrix in a text file
void mtxsave(char *filename, mtx matrix){
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
//matrix subtraction
mtx mtxsub(mtx matrix1, mtx matrix2){
 mtx sum;
 sum.nrows=matrix1.nrows;
 sum.ncols=matrix1.ncols;
 long double **answ=(long double **)malloc(matrix1.nrows*sizeof(long double *));
 int j;
 for(int i=0; i<matrix1.nrows; i++){
  answ[i]=(long double *)malloc(matrix1.ncols*sizeof(long double));
  for(j=0;j<matrix1.ncols;j++){
   answ[i][j]=matrix1.data[i][j]-matrix2.data[i][j];
  }
 }
 sum.data=answ;
 return sum;
}
//mutiply each term of matrix1 by the corresponding
//will not verify the dimensions
mtx mtxtermsmult(mtx matrix1, mtx matrix2){
 int j;
 mtx result=nullmatrix(matrix1.nrows, matrix1.ncols);
 for(int i=0; i<matrix1.nrows; i++){
  for(j=0; j<matrix1.ncols; j++){
   result.data[i][j]=matrix1.data[i][j]*matrix2.data[i][j];
  }
 }
 return result;
}
//multiply each term of matrix by the value of q
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
//matrix multiplication
//will not verify the dimension of matrix
//its not the most efficient method
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
//return a copy of matrix A
mtx mtxclone(mtx A){
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

//it will cut the matrix starting in position srowcut:scolumncut with nrows lines and ncols columns ahead
mtx mtxcut(mtx matrix, int srowcut, int nrows, int scolumncut, int ncols){
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

//it will multiply each term by the corresponding and sum
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

//it will copy matrix2 to the address of matrix1
//it will realloc matrix1
void mtxcopy(mtx *matrix1, mtx matrix2){
 int j, difsize=0;
 if(matrix1->nrows != matrix2.nrows){
  matrix1->data=(long double **)realloc(matrix1->data, matrix2.nrows*sizeof(long double *));
  matrix1->nrows=matrix2.nrows;
 }
 if(matrix1->ncols != matrix2.ncols){
  difsize=1;
  matrix1->ncols=matrix2.ncols;
 }
 for(int i=0; i<matrix2.nrows; i++){
  if(difsize){
   matrix1->data[i]=(long double *)realloc(matrix1->data[i], matrix2.ncols*sizeof(long double));
  }
  for(j=0; j<matrix2.ncols; j++){
   matrix1->data[i][j]=matrix2.data[i][j];
  }
 }
}
//it will free the buffer of matrix
void mtxfree(mtx *matrix){
 for(int i=0; i<matrix->nrows; i++){
  free(matrix->data[i]);
 }
 free(matrix->data);
}
//it will replace the line "line" of matrix by the line matrix "line"
//will not verify the dimensions
void putline(mtx *matrix, mtx line, int posi){
 int ncols=matrix->ncols;
 for(int i=0; i<ncols; i++){
  matrix->data[posi][i]=line.data[0][i];
 }
}
//it will fit the line matrix "line" in the position posi of matrix
void addline(mtx *matrix, mtx line, int posi){
 mtx bak;
 bak=mtxclone(*matrix);
 int j=0;
 mtxfree(matrix);
 *matrix=nullmatrix(bak.nrows+1, bak.ncols);
 for(int i=0; i<bak.nrows+1; i++){
  if(i==posi){
   putline(matrix, line, i);
  }
  else{
   putline(matrix, mtxcut(bak, j, 1, 0, bak.ncols), i);
   j++;
  }
 }
}
//do the same with columns
void putcol(mtx *matrix, mtx col, int posi){
 int nrows=matrix->nrows;
 for(int i=0; i<nrows; i++){
  matrix->data[i][posi]=col.data[i][0];
 }
}
void addcol(mtx *matrix, mtx col, int posi){
 mtx bak;
 bak=mtxclone(*matrix);
 int j=0;
 mtxfree(matrix);
 *matrix=nullmatrix(bak.nrows, bak.ncols+1);
 for(int i=0; i<bak.ncols+1; i++){
  if(i==posi){
   putcol(matrix, col, i);
  }
  else{
   putcol(matrix, mtxcut(bak, 0, bak.nrows, j, 1), i);
   j++;
  }
 }
}