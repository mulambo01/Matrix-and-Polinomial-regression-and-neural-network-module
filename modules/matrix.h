typedef struct
{
 int nrows;
 int ncols;
 long double **data;
}mtx;

mtx nullmatrix(int nrows, int ncols);
mtx randmatrix(int nrows, int ncols);
void draw(mtx matrix);
mtx load(char *filename, int nrows, int ncols);
void savematrix(char *filename, mtx matrix);
mtx mtxsum(mtx matrix1, mtx matrix2);
mtx mtxmult(mtx matrix, long double q);
mtx mtxprod(mtx matrix1, mtx matrix2);
mtx transpose(mtx matrix);
mtx matrixcopy(mtx A);
mtx cutmatrix(mtx matrix, int srowcut, int nrows, int scolumncut, int ncols);
