typedef struct
{
 int nrows;
 int ncols;
 long double **data;
}mtx;

mtx crystalmatrix(int nrows, int ncols, long double value);
mtx nullmatrix(int nrows, int ncols);
mtx randmatrix(int nrows, int ncols);
void draw(mtx matrix);
mtx mtxload(char *filename, int nrows, int ncols);
void mtxsave(char *filename, mtx matrix);
mtx mtxsum(mtx matrix1, mtx matrix2);
mtx mtxsub(mtx matrix1, mtx matrix2);
mtx mtxtermsmult(mtx matrix1, mtx matrix2);
mtx mtxmult(mtx matrix, long double q);
mtx mtxprod(mtx matrix1, mtx matrix2);
mtx transpose(mtx matrix);
mtx mtxclone(mtx A);
mtx mtxcut(mtx matrix, int srowcut, int nrows, int scolumncut, int ncols);
long double vectprod(mtx vec1, mtx vec2);
void mtxcopy(mtx *matrix1, mtx matrix2);
void mtxfree(mtx *matrix);
void putline(mtx *matrix, mtx line, int posi);
void addline(mtx *matrix, mtx line, int posi);
void putcol(mtx *matrix, mtx col, int posi);
void addcol(mtx *matrix, mtx col, int posi);
mtx *mtxsplitlines(mtx matrix);
