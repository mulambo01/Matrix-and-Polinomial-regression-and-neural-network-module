#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../modules/matrix.h"

int g(mtx w, mtx input){
 int qtw=w.ncols;
 long double u;
 u=vectprod(input, w);
 if(u>=0.0)return 1;
 else return -1;
}
void main(){
 int qtw, qtspl;
 mtx w, input;
 qtw=5;
 qtspl=35;
 input=nullmatrix(1,5);
 w=mtxload("neuron.dat", 1, qtw);
 input.data[0][0]=-1.0;
 for(int i=1; i<qtw; i++){
  scanf("%Lf", &input.data[0][i]);
 }
 printf("\n%d\n", g(w, input));
}