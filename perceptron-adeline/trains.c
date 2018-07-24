#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../modules/matrix.h"

long double meansqrerr(mtx w, mtx d, mtx samples){
 long double result=0.0;
 int qtspl=samples.nrows;
 int qtw=w.ncols;
 mtx u;
 u=mtxprod(samples, transpose(w));
 for(int i=0; i<qtspl; i++){
  result=result+powl((d.data[i][0]-u.data[i][0]),2);
 }
 result=result/qtspl;
 return result;
}

void main(){
 int qtw, qtspl, steps,i;
 long double N, precis, E, lastE, u;
 mtx samples, d, w, x;
 qtw=5;
 qtspl=35;
 N=0.0025;
 precis=powl(2.0,-65);
 samples=mtxload("Table.dat", qtspl, qtw+1);
 d=mtxcut(samples,0,qtspl,qtw,1);
 samples=mtxcut(samples,0,qtspl,0,qtw);
 w=randmatrix(1,qtw, 1);
// w=mtxload("w2.dat",1,qtw);
 steps=0;
 lastE=0.0;
 E=meansqrerr(w,d,samples);
 while(fabsl(E-lastE)>precis){
  lastE=E;
  for(i=0; i<qtspl; i++){
   x=mtxcut(samples, i, 1, 0, qtw);
   u=vectprod(x, w);
   x=mtxmult(x, N*(d.data[i][0]-u));
   w=mtxsum(w,x);
  }
  E=meansqrerr(w,d,samples);
  steps=steps+1;
 }
 mtxsave("neuron.dat", w);
 printf("number of steps: %d\n", steps);
}