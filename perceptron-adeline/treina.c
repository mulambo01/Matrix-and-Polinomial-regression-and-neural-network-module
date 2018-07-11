#include "../modules/linsys.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

long double meansqrerr(mtx w, mtx d, mtx samples){
 long double result=0.0;
 int qtspl=samples.nrows;
 int qtw=w.ncols;
 mtx u;
 u=mtxprod(samples, transpose(w));
 printf("\n");
 draw(u);
 printf("\n");
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
 precis=10.e-50;
 samples=load("apendice2.txt", qtspl, qtw+1);
 d=cutmatrix(samples,0,qtspl,qtw,1);
 samples=cutmatrix(samples,0,qtspl,0,qtw);
 w=randmatrix(1,qtw);
// w=load("w2.dat",1,qtw);
 steps=0;
 lastE=0.0;
 E=meansqrerr(w,d,samples);
 while(fabsl(E-lastE)>precis){
  lastE=E;
  for(i=0; i<qtspl; i++){
   x=cutmatrix(samples, i, 1, 0, qtw);
   u=vectprod(x, w);
   x=mtxmult(x, N*(d.data[i][0]-u));
   w=mtxsum(w,x);
  }
  E=meansqrerr(w,d,samples);
  steps=steps+1;
 }
 savematrix("pesos.txt", w);
 printf("%d\n", steps);
}