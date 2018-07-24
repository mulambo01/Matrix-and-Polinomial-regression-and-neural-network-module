#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../modules/pmc.h"

void main(){
 mtx samples, ds, **w, output, sample, d;
 char filename[100];
 int i, j, qtspl, qtdatabyrow, *qtneurons, qtinput, qtlayers, steps, *ftype;
 long double Eqm, lastEqm, precis, lrn;
//start protocols
 precis=10.e-6;
 lrn=0.1;
 qtspl=200;
 qtlayers=2;
 ftype=(int *)malloc(qtlayers*sizeof(int));
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtinput=4;//3+1 input + bias
 qtneurons[0]=10;
 qtneurons[1]=1;
 qtdatabyrow=qtinput+qtneurons[1];
 qtdatabyrow=qtdatabyrow-1; //bias is not included in file
 samples=mtxload("Table.dat", qtspl, qtdatabyrow);
 addcol(&samples, crystalmatrix(qtspl, 1, -1.0), 0); //including bias
 qtdatabyrow=qtdatabyrow+1; //now it is
 w=createneurons(qtneurons, qtlayers, qtinput, 1);
 ds=mtxcut(samples,0, qtspl, qtinput, qtdatabyrow-qtinput);
 samples=mtxcut(samples, 0, qtspl, 0, qtinput);

//will repeat the adjust process until reachs the precision value
 ftype[0]=SIGM;
 ftype[1]=LINEAR;
 d=nullmatrix(1,1);
 sample=nullmatrix(1,1);
 steps=0;
 lastEqm=0.0;
 Eqm=meansqrerr(samples, w, qtneurons, qtlayers, ds, ftype);
 while(fabsl(Eqm-lastEqm)>precis){
  lastEqm=Eqm;
  for(i=0; i<samples.nrows; i++){
   mtxcopy(&sample,mtxcut(samples, i, 1, 0, qtinput));
   mtxcopy(&d,mtxcut(ds, i, 1, 0, ds.ncols));
   adjust(sample, w, qtneurons, qtlayers, d, lrn, ftype);
  }
  Eqm=meansqrerr(samples, w, qtneurons, qtlayers, ds, ftype);
  steps++;
 }
printf("number of steps: %d\n", steps);
savenet("layers", "layer",w, qtneurons, qtlayers);
}

