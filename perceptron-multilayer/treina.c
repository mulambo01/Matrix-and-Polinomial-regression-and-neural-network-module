
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../modules/pmc.h"


void main(){
 mtx samples, ds, **w, output, sample, d;
 char filename[100];
 int i, j, qtspl, qtdatabyrow, *qtneurons, qtinput, qtlayers, steps, *ftype;
 long double Eqm, lastEqm, precis, lrn;
//start protocols
 precis=10.e-6;
 lrn=0.1;
 qtspl=120;
 qtlayers=2;
 ftype=(int *)malloc(qtlayers*sizeof(int));
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtinput=5;//4 +1
 qtneurons[0]=15;
 qtneurons[1]=3;
 qtdatabyrow=qtinput+qtneurons[1];
 samples=mtxload("Tabela3.txt", qtspl, qtdatabyrow);
 w=createneurons(qtneurons, qtlayers, qtinput);
 ds=mtxcut(samples,0, qtspl, qtinput, qtdatabyrow-qtinput);
 samples=mtxcut(samples, 0, qtspl, 0, qtinput);

//will repeat the adjust process until reachs the precision value

 ftype[0]=SIGM;
 ftype[1]=SIGM;
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
savenet("layers", "layer",w, qtneurons, qtlayers);
}

