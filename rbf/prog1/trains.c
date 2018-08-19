#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../modules/rbf.h"

void main(){
 rbfnet net;
 mtx samples, sample, ds, d, copy;
 int i, qtneurons1, qtsamples, qtw1, qtneurons2[1], qtdatabyrow, steps, ftype[1], qttaillayers;
 long double *o2, Eqm, lastEqm, lrn, precis;

 lrn=0.1;
 precis=10.0e-30;
 qtsamples=150;
 qtneurons1=2;
 qtw1=3;
 qttaillayers=1;
 qtneurons2[0]=1;

//selecting and straightening samples
 qtdatabyrow=qtw1+qtneurons2[0];
 samples=mtxload("Table.dat", qtsamples, qtdatabyrow);
 ds=mtxcut(samples, 0, qtsamples, qtw1, qtneurons2[0]);
 copy=mtxcut(samples, 0, qtsamples, 0, qtw1);
 mtxcopy(&samples, copy);
 mtxfree(&copy);

//generating rbf first layer and updating samples
 net=rbfcreatenet(samples, qtneurons1, qtneurons2, qttaillayers);
 rbfsmplsprocess(&samples, net);

//pmc standard adjust
 ftype[0]=LINEAR;
 d=nullmatrix(1,1);
 sample=nullmatrix(1,1);
 steps=0;
 lastEqm=0.0;
 Eqm=meansqrerr(samples, net.taillayers, ds, ftype);
 while(fabsl(Eqm-lastEqm)>precis){
  lastEqm=Eqm;
  for(i=0; i<samples.nrows; i++){
   copy=mtxcut(samples, i, 1, 0, samples.ncols);
   mtxcopy(&sample,copy);
   mtxfree(&copy);
   copy=mtxcut(ds, i, 1, 0, ds.ncols);
   mtxcopy(&d, copy);
   mtxfree(&copy);
   adjust(sample, net.taillayers, d, lrn, ftype);
  }
  Eqm=meansqrerr(samples, net.taillayers, ds, ftype);
  steps++;
 }
 printf("number of steps: %d\n", steps);
 rbfsavenet("layers", net);
 
 free(o2);
}
