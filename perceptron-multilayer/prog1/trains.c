#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../modules/pmc.h"

void main(){
 mtx samples, ds, output, sample, d, copy;
 pmcnet net, oldnet;
 char filename[100];
 int i, j, qtspl, qtdatabyrow, *qtneurons, qtw1, qtlayers, steps, *ftype;
 long double Eqm, lastEqm, precis, lrn, momentum;
//start protocols
 precis=10.e-6;
 lrn=0.1;
 momentum=0.9;
 qtspl=200;
 qtlayers=2;
 ftype=(int *)malloc(qtlayers*sizeof(int));
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtw1=4;//3+1 input + bias
 qtneurons[0]=10;
 qtneurons[1]=1;
 qtdatabyrow=qtw1+qtneurons[1];
 qtdatabyrow=qtdatabyrow-1; //bias is not included in the file
 samples=mtxload("Table.dat", qtspl, qtdatabyrow);
 copy=crystalmatrix(qtspl, 1, -1.0);
 addcol(&samples, copy, 0); //including bias
 mtxfree(&copy);
 qtdatabyrow=qtdatabyrow+1; //now it is
 ftype[0]=SIGM;
 ftype[1]=LINEAR;
 srand(1);
 net=pmccreatenet(qtneurons, qtlayers, qtw1, ftype);
 ds=mtxcut(samples,0, qtspl, qtw1, qtdatabyrow-qtw1);
 samples=mtxcut(samples, 0, qtspl, 0, qtw1);
 oldnet=clonenet(net);

//will repeat the adjust process until reachs the precision value
 d=nullmatrix(1,1);
 sample=nullmatrix(1,1);
 steps=0;
 lastEqm=0.0;
 Eqm=meansqrerr(samples, net, ds);
 while(fabsl(Eqm-lastEqm)>precis){
  lastEqm=Eqm;
  for(i=0; i<samples.nrows; i++){
   copy=mtxcut(samples, i, 1, 0, qtw1);
   mtxcopy(&sample,copy);
   mtxfree(&copy);
   copy=mtxcut(ds, i, 1, 0, ds.ncols);
   mtxcopy(&d,copy);
   mtxfree(&copy);
   adjustbymomentum(&net, &oldnet, momentum);
   adjust(sample, net, d, lrn);
  }
  Eqm=meansqrerr(samples, net, ds);
  steps++;
 }
printf("number of steps: %d\n", steps);
//savenet("layers", "layer",net);
pmcsavenet("layers", net);

}

