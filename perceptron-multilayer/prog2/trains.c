#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../modules/pmc.h"

#include <unistd.h> 

void main(){
 mtx samples, ds, output, sample, d, copy;
 pmcnet net, oldnet1, oldnet2;
 char filename[100];
 int i, j, qtspl, qtdatabyrow, *qtneurons, qtw1, qtlayers, steps, *ftype;
 long double Eqm, lastEqm, precis, lrn, momentum;
//start protocols
 precis=10.e-6;
 lrn=0.1;
 momentum=0.9;
 qtspl=120;
 qtlayers=2;
 ftype=(int *)malloc(qtlayers*sizeof(int));
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtw1=5;//4+1 input + bias
 qtneurons[0]=15;
 qtneurons[1]=3;
 qtdatabyrow=qtw1+qtneurons[1];
 samples=mtxload("Table.dat", qtspl, qtdatabyrow);
 srand(2);
 net=pmccreatenet(qtneurons, qtlayers, qtw1);
 ds=mtxcut(samples,0, qtspl, qtw1, qtdatabyrow-qtw1);
 samples=mtxcut(samples, 0, qtspl, 0, qtw1);

 oldnet1=clonenet(net);
 oldnet2=clonenet(net);

//will repeat the adjust process until reachs the precision value
 ftype[0]=SIGM;
 ftype[1]=SIGM;
 d=nullmatrix(1,1);
 sample=nullmatrix(1,1);
 steps=0;
 lastEqm=0.0;
 Eqm=meansqrerr(samples, net, ds, ftype);
 while(fabsl(Eqm-lastEqm)>precis){
  lastEqm=Eqm;
  for(i=0; i<samples.nrows; i++){
   copy=mtxcut(samples, i, 1, 0, qtw1);
   mtxcopy(&sample,copy);
   mtxfree(&copy);
   copy=mtxcut(ds, i, 1, 0, ds.ncols);
   mtxcopy(&d, copy);
   mtxfree(&copy);
   adjustbymomentum(&net, &oldnet1, &oldnet2, momentum);
   adjust(sample, net, d, lrn, ftype);
  }
 sleep(1);
  Eqm=meansqrerr(samples, net, ds, ftype);
  steps++;
 }
printf("number of steps: %d\n", steps);
savenet("layers", "layer", net);
}

