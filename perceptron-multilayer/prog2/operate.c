#include <stdio.h>
#include <stdlib.h>
#include "../../modules/pmc.h"

void main(){
 mtx input, y;
 pmcnet net;
 int *qtneurons, qtlayers, qtw1, *ftype;
//quantity of input including the threshold
 qtw1=5;
 qtlayers=2;
 ftype=(int *)malloc(qtlayers*sizeof(int));
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtneurons[0]=15;
 qtneurons[1]=3;
//load a network
 net=loadnet("layers", "layer", qtneurons, qtlayers, qtw1);
 ftype[0]=SIGM;
 ftype[1]=SIGM;

 input=nullmatrix(1,qtw1);
 input.data[0][0]=-1.0;
 for(int i=1; i<input.ncols; i++){
  scanf("%Lf", &input.data[0][i]);
 }

 y=netansw(input, net, ftype);
 draw(y);
}
