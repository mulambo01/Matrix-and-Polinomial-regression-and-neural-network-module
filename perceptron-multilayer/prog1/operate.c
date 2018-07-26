#include <stdio.h>
#include <stdlib.h>
#include "../../modules/pmc.h"

void main(){
 mtx input, y;
 pmcnet net;
 int *qtneurons, qtlayers, qtw1, *ftype;
//quantity of input including the bias
 qtw1=4;
 qtlayers=2;
 ftype=(int *)malloc(qtlayers*sizeof(int));
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtneurons[0]=10;
 qtneurons[1]=1;
//load a network
//the quantity of inputs dont include the bias
 net=loadnet("layers", "layer", qtneurons, qtlayers, qtw1);
 ftype[0]=SIGM;
 ftype[1]=LINEAR;

 input=nullmatrix(1,qtw1);
 input.data[0][0]=-1.0;
 for(int i=1; i<input.ncols; i++){
  scanf("%Lf", &input.data[0][i]);
 }

 y=netthink(input, net, ftype);
 draw(y);
}
