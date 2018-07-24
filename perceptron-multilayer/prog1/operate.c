#include <stdio.h>
#include <stdlib.h>
#include "../../modules/pmc.h"

void main(){
 mtx **w, input, y;
 int *qtneurons, qtlayers, qtinput, *ftype;
//quantity of input including the bias
 qtinput=4;
 qtlayers=2;
 ftype=(int *)malloc(qtlayers*sizeof(int));
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtneurons[0]=10;
 qtneurons[1]=1;
//load a network
//the quantity of inputs dont include the bias
 w=loadnet("layers", "layer", qtneurons, qtlayers, qtinput-1);
 ftype[0]=SIGM;
 ftype[1]=LINEAR;

 input=nullmatrix(1,qtinput);
 input.data[0][0]=-1.0;
 for(int i=1; i<input.ncols; i++){
  scanf("%Lf", &input.data[0][i]);
 }

 mtx d;
 d=nullmatrix(1,qtneurons[qtlayers-1]);

 y=netthink(input, w, qtneurons, qtlayers, ftype);
 draw(y);
}
