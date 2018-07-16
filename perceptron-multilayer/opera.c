#include <stdio.h>
#include <stdlib.h>
#include "../modules/pmc.h"

void main(){
 mtx **w, input, y;
 int *qtneurons, qtlayers, qtinput, ftype;
 ftype=0;
//quantity of input including the threshold
 qtinput=5;
 qtlayers=2;
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtneurons[0]=15;
 qtneurons[1]=3;
//load a network
//the quantity of inputs dont include the threshold
 w=loadnet("layers", "layer", qtneurons, qtlayers, qtinput-1);

 input=nullmatrix(1,qtinput);
 input.data[0][0]=-1.0;
 for(int i=1; i<input.ncols; i++){
  scanf("%Lf", &input.data[0][i]);
 }

 mtx d;
 d=nullmatrix(1,qtneurons[qtlayers-1]);

 y=netansw(input, w, qtneurons, qtlayers, ftype);
 draw(y);
}
