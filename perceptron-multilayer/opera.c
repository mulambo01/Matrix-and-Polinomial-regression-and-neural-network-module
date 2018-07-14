#include <stdio.h>
#include <stdlib.h>
#include "../modules/pmc.h"

void main(){
 mtx **w, input, y;
 int *qtneurons, qtlayers, qtinput;
 qtinput=4;
 qtlayers=2;
 qtneurons=(int *)malloc(qtlayers*sizeof(int));
 qtneurons[0]=15;
 qtneurons[1]=3;
w=loadnet("layers", "layer", qtneurons, qtlayers, qtinput);
 input=nullmatrix(1,5);
 input.data[0][0]=-1.0;
 for(int i=1; i<input.ncols; i++){
  scanf("%Lf", &input.data[0][i]);
 }
 y=netansw(input, w, qtneurons, qtlayers);
draw(y);
}
