#include <stdio.h>
#include <stdlib.h>
#include "../../modules/pmc.h"

void main(){
 mtx input, y;
 pmcnet net;
 int qtw1, *ftype;
//quantity of input including the bias
//load a network
//the quantity of inputs dont include the bias
 net=pmcloadnet("layers");
 qtw1=net.layer[0].w[0].ncols;
 ftype=net.ftype;

 input=nullmatrix(1,qtw1);
 input.data[0][0]=-1.0;
 for(int i=1; i<input.ncols; i++){
  scanf("%Lf", &input.data[0][i]);
 }

 y=netthink(input, net);
 draw(y);
}
