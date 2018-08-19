#include <stdio.h>
#include <stdlib.h>
#include "../../modules/pmc.h"

void main(){
 mtx input, y;
 pmcnet net;
 int qtw1;
//quantity of input including the threshold
//load a network
 net=pmcloadnet("layers");
 qtw1=net.layer[0].w[0].ncols;
 input=nullmatrix(1,qtw1);
 input.data[0][0]=-1.0;
 for(int i=1; i<input.ncols; i++){
  scanf("%Lf", &input.data[0][i]);
 }
net.ftype[1]=SIGM;
 y=netthink(input, net);
 draw(y);
}
