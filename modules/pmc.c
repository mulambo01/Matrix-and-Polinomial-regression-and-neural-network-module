#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "matrix.h"


long double sigm(long double u){
 return 1.0/(1.0+expl(-u));
}
long double sigm2(long double u){
 return (1.0-powl(sigm(u),2));
}

long double func(long double u){
 return sigm(u);
}

long double func2(long double u){
 return sigm2(u);
}

mtx **createneurons(int *qtneurons, int qtlayers, int qtinput){
 int j, qtw;
 mtx **net=(mtx **)malloc(qtlayers*sizeof(mtx *));
 qtw=qtinput;
 for(int i=0; i<qtlayers; i++){
  net[i]=(mtx *)malloc(qtneurons[i]*sizeof(mtx));
  for(j=0; j<qtneurons[i]; j++){
   net[i][j]=randmatrix(1,qtw);
  }
  qtw=qtneurons[i]+1;//+1 pcausa do -1
 }
 return net;
}


long double neuronthink(mtx x, mtx w){
 long double input, y;
 input=vectprod(x,w);
 return input;
}
long double neuronansw(long double input){
 return func(input);
}

mtx layerthink(mtx x, mtx *w, int qtneurons){
 mtx input;
 input=nullmatrix(1, qtneurons);
 for(int i=0; i<qtneurons; i++){
  input.data[0][i]=neuronthink(x, w[i]);
 }
 return input;
}

mtx layeransw(mtx input){
 mtx y;
 y=nullmatrix(1, input.ncols);
 for(int i=0; i<input.ncols; i++){
  y.data[0][i]=neuronansw(input.data[0][i]);
 }
 return y;
}

mtx netthink(mtx x, mtx **w, int *qtneurons, int qtlayers){
 mtx y, input, x0;
 x0=crystalmatrix(1,1,-1.0);
 input=mtxclone(x);
 y=nullmatrix(1,1);
 for(int i=0; i<qtlayers; i++){
  mtxcopy(&input, layerthink(input, w[i], qtneurons[i]));
  mtxcopy(&y, layeransw(input));
  mtxcopy(&input, y);
  addcol(&input, x0, 0);
 }
 return y;
}

mtx netansw(mtx x, mtx **w, int *qtneurons, int qtlayers){
 mtx answ, y;
 y=netthink(x, w, qtneurons, qtlayers);
 answ=mtxclone(y);
 for(int i=0; i<answ.ncols; i++){
  if(answ.data[0][i]>0.5){
   answ.data[0][i]=1.0;
  }
  else{
   answ.data[0][i]=0.0;
  }
 }
 return answ;
}


void fitbydelta(mtx **w, mtx *input, mtx x, mtx *y, mtx d, int qtlayers, long double lrn){
 mtx *delta, change;
 int qtneurons, qtnfrwrd, i, j, layer;
 long double der, del;
 delta=(mtx *)malloc(qtlayers*sizeof(mtx));
//last layer
 layer=qtlayers-1;
 qtneurons=input[layer].ncols;
 change=nullmatrix(1,1);
 delta[layer]=nullmatrix(1,qtneurons);
 for(i=0; i<qtneurons; i++){
  der=func2(input[layer].data[0][i]);
  del=(d.data[0][i]-y[layer].data[0][i])*der;
  delta[layer].data[0][i]=del;
  mtxcopy(&change,mtxmult(y[layer-1], del*lrn));
  mtxcopy(&w[layer][i],mtxsum(w[layer][i],change));
 }
//middle layers
 for(layer--;layer>0;layer--){
  qtnfrwrd=qtneurons;
  qtneurons=input[layer].ncols;
  delta[layer]=nullmatrix(1,qtneurons);
  for(i=0; i<qtneurons; i++){
   der=func2(input[layer].data[0][i]);
   del=0.0;
   for(j=0; j<qtnfrwrd; j++){
    del=del+delta[layer+1].data[0][j]*w[layer+1][j].data[0][i];
   }
   del=-del*der;
   delta[layer].data[0][i]=del;
   mtxcopy(&change, mtxmult(y[layer-1], del*lrn));
   mtxcopy(&w[layer][i], mtxsum(w[layer][i],change));
  }
 }
//first layer
 qtnfrwrd=qtneurons;
 qtneurons=input[0].ncols;
 delta[0]=nullmatrix(1,qtneurons);
 for(i=0; i<qtneurons; i++){
  der=func2(input[0].data[0][i]);
  del=0.0;
  for(j=0; j<qtnfrwrd; j++){
   del=del+delta[1].data[0][j]*w[1][j].data[0][i];
  }
  del=-del*der;
  delta[0].data[0][i]=del;
  mtxcopy(&change, mtxmult(x, del*lrn));
  mtxcopy(&w[0][i], mtxsum(w[0][i],change));
 }
//end
}

mtx adjust(mtx x, mtx d, mtx **w, int *qtneurons, int qtlayers, long double lrn){
 mtx *y, *input, yy, ii, x0;
 x0=crystalmatrix(1,1,-1.0);
 ii=mtxclone(x);
 yy=nullmatrix(1,1);
 y=(mtx *)malloc(qtlayers*sizeof(mtx));
 input=(mtx *)malloc(qtlayers*sizeof(mtx));

 for(int i=0; i<qtlayers; i++){
  mtxcopy(&ii, layerthink(ii, w[i], qtneurons[i]));
  mtxcopy(&yy, layeransw(ii));
  input[i]=mtxclone(ii);
  if(i!=qtlayers-1){
   addcol(&yy, x0, 0);
  }
  y[i]=mtxclone(yy);
  mtxcopy(&ii, yy);
 }
 fitbydelta(w, input, x, y, d, qtlayers, lrn);
}

long double meansqrerr(mtx samples, mtx d, mtx **w, int *qtneurons, int qtlayers){
 long double result=0.0;
 int qtspl=samples.nrows;
 int qtinput=samples.ncols;
 mtx u, x, output, error;
 u=nullmatrix(d.nrows, d.ncols);
 x=nullmatrix(1, qtinput);
 output=nullmatrix(1, d.ncols);
 for(int i=0; i<qtspl; i++){
  mtxcopy(&x,mtxcut(samples, i, 1, 0, qtinput));
  mtxcopy(&output, netthink(x, w, qtneurons, qtlayers));
  putline(&u, output, i);
 }
 mtxcopy(&u, mtxmult(u, -1.0));
 error=mtxsum(u, d);
 result=vectprod(error, error)/qtspl;

 return result;
}

void savelayer(char *filename, mtx *w, int qtneurons){
 mtx wbylayer;
 wbylayer=nullmatrix(qtneurons, w[0].ncols);
 for(int i=0; i<qtneurons; i++){
  putline(&wbylayer, w[i], i);
 }
 mtxsave(filename, wbylayer);
}

void savenet(char *dirname, char *fileprefix, mtx **w, int *qtneurons, int qtlayers){
 int namelen, sufixlen;
 sufixlen=1+(int)log10(qtlayers-1);
 namelen=strlen(dirname)+strlen(fileprefix)+sufixlen+10;
 char filename[namelen], sufix[sufixlen];
 mkdir(dirname, S_IRWXU);
 for(int i=0; i<qtlayers; i++){
  sprintf(sufix, "%d", (i+(int)pow(10.0,sufixlen)));
  sufix[0]='-';
  sprintf(filename, "%s/%s%s.dat",dirname, fileprefix, sufix);
  savelayer(filename, w[i], qtneurons[i]);
 }
}

mtx* loadlayer(char *filename, int qtneurons, int qtw){
 mtx layer, *w;
 w=(mtx *)malloc(qtneurons*sizeof(mtx));
 layer=mtxload(filename, qtneurons, qtw);
 for(int i=0; i<qtneurons; i++){
  w[i]=mtxcut(layer, i, 1, 0, qtw);
 } 
 return w;
}

mtx** loadnet(char *dirname, char *fileprefix, int *qtneurons, int qtlayers, int qtinput){
 mtx **w;
 int namelen, sufixlen, qtw;
 sufixlen=1+(int)log10(qtlayers-1);
 namelen=strlen(dirname)+strlen(fileprefix)+sufixlen+10;
 char filename[namelen], sufix[sufixlen];
 w=(mtx **)malloc(qtlayers*sizeof(mtx *));
 qtw=qtinput+1; //to include the limiar
 for(int i=0; i<qtlayers; i++){
  sprintf(sufix, "%d", (i+(int)pow(10.0,sufixlen)));
  sufix[0]='-';
  sprintf(filename, "%s/%s%s.dat",dirname, fileprefix, sufix);
  w[i]=loadlayer(filename, qtneurons[i], qtw);
  qtw=qtneurons[i]+1;
 }
 return w;
}

