#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "matrix.h"

//types of functions
//sigmoid function
#define SIGM 0
//hyberbolic tangent
#define TGHYP 1
//linear function
#define LINEAR 2

//return the sigmoid function
long double sigm(long double u){
 return 1.0/(1.0+expl(-u));
}
//return the derivative of sigmoid function
long double sigm1(long double u){
 return (1.0-powl(sigm(u),2));
}
long double tghyp(long double u){
 return (tanh(u));
}

long double tghyp1(long double u){
 return (1.0/(powl(cosh(u),2)));
}

long double linear(long double u){
 return (u);
}

long double linear1(long double u){
 return (1.0);
}

//return the value of function of type "ftype"
long double func(long double u, int ftype){
 if(ftype==SIGM)return sigm(u);
 else if(ftype==TGHYP)return tghyp(u);
}
//return the value of derivative function of type "ftype"
long double func1(long double u, int ftype){
 if(ftype==SIGM)return sigm1(u);
 else if(ftype==TGHYP)return tghyp1(u);
}

//create a neural network
//qtneurons is an array with the quantity of neurons by layer
//qtlayer defines the quantity of layers
//qtinput define the quantity of data input without the threshold
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

//return the linear combination with values of the synaptic weights and inputs
long double neuronthink(mtx x, mtx w){
 long double input, y;
 input=vectprod(x,w);
 return input;
}
//return the output of a neuron
//input is the answer of function neuronthink
long double neuronansw(long double input, int ftype){
 return func(input, ftype);
}
//return the linear combinations of each neuron with the sample x
mtx layerthink(mtx x, mtx *w, int qtneurons){
 mtx input;
 input=nullmatrix(1, qtneurons);
 for(int i=0; i<qtneurons; i++){
  input.data[0][i]=neuronthink(x, w[i]);
 }
 return input;
}
//return the output of a neural layer
//input is the answer of function layerthink
mtx layeransw(mtx input, int ftype){
 mtx y;
 y=nullmatrix(1, input.ncols);
 for(int i=0; i<input.ncols; i++){
  y.data[0][i]=neuronansw(input.data[0][i], ftype);
 }
 return y;
}

//return the answer of a neural network
//array w stores the synaptic weights matrix indexing by the layer and position of neuron
//array qtneurons stores the quantity of neurons in each layer
//qtlayers is the quantity of layers
//ftype is the type of activation function
mtx netthink(mtx x, mtx **w, int *qtneurons, int qtlayers, int *ftype){
 mtx y, input, x0;
 x0=crystalmatrix(1,1,-1.0);
 input=mtxclone(x);
 y=nullmatrix(1,1);
 for(int i=0; i<qtlayers; i++){
  mtxcopy(&input, layerthink(input, w[i], qtneurons[i]));
  mtxcopy(&y, layeransw(input, ftype[i]));
  mtxcopy(&input, y);
  addcol(&input, x0, 0);
 }
 return y;
}

//just call the function netthink and process the output returning the sample class
mtx netansw(mtx x, mtx **w, int *qtneurons, int qtlayers, int *ftype){
 mtx answ, y;
 y=netthink(x, w, qtneurons, qtlayers, ftype);
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

//it will adjust each weight of a network using the inputs and outputs
void fitbydelta(mtx x, mtx **w, int qtlayers, mtx d, mtx *input, mtx *y, long double lrn, int *ftype){
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
  der=func1(input[layer].data[0][i], ftype[layer]);
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
   der=func1(input[layer].data[0][i], ftype[layer]);
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
 qtneurons=input[layer].ncols;
 delta[layer]=nullmatrix(1,qtneurons);
 for(i=0; i<qtneurons; i++){
  der=func1(input[layer].data[0][i], ftype[layer]);
  del=0.0;
  for(j=0; j<qtnfrwrd; j++){
   del=del+delta[layer+1].data[0][j]*w[layer+1][j].data[0][i];
  }
  del=-del*der;
  delta[layer].data[0][i]=del;
  mtxcopy(&change, mtxmult(x, del*lrn));
  mtxcopy(&w[layer][i], mtxsum(w[layer][i],change));
 }
//end
}

//it will genearate the arrays input and y and pass them to the function fitbydelta
mtx adjust(mtx x, mtx **w, int *qtneurons, int qtlayers, mtx d, long double lrn, int *ftype){
 mtx *y, *input, yy, ii, x0;
 x0=crystalmatrix(1,1,-1.0);
 ii=mtxclone(x);
 yy=nullmatrix(1,1);
 y=(mtx *)malloc(qtlayers*sizeof(mtx));
 input=(mtx *)malloc(qtlayers*sizeof(mtx));
 for(int i=0; i<qtlayers; i++){
  mtxcopy(&ii, layerthink(ii, w[i], qtneurons[i]));
  mtxcopy(&yy, layeransw(ii, ftype[i]));
  input[i]=mtxclone(ii);
  if(i!=qtlayers-1){
   addcol(&yy, x0, 0);
  }
  y[i]=mtxclone(yy);
  mtxcopy(&ii, yy);
 }
 fitbydelta(x, w, qtlayers, d, input, y, lrn, ftype);
}

//calculate the mean square error of the samples
long double meansqrerr(mtx samples, mtx **w, int *qtneurons, int qtlayers, mtx d, int *ftype){
 long double result=0.0;
 int qtspl=samples.nrows;
 int qtinput=samples.ncols;
 mtx u, x, output, error;
 u=nullmatrix(d.nrows, d.ncols);
 x=nullmatrix(1, qtinput);
 output=nullmatrix(1, d.ncols);
 for(int i=0; i<qtspl; i++){
  mtxcopy(&x,mtxcut(samples, i, 1, 0, qtinput));
  mtxcopy(&output, netthink(x, w, qtneurons, qtlayers, ftype));
  putline(&u, output, i);
 }
 mtxcopy(&u, mtxmult(u, -1.0));
 error=mtxsum(u, d);
 result=vectprod(error, error)/qtspl;

 return result;
}

//save a layer in a text file
//receive the array of weights by neuron, put it in a matrix and save it
void savelayer(char *filename, mtx *w, int qtneurons){
 mtx wbylayer;
 wbylayer=nullmatrix(qtneurons, w[0].ncols);
 for(int i=0; i<qtneurons; i++){
  putline(&wbylayer, w[i], i);
 }
 mtxsave(filename, wbylayer);
}

//save a network in a text file
//dirname is the name of directory
//fileprefix is the the prefix that each layer file will receive
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

//load a neural layer by a text file
mtx* loadlayer(char *filename, int qtneurons, int qtw){
 mtx layer, *w;
 w=(mtx *)malloc(qtneurons*sizeof(mtx));
 layer=mtxload(filename, qtneurons, qtw);
 for(int i=0; i<qtneurons; i++){
  w[i]=mtxcut(layer, i, 1, 0, qtw);
 } 
 return w;
}

//load a neural network by the files of a directory with a prefix
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

