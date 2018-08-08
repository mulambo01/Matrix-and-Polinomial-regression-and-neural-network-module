#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pmc.h"

pmclayer starthiddenlayer(mtx samples, int qtneurons){
 pmclayer layer;
 int qtw=samples.ncols;
 mtx *w;
 w=(mtx *)malloc(qtneurons*sizeof(mtx));
 for(int i=0; i<qtneurons; i++){
  w[i]=mtxcut(samples, i, 1, 0, qtw);
 }
 layer.w=w;
 layer.qtneurons=qtneurons;
 return layer;
}

long double distance(mtx vec1, mtx vec2){
 int dimension=vec1.ncols;
 long double dist=0.0;
 for(int i=0; i<dimension; i++){
  dist=dist+powl(vec1.data[0][i]-vec2.data[0][i],2);
 }
 dist=sqrtl(dist);
 return dist;
}

void vectorquantization(pmclayer *layer, mtx samples, int *class){
 int i, j, qtsamples, qtclasses, qtw, *qtsmpbyclass;
 long double posi;
 mtx sample, copy;
 qtsamples=samples.nrows;
 qtw=samples.ncols;
 qtclasses=layer->qtneurons;
 qtsmpbyclass=(int *)calloc(qtclasses,sizeof(int));
 for(i=0; i<qtclasses; i++){
  copy=nullmatrix(1, qtw);
  mtxcopy(&(layer->w[i]), copy);
  mtxfree(&copy);
 }
 sample=nullmatrix(1,qtw);
 for(i=0; i<qtsamples; i++){
  copy=mtxcut(samples, i, 1, 0, qtw);
  mtxcopy(&sample, copy);
  mtxfree(&copy);
  for(j=0; j<qtclasses; j++){
   if(j==class[i]){
    qtsmpbyclass[j]++;
    copy=mtxsum(layer->w[j], sample);
    mtxcopy(&(layer->w[j]), copy);
    mtxfree(&copy);
    break;
   }
  }
 }
 for(j=0; j<qtclasses; j++){
  copy=mtxmult(layer->w[j], 1.0/qtsmpbyclass[j]);
  mtxcopy(&(layer->w[j]), copy);
  mtxfree(&copy);
 }

 free(qtsmpbyclass);
 mtxfree(&sample);
}

long double* createhiddenlayer(pmclayer *layer, mtx samples, int qtneurons){
 long double *o2;
 int *class, change, i, j, qtsamples, group, qtw, *qtsmpbyclass;
 mtx sample;
 o2=(long double *)calloc(qtneurons,sizeof(long double));
 qtsmpbyclass=(int *)calloc(qtneurons,sizeof(int));
 qtsamples=samples.nrows;
 qtw=samples.ncols;

 *layer=starthiddenlayer(samples, qtneurons);

 class=(int *)calloc(qtsamples,sizeof(int));
 change=1;
 while(change!=0){
  change=0;
  
  for(i=0; i<qtsamples; i++){
   sample=mtxcut(samples,i,1,0,qtw);
   group=class[i];
   for(j=0; j<qtneurons; j++){
    if(distance(sample, layer->w[j])<distance(sample, layer->w[group])){
     group=j;
     class[i]=j;
     change=1;
     break;
    }
   }
   mtxfree(&sample);
  }
  vectorquantization(layer, samples, class);
 }
 for(i=0; i<qtsamples; i++){
  sample=mtxcut(samples,i,1,0,qtw);
  group=class[i];
  qtsmpbyclass[group]++;
 
  o2[group]=o2[group]+powl(distance(sample, layer->w[group]), 2);
  mtxfree(&sample);
 }
 for(i=0; i<qtneurons; i++){
  o2[i]=o2[i]/qtsmpbyclass[i];
 }
 free(qtsmpbyclass);
 free(class);
 return o2;
}

long double rbfneuronthink(mtx x, mtx w, long double o2){
 long double answ=0.0;
 for(int i=0; i<w.ncols; i++){
  answ=answ+powl((x.data[0][i]-w.data[0][i]),2);
 }
 answ=answ/(2.0*o2);
 answ=exp(-answ);
 return answ;
}

mtx rbflayerthink(mtx x, pmclayer layer, long double *o2){
 mtx answ;
 int qtneurons=layer.qtneurons;
 long double **w=(long double **)malloc(sizeof(long double *));
 answ.nrows=1;
 answ.ncols=qtneurons;
 w[0]=(long double *)malloc(qtneurons*sizeof(long double));
 for(int i=0; i<qtneurons; i++){
  w[0][i]=rbfneuronthink(x, layer.w[i], o2[i]);
 }
 answ.data=w;
 return answ;
}


void updatesamples(mtx *samples, pmclayer layer, long double *o2){

 mtx sample, newsamples, copy;
 int nrows, qtneurons;
 nrows=samples->nrows;
 qtneurons=layer.qtneurons;
 newsamples=nullmatrix(nrows, qtneurons);

 for(int i=0; i<nrows; i++){
  sample=mtxcut(*samples, i, 1, 0, samples->ncols);
  copy=rbflayerthink(sample, layer, o2);
  putline(&newsamples, copy, i);
  mtxfree(&sample);
  mtxfree(&copy);
 }
 mtxcopy(samples, newsamples);
 mtxfree(&newsamples);
}

void main(){
 pmclayer layer1;
 mtx samples, sample, ds, copy;
 int qtneurons1, qtsamples, qtw1, qtneurons2, qtdatabyrow;
 long double *o2;
 qtneurons1=2;
 qtsamples=150;
 qtw1=3;
 qtneurons2=1;
 qtdatabyrow=qtw1+qtneurons2;
 samples=mtxload("Table.dat", qtsamples, qtdatabyrow);
 ds=mtxcut(samples, 0, qtsamples, qtw1, qtneurons2);
 copy=mtxcut(samples, 0, qtsamples, 0, qtw1);
 mtxcopy(&samples, copy);
 mtxfree(&copy);
while(1){
 o2=createhiddenlayer(&layer1, samples, qtneurons1);
pmclayerfree(&layer1);
//draw(samples);
//updatesamples(&samples, layer1, o2);
//draw(samples);
 free(o2);
}
}
