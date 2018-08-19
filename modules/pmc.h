#include "matrix.h"

typedef struct
{
 int qtneurons;
 mtx *w;
}pmclayer;

typedef struct
{
 int qtlayers;
 int *ftype;
 pmclayer *layer ;
}pmcnet;

#define SIGM 0
#define TGHYP 1
#define LINEAR 2

long double sigm(long double u);
long double sigm2(long double u);
long double func(long double u, int ftype);
long double func2(long double u, int ftype);
pmclayer pmccreatelayer(int qtneurons, int qtw);
pmcnet pmccreatenet(int *qtneurons, int qtlayers, int qtw1, int *ftype);
long double neuronthink(mtx x, mtx w);
long double neuronansw(long double input, int ftype);
mtx layerthink(mtx x, pmclayer layer);
mtx layeransw(mtx input, int ftype);
mtx netthink(mtx x, pmcnet net);
mtx netansw(mtx x, pmcnet net);
void fitbydelta(mtx x, pmcnet net, mtx d, mtx *input, mtx *y, long double lrn);
void adjust(mtx x, pmcnet net, mtx d, long double lrn);
void adjustbymomentum(pmcnet *net, pmcnet *oldnet, long double momentum);
long double meansqrerr(mtx samples, pmcnet net, mtx d);
void savelayer(char *filename, pmclayer layer);
void savenet(char *dirname, char *fileprefix, pmcnet net);
void pmcsavenet(char *dirname, pmcnet net);
pmclayer loadlayer(char *filename, int qtneurons, int qtw);
pmcnet loadnet(char *dirname, char *fileprefix, int *qtneurons, int qtlayers, int qtw1);
pmcnet pmcloadnet(char *dirname);
pmclayer clonelayer(pmclayer layer);
pmcnet clonenet(pmcnet net);
void pmclayerfree(pmclayer *layer);
void pmclayercopy(pmclayer *layer, pmclayer copy);