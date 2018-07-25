#include "matrix.h"

typedef struct
{
 int qtneurons;
 mtx *w;
}pmclayer;

typedef struct
{
 int qtlayers;
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
pmcnet pmccreatenet(int *qtneurons, int qtlayers, int qtw1);
long double neuronthink(mtx x, mtx w);
long double neuronansw(long double input, int ftype);
mtx layerthink(mtx x, pmclayer layer);
mtx layeransw(mtx input, int ftype);
mtx netthink(mtx x, pmcnet net, int *ftype);
mtx netansw(mtx x, pmcnet net, int *ftype);
void fitbydelta(mtx x, pmcnet net, mtx d, mtx *input, mtx *y, long double lrn, int *ftype);
void adjust(mtx x, pmcnet net, mtx d, long double lrn, int *ftype);
void adjustbymomentum(pmcnet *net, pmcnet *oldnet1, pmcnet *oldnet2, long double momentum);
long double meansqrerr(mtx samples, pmcnet net, mtx d, int *ftype);
void savelayer(char *filename, pmclayer layer);
void savenet(char *dirname, char *fileprefix, pmcnet net);
pmclayer loadlayer(char *filename, int qtneurons, int qtw);
pmcnet loadnet(char *dirname, char *fileprefix, int *qtneurons, int qtlayers, int qtw1);
pmclayer clonelayer(pmclayer layer);
pmcnet clonenet(pmcnet net);

