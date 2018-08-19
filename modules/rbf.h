#include "pmc.h"

typedef struct
{
 long double *o2;
 pmclayer layer1;
 pmcnet taillayers;
}rbfnet;

pmclayer starthiddenlayer(mtx samples, int qtneurons);
long double distance(mtx vec1, mtx vec2);
void vectorquantization(pmclayer *layer, mtx samples, int *class);
long double* createhiddenlayer(pmclayer *layer, mtx samples, int qtneurons);
rbfnet rbfcreatenet(mtx samples, int qtneurons1, int *qtneuronstail, int qtlayerstail);
long double rbfneuronl1think(mtx x, mtx w, long double o2);
mtx rbfl1think(mtx x, pmclayer layer, long double *o2);
void rbfsmplsprocess(mtx *samples, rbfnet net);

void rbfsavenet(char *dirname, rbfnet net);
