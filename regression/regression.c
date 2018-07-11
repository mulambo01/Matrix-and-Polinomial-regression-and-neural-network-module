#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../modules/linsys.h"

void main(int argc, char *argv[]){
 if(argc!=4){
  printf("Use: %s [FILE NAME] [DEGREE] [NUMBER OF POINTS]\n", argv[0]);
 }
 else{
  char *filename=argv[1];
  int degree=atoi(argv[2]);
  int nrows=atoi(argv[3]);
  mtx points=load(filename, nrows, 2);
  mtx coefs=regression(points, degree);
  draweq(coefs);
  printf("r2= %.17Le\n",r2(coefs, points, degree));
  printf("residual variance= %.17Le\n",residvariance(coefs, points, degree));
 }
}
