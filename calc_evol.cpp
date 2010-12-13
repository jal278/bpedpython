#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;

#define POPSIZE 500 //TODO: 500 for biped, 250 for maze
#define TRIALS 300
#define SAMPLES_PER 20
#define MAXDIM 30
#define DELTA 0.01 //was 0.01

double dist(double *x,double *y,int dim) {
 double sum=0.0; 
 for (int k=0;k<dim;k++) {
  double d= x[k]-y[k];
  sum+=d*d;
 }
 return sum;
}

int distinct(double *pts,int len,int dim) {
double points[POPSIZE*MAXDIM];
int size=0;

for(int x=0;x<len;x++) {
 bool pass=true;
 for(int y=0;y<size;y++) 
  if(dist(&points[y*dim],&pts[x*dim],dim)<DELTA) {
   pass=false;
   break;
  }
  
  if(pass) {
   for(int i=0;i<dim;i++) {
   points[size*dim+i]=pts[x*dim+i];
   }
   size++;
  }
 }
 return size;
}

