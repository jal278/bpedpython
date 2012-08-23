#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
using namespace std;

class plot {
public:
 FILE *gp;
 plot() {
  gp=popen("gnuplot -persist","w");
 }

 void plot_data(vector<float> x,vector<float> y) {
  fprintf(gp,"plot '-' \n");
  for(int k=0;k<x.size();k++) {
   fprintf(gp,"%f %f \n",x[k],y[k]);
  }
  fprintf(gp,"e\n");
  fflush(gp);
 }

 ~plot() {
  pclose(gp);
 }
 
};
/*
int main(int argc,char **argv) {
 plot p;
 for(int a=0;a<100;a++) {
 vector<float> x,y;
 for(int k=a;k<a+200;k++) {
   x.push_back(k);
   y.push_back(k*k);
 }
  p.plot_data(x,y);
 }
 return 0;
}
*/
