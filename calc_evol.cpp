#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;

#define POPSIZE 250 //TODO: 500 for biped 
#define TRIALS 300
#define SAMPLES_PER 20
#define DELTA 0.01
double coords[POPSIZE][TRIALS][2];

double dist(double x[2],double y[2]) {
 double dx= x[0]-y[0];
 double dy= x[1]-y[1];
 return dx*dx+dy*dy;
}

int distinct(double *pts,int len) {
double points[POPSIZE*2];
int size=0;

for(int x=0;x<len;x++) {
 bool pass=true;
 for(int y=0;y<size;y++) 
  if(dist(&points[y*2],&pts[x*2])<DELTA) {
   pass=false;
   break;
  }
  
  if(pass) {
   points[size*2]=pts[x*2];
   points[size*2+1]=pts[x*2+1];
   size++;
  }
 }
 return size;
}


double test_indiv(double *pts,int len) {
  double sum=0.0;

  for(int x=0;x<SAMPLES_PER;x++) {
   for(int y=0;y<SAMPLES_PER;y++) {
   double samp[2];
   samp[0]=((float)x/(float)SAMPLES_PER);
   samp[1]=((float)y/(float)SAMPLES_PER);
   double mindist = dist(&pts[0],samp);
   for(int y=1;y<len;y++) {
    double newdist = dist(&pts[y*2],samp);
    if(newdist<mindist) mindist=newdist;
   }
   sum+=mindist;
  }
  }
 return sum;
}
/*
double test_all() {
 double sum=test_indiv(coords[0]);
 double min=sum;
 for(int x=1;x<POPSIZE;x++) {
  double newval=test_indiv(coords[x]);
  if(newval<min) min=newval;
  sum+=newval;
 }
 cout << sum << endl;
 cout << min << endl;
 return min;
}
*/
void read_file(char* fn) {
 ifstream file(fn);
 int cnt=-1;
 int coord_cnt=0;
 int xy_cnt=0;
 while(!file.eof()) {
  char line[80];
  file >> line;
  if (strcmp(line,"---")==0) {
   cnt+=1;
   coord_cnt=0;
  }
  else { 
   coords[cnt][coord_cnt][xy_cnt]=atof(line);
   xy_cnt+=1; 
   if(xy_cnt==2) {
    xy_cnt=0;
    coord_cnt+=1;
   }
  }
 }
}

/*
int main(int argc, char** argv) {
char fn[100]="evolnov__evolvability1.dat";
srand(time(NULL));
read_file(argv[1]);
test_all();
return 0;
}
*/