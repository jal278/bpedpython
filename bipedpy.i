%include "std_vector.i"
%include "carrays.i"
%array_class(float,floatArray);

%newobject *::copy;
%module bipedpy
%{
#include "biped_py.h"
%}

namespace std {
 %template(vectorf) vector<float>;
};
using namespace std;

class bipedsim {
public: 
 void make_random();

 static void seed(int sd);
 static void random_seed();
 bipedsim();

  vector<float> get_behavior();

  bipedsim* copy(bool z=false);
  int complexity();

  void map();
  static void initmaze(const char* mazefile);
  void mutate();
  bool isvalid();
  void clear();
  double distance(bipedsim* other);
  void init_rand();
  bool viable();
 
  void save(const char *fname);
  void load_new(const char *fname);

  float get_x(); 
  float get_y();
  bool solution();
  ~bipedsim();

};
