#ifndef POPSTATE_H
#define POPSTATE_H
class population_state {
 public:

  population_state(Population* _pop,bool _n,noveltyarchive* _arc) {
   pop=_pop;
   novelty=_n;
   archive=_arc;
   best_fitness=0;
   best_secondary = -100000.0;
  }
 
 
  double best_fitness; 
  double best_secondary; 
  int max_age; 
  noveltyarchive* archive;
  Population* pop;
  data_rec Record;
  vector<Organism*> measure_pop;
  bool novelty;

  void clear_measure_pop() {
  vector<Organism*>::iterator curorg;
 for (curorg=measure_pop.begin(); curorg!=measure_pop.end(); curorg++) delete (*curorg);
        //clear the old population
        measure_pop.clear();
  } 
};
#endif
