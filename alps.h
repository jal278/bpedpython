#ifndef ALPS_H
#define ALPS_H

#include "experiments.h"
#include "population.h"
#include "population_state.h"
#include "graph.h"
#define ALPS_MAX_LAYERS 10
class alps {

 public:
  long evals;
  int generation;
  int scale;
  int max_layers;
  int clayers;
  int layer_ceiling;
  bool novelty;
  int popsize;
  plot ageplot;
  noveltyarchive *protoarc;
  Genome* start_genome;
  population_state* layers[ALPS_MAX_LAYERS];
  evaluatorfunc evalf;

  vector< vector< float > > age_record;

  alps(int _layers, int _scale,Genome* sg,population_state* initial,long maxevals=1000000) {
   start_genome =sg->duplicate(0);
   scale=_scale;
   max_layers=_layers;
   clayers=1;
   generation=0;
   novelty=initial->novelty;
   protoarc=new noveltyarchive(initial->archive);
   evalf=initial->pop->evaluator;
   layers[0]=initial;
   layers[0]->max_age=scale;
   update_layer_cnt();
   evals=0;
   popsize=initial->pop->organisms.size();
   vector<float> k;
   age_record.push_back(k);

  }

 void update_layer_cnt() {
   layer_ceiling = clayers*clayers*scale;
  }

 void reproduce_alps(population_state* r,population_state* bef,population_state* aft,int psize,bool onlybefore=false) {
  int added=0;

  while(added<psize) {  
   Organism* baby;

   if(!onlybefore &&(bef==NULL || randfloat()<0.5))
	baby=r->pop->species[0]->reproduce_simple(0,r->pop);
   else {
        baby=bef->pop->species[0]->reproduce_simple(0,r->pop);    
   }
   evals++;
   if(baby->age < r->max_age) {
    added++;
    r->pop->species[0]->add_Organism(baby);
    baby->species=r->pop->species[0];  //Point baby to its species
   }
   else
    continue;

  }
  r->pop->rebuild();
 }

 void reproduce_layer(int i) {
  population_state *before=NULL, *after=NULL;
  if(i>0)
	before=layers[i-1];
  if(i<(clayers-1))
	after=layers[i+1];

  reproduce_alps(layers[i],before,after,popsize);
 }

population_state* create_new_layer(int from_layer) {

 int psize=layers[clayers-1]->pop->organisms.size();

 noveltyarchive *new_arc =new noveltyarchive(protoarc);
 Population* new_pop=new Population(start_genome,1);
 new_pop->set_evaluator(evalf);

 population_state* np = new population_state(new_pop,novelty,new_arc);
 np->max_age=layer_ceiling*2;
 reproduce_alps(np,layers[clayers-1],NULL,psize,true);
 np->pop->evaluate_all();
 cout << "POPSIZE: " << np->pop->organisms.size();
 cout << "NEWLAYERAGE:" << endl;
 np->pop->print_avg_age();

 evals+=np->pop->organisms.size();
 return np;

}

void repopulate(population_state* r) {

  int psize = r->pop->organisms.size();
  int added=0;

  while(added<psize) {  
   		Organism* baby;

		Genome* new_genome=start_genome->duplicate(added); 		
		//new_genome->mutate_link_weights(1.0,1.0,GAUSSIAN);
		new_genome->mutate_link_weights(2.0,1.0,COLDGAUSSIAN);
		new_genome->mutate_node_parameters(0.01,1.0,3.0,1.0,true);
		//new_genome->mutate_node_parameters(3.0,1.0,4.0,1.0,true);
		new_genome->randomize_traits();
		baby=new Organism(0.0,new_genome,1);
    r->pop->species[0]->add_Organism(baby);
    baby->species=r->pop->species[0];  //Point baby to its species
    added++;
  }
  r->pop->rebuild();
  r->pop->evaluate_all();
  evals+=r->pop->organisms.size();
} 


void do_alps() {
  
 while(true) {
  cout << "alps generation " << generation << endl;
  cout << "layers " << clayers << " ceiling" << layer_ceiling << endl; 
 //reproduce layers
  for(int i=0;i<clayers;i++) {
   maze_generational_epoch(layers[i],generation);
   age_record[i].push_back(layers[i]->pop->avgage);
   reproduce_layer(i);
  }
  ageplot.plot_data(age_record);
  //if time to add new layer, create if from previous layer
  if(clayers<max_layers && ((generation+1)==layer_ceiling)) {
    cout << "Adding layer " << clayers << endl;
    //create new layer, set age limit, etc
    layers[clayers]=create_new_layer(clayers);
    clayers++;
    vector<float> k;
    age_record.push_back(k);
    update_layer_cnt();
   
  } 
  
  //if multiple of scale, recreate first layer from start genome
  if((generation+1)%scale==0) {
	cout << "regenerating first layer" << endl;
	repopulate(layers[0]);
        layers[0]->clear_measure_pop();	
  }

  generation++;
 }
 }
};
#endif
