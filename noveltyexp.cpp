#include "experiments.h"
#include "noveltyset.h"

#include "datarec.h"
#include "maze.h"

#include "histogram.h"

//#define DEBUG_OUTPUT 1
#include <algorithm>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>

double evolvability(Organism* org,char* fn);
using namespace std;
enum novelty_measure_type { novelty_sample, novelty_accum, novelty_sample_free };
static novelty_measure_type novelty_measure = novelty_sample;

enum fitness_measure_type { fitness_goal, fitness_drift, fitness_std,fitness_rnd };
static fitness_measure_type fitness_measure = fitness_goal;

static int number_of_samples = 1;
static int simulated_timesteps = 400;
static bool seed_mode = false;
static char seed_name[40]="";
static char mc_mazes[40]="";
static bool minimal_criteria=false;
static bool goal_attract=true;

static bool activity_stats=false;
static bool constraint_switch=false;
static bool area_of_interest=false;
static bool rand_repl=false;

void set_random_replace(bool val)
{
rand_repl = val;
}

void set_aoi(bool val)
{
area_of_interest=val;
}

void  set_constraint_switch(bool val)
{
constraint_switch=val;
}
void set_minimal_criteria(bool mc)
{
 minimal_criteria=mc;
}

void set_goal_attract(bool ga)
{
 goal_attract=ga;
}

void set_samples(int s)
{
 number_of_samples=s;
}

void set_timesteps(int s)
{
 simulated_timesteps=s;
}

void set_mcmaze(string s)
{
strcpy(mc_mazes,s.c_str());
}

void set_seed(string s)
{
strcpy(seed_name,s.c_str());
if(strlen(seed_name)>0)
    seed_mode=true;
}

void set_fit_measure(string m)
{
if(m=="rnd")
   fitness_measure=fitness_rnd;
if(m=="std")
   fitness_measure=fitness_std;
if(m=="drift")
   fitness_measure=fitness_drift;
if(m=="goal")
   fitness_measure=fitness_goal;
cout << "Fitness measure " << fitness_measure << endl;
}

void set_nov_measure(string m)
{
if(m=="std" || m=="sample")
   novelty_measure=novelty_sample;
if(m=="accum")
   novelty_measure=novelty_accum;
if(m=="sample_free")
   novelty_measure=novelty_sample_free;
cout << "Novelty measure " << novelty_measure << endl;
}

static char output_dir[30]="";

static Environment* env;
static vector<Environment*> envList;
static vector<Environment*> mcList;

static int param=-1;
static int push_back_size = 20;

//used for discretization studies
double discretize(double x,long bins,double low, double high)
{
	double norm = x-low;
	double binsize = (high-low)/bins;
	int bin = (int)(norm/binsize);
	if(bin==bins)
		bin--;
        double result = (double)binsize*bin+binsize/2.0+low;
	return result;
}

long powerof2(int num)
{
long x=1;
if(num==0) return 1;
for(int y=0;y<num;y++)
	x*=2;
return x;
}

//novelty metric for maze simulation
float maze_novelty_metric(noveltyitem* x,noveltyitem* y)
{
	float diff = 0.0;
	for(int k=0;k<(int)x->data.size();k++) 
	{
		diff+=hist_diff(x->data[k],y->data[k]);
	}
	return diff;
}

static void read_in_environments(const char* mazefile, vector<Environment*>& envLst)
{
ifstream listfile(mazefile);

while(!listfile.eof())
{
string filename;
getline(listfile,filename);
if(filename.length() == 0)
 break;
cout << "Reading maze: " << filename << endl;
Environment* new_env = new Environment(filename.c_str());
envLst.push_back(new_env);
}

}

//novelty maze navigation run
Population *maze_novelty_realtime(char* outputdir,const char* mazefile,int par,const char* genes,bool novelty) {
	
    Population *pop;
    Genome *start_genome;
    char curword[20];

    int id;

	
	//crgate new maze environment
	//env=new Environment(mazefile);
        //read in environments
        read_in_environments(mazefile,envList);	

	if(strlen(mc_mazes)>0)
	 read_in_environments(mc_mazes,mcList);
        //param=par;
	push_back_size=par;
	if(outputdir!=NULL) strcpy(output_dir,outputdir);
		
    if(!seed_mode)
        strcpy(seed_name,genes);
	//starter genes file
    ifstream iFile(seed_name,ios::in);
	
    cout<<"START MAZE NAVIGATOR NOVELTY REAL-TIME EVOLUTION VALIDATION"<<endl;
if(!seed_mode)
 {
    cout<<"Reading in the start genome"<<endl;
}
else
    cout<<"Reading in the seed genome" <<endl;
    
  //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    iFile.close();

    cout<<"Start Genome: "<<start_genome<<endl;

    //Spawn the Population from starter gene
    cout<<"Spawning Population off Genome"<<endl;
    if(!seed_mode)
    pop=new Population(start_genome,NEAT::pop_size);
    else
    {
    pop=new Population(seed_name);//start_genome,NEAT::pop_size,0.0);   
    }
    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
      
    //Start the evolution loop using rtNEAT method calls 
    maze_novelty_realtime_loop(pop,novelty);

	//clean up
	delete env;
    return pop;
}

//actual rtNEAT loop for novelty maze navigation runs
int maze_novelty_realtime_loop(Population *pop,bool novelty) {
	bool firstflag=false; //indicates whether the maze has been solved yet
        bool weakfirst=false;	
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
  vector<Species*>::iterator curspec; //used in printing out debug info                                                         

  vector<Species*> sorted_species;  //Species sorted by max fit org in Species 

//was 1.0*number_of_samples+1.0 for earlier results...
   float archive_thresh=(1.0*number_of_samples+1.0);// * 20.0 * envList.size(); //initial novelty threshold
  //if(!minimal_criteria)
  //	archive_thresh*=20;
 //if(constraint_switch)
   //archive_thresh/=200.0;
  cout << "Archive threshold: " << archive_thresh << endl;
  //archive of novel behaviors
  noveltyarchive archive(archive_thresh,*maze_novelty_metric,true,push_back_size,minimal_criteria);
	
  data_rec Record; //stores run information
	
  int count=0;
  int pause;

  //Real-time evolution variables                                                                                             
  int offspring_count;
  Organism *new_org;

  //We try to keep the number of species constant at this number                                                    
  int num_species_target=NEAT::pop_size/15;
  
  //This is where we determine the frequency of compatibility threshold adjustment
  int compat_adjust_frequency = NEAT::pop_size/20;
  if (compat_adjust_frequency < 1)
    compat_adjust_frequency = 1;

//activity stat log file
  char asfn[100];
  sprintf(asfn,"%s_activitystats.dat",output_dir);
  ofstream activity_stat_file(asfn);
  if(activity_stats)
     reset_activity();  
  //Initially, we evaluate the whole population                                                                               
  //Evaluate each organism on a test                   
  int indiv_counter=0;  
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {

    //shouldn't happen                                                                                                        
    if (((*curorg)->gnome)==0) {
      cout<<"ERROR EMPTY GEMOME!"<<endl;
      cin>>pause;
    }

	//evaluate each individual
	(*curorg)->noveltypoint = maze_novelty_map((*curorg));
	(*curorg)->noveltypoint->indiv_number=indiv_counter;
	(*curorg)->fitness = (*curorg)->noveltypoint->fitness;
        //make sure that we are resetting novelty metric of
	//individuals that fail to meet goal?
	indiv_counter++;
  }

  if(novelty) {
  //assign fitness scores based on novelty
  archive.evaluate_population(pop,true);
  //add to archive
  archive.evaluate_population(pop,false);
  }

  if(novelty && minimal_criteria)
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) 
{
(*curorg)->fitness = SNUM/1000.0;
}  
//Get ready for real-time loop
  //Rank all the organisms from best to worst in each species
  pop->rank_within_species();                                                                            

  //Assign each species an average fitness 
  //This average must be kept up-to-date by rtNEAT in order to select species probabailistically for reproduction
  pop->estimate_all_averages();

  cout <<"Entering real time loop..." << endl;

  //Now create offspring one at a time, testing each offspring,                                                               
  // and replacing the worst with the new offspring if its better
  for 
(offspring_count=0;offspring_count<NEAT::pop_size*1001;offspring_count++) 
{
//fix compat_threshold, so no speciation...
//      NEAT::compat_threshold = 1000000.0;
	//only continue past generation 1000 if not yet solved
	//if(offspring_count>=pop_size*1000 && firstflag)
        // if(firstflag)
	// break;

int evolveupdate=50000;
if(offspring_count % evolveupdate ==0) {
   char fn[100];
   sprintf(fn,"%s_evolvability%d.dat",output_dir,offspring_count/evolveupdate);
   for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
    evolvability(*curorg,fn);
   }
}

if(activity_stats&& offspring_count % 10000 == 0)
{
  pop->update_statistics();
 activity_stat_file << offspring_count << " "  << calculate_diversity() 
                   << " " << calculate_cumulative_activity() << " " <<
	                     calculate_average_activity() << endl; 
}	
	
	//end of generation
    if(offspring_count % (NEAT::pop_size*1) == 0)
	{
                 if(novelty) {
			archive.end_of_gen_steady(pop);
			//archive.add_randomly(pop);
			archive.evaluate_population(pop,false);
			cout << "ARCHIVE SIZE:" << 
			archive.get_set_size() << endl;
                 }
                 cout << "GEN" << offspring_count/NEAT::pop_size << endl;
	}

	//write out current generation and fittest individuals
    if( offspring_count % (NEAT::pop_size*NEAT::print_every) == 0 )
	{
        	cout << offspring_count << endl;
			char fname[100];
			sprintf(fname,"%sarchive.dat",output_dir);
			archive.Serialize(fname);		
	
			sprintf(fname,"%sfittest_%d",output_dir,offspring_count/NEAT::pop_size);
			archive.serialize_fittest(fname);

			sprintf(fname,"%sgen_%d",output_dir,offspring_count/NEAT::pop_size);
			pop->print_to_file_by_species(fname);
	

                        sprintf(fname,"%srecord.dat",output_dir);
                        Record.serialize(fname);
      }
	
    //Every pop_size reproductions, adjust the compat_thresh to better match the num_species_targer
    //and reassign the population to new species                                              
    if (offspring_count % compat_adjust_frequency == 0) {
		count++;
		#ifdef DEBUG_OUTPUT
		cout << "Adjusting..." << endl;
		#endif
          if(novelty) {	
	   //update fittest individual list		
	   archive.update_fittest(pop);
	   //refresh generation's novelty scores
	   archive.evaluate_population(pop,true);
          }	
      int num_species = pop->species.size();
      double compat_mod=0.1;  //Modify compat thresh to control speciation                                                     
      // This tinkers with the compatibility threshold 
      if (num_species < num_species_target) {
	NEAT::compat_threshold -= compat_mod;
      }
      else if (num_species > num_species_target)
	NEAT::compat_threshold += compat_mod;

      if (NEAT::compat_threshold < 0.3)
	NEAT::compat_threshold = 0.3;
		#ifdef DEBUG_OUTPUT
      cout<<"compat_thresh = "<<NEAT::compat_threshold<<endl;
		#endif
	  
      //Go through entire population, reassigning organisms to new species                                                  
      for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
			pop->reassign_species(*curorg);
      }
    }
    

    //For printing only
	#ifdef DEBUG_OUTPUT
    for(curspec=(pop->species).begin();curspec!=(pop->species).end();curspec++) {
      cout<<"Species "<<(*curspec)->id<<" size"<<(*curspec)->organisms.size()<<" average= "<<(*curspec)->average_est<<endl;
    }

    cout<<"Pop size: "<<pop->organisms.size()<<endl;
	#endif
	
    //Here we call two rtNEAT calls: 
    //1) choose_parent_species() decides which species should produce the next offspring
    //2) reproduce_one(...) creates a single offspring fromt the chosen species
    new_org=(pop->choose_parent_species())->reproduce_one(offspring_count,pop,pop->species);

    //Now we evaluate the new individual
    //Note that in a true real-time simulation, evaluation would be happening to all individuals at all times.
    //That is, this call would not appear here in a true online simulation.
	#ifdef DEBUG_OUTPUT
    cout<<"Evaluating new baby: "<<endl;
	#endif
	
	data_record* newrec=new data_record();
	newrec->indiv_number=indiv_counter;
	//evaluate individual, get novelty point
	new_org->noveltypoint = maze_novelty_map(new_org,newrec);
	new_org->noveltypoint->indiv_number = indiv_counter;
	new_org->fitness=new_org->noveltypoint->fitness;
        //calculate novelty of new individual
	if(novelty) {
        archive.evaluate_individual(new_org,pop->organisms);
	//newrec->ToRec[5] = archive.get_threshold();
	newrec->ToRec[6] = archive.get_set_size();
	newrec->ToRec[RECSIZE-2] = new_org->noveltypoint->novelty;
        }
	if( !new_org->noveltypoint->viable && minimal_criteria)
	{
		new_org->fitness = SNUM/1000.0;
                //new_org->novelty = 0.00000001;
                //reset behavioral characterization
                new_org->noveltypoint->reset_behavior();
                //cout << "fail" << endl;
               // cout << " :( " << endl;
	}	
        else
        {
             //cout << ":)" << new_org->noveltypoint->indiv_number << endl;
        }	
        //add record of new indivdual to storage
	Record.add_new(newrec);
	indiv_counter++;
	
	//update fittest list
	archive.update_fittest(new_org);
	#ifdef DEBUG_OUTPUT
	cout << "Fitness: " << new_org->fitness << endl;
	cout << "Novelty: " << new_org->noveltypoint->novelty << endl;
	cout << "RFit: " << new_org->noveltypoint->fitness << endl;
    #endif
	
    //Now we reestimate the baby's species' fitness
    new_org->species->estimate_average();
        if(!weakfirst && (newrec->ToRec[3]>=envList.size())) {
		weakfirst=true;
		char filename[100];
		sprintf(filename,"%srtgen_weakfirst",output_dir);
		pop->print_to_file_by_species(filename);
		cout << "Maze weakly solved by indiv# " << indiv_counter << endl;	
//disable quit for now
          if(fitness_measure == fitness_goal && false) 
           firstflag=true; 
        }
	//write out the first individual to solve maze
	if(!firstflag && (newrec->ToRec[3]>=envList.size() && newrec->ToRec[4]>=envList.size())) {
		firstflag=true;
		char filename[100];
		sprintf(filename,"%srtgen_first",output_dir);
		pop->print_to_file_by_species(filename);
		cout << "Maze solved by indiv# " << indiv_counter << endl;	
		//break;
	}

    //Remove the worst organism                                                                                               
    if(rand_repl || fitness_measure ==fitness_rnd)
     pop->remove_random();
    else    
     pop->remove_worst();

  }
  
  //write out run information, archive, and final generation
  cout << "COMPLETED...";
  char filename[100];
  sprintf(filename,"%srecord.dat",output_dir);
  char fname[100];
  sprintf(fname,"%srtarchive.dat",output_dir);
  archive.Serialize(fname);
  Record.serialize(filename);
  
  sprintf(fname,"%sfittest_final",output_dir);
  archive.serialize_fittest(fname);

  sprintf(fname,"%srtgen_final",output_dir);
  pop->print_to_file_by_species(fname);
  exit(0);
  return 0;
}
  
//initialize the maze simulation
Environment* mazesimIni(Environment* tocopy,Network *net, vector< vector<float> > &dc)
  {
    double inputs[20];
    Environment *newenv= new Environment(*tocopy);
	 
	//flush the neural net
	net->flush();
	//update the maze
	newenv->Update();
	//create neural net inputs
	newenv->generate_neural_inputs(inputs);
	//load into neural net
	net->load_sensors(inputs);
	
	//propogate input through net
    for(int i=0;i<10;i++)
		net->activate();
	
	return newenv;
  }
  
  //execute a timestep of the maze simulation evaluation
  double mazesimStep(Environment* newenv,Network *net,vector< vector<float> > &dc)
  {
	  double inputs[20];
	  
	  
		newenv->generate_neural_inputs(inputs);
		net->load_sensors(inputs);
		net->activate();
	
	  	//use the net's outputs to change heading and velocity of navigator
		newenv->interpret_outputs(net->outputs[0]->activation,net->outputs[1]->activation,0.0); //net->outputs[2]->activation);
	  	//update the environment
		newenv->Update();
	  newenv->distance_to_poi(); 
	  double dist = newenv->distance_to_target();
	  if(dist<=1) dist=1;
	  double fitness = 5.0/dist; //used for accumulated fitness (obselete)
	  
	  return fitness;
  }
double mazesim(Network* net, vector< vector<float> > &dc, data_record *record,Environment* the_env)
{
	vector<float> data;
	
	int timesteps=the_env->steps; //simulated_timesteps;
	int stepsize=10000;
	
	double fitness=0.0;
	Environment *newenv;
        position_accumulator *accum; 
	
	newenv=mazesimIni(the_env,net,dc);
        newenv->goalattract = goal_attract;	
	//data collection vector initialization
	//dc.clear();

	if(novelty_measure == novelty_sample || 
           novelty_measure ==novelty_sample_free)
		data.reserve(timesteps/stepsize);
	if(novelty_measure == novelty_accum)
	{
		data.reserve(100);	
		float minx,miny,maxx,maxy;
		newenv->get_range(minx,miny,maxx,maxy);
		vector<int> dims;
		dims.push_back(10);
		dims.push_back(10);
		accum=new position_accumulator(dims,minx,miny,maxx,maxy);
	}

	/*ENABLE FOR ADDT'L INFO STUDIES*/
        if(number_of_samples>0)	
	 stepsize=timesteps/number_of_samples;
		
	for(int i=0;i<timesteps;i++)
	{
		fitness+=mazesimStep(newenv,net,dc);
		//if taking additional samples, collect during run
		if(novelty_measure==novelty_sample ||
                   novelty_measure==novelty_sample_free)
		if((timesteps-i-1)%stepsize==0)
		{
				data.push_back(newenv->hero.location.x);
				data.push_back(newenv->hero.location.y);
		}
		
		float loc[2]={newenv->hero.location.x,newenv->hero.location.y};
		if(novelty_measure==novelty_accum)
		{
			accum->add_point(loc);
		}
	}

	//calculate fitness of individual as closeness to target
	if(fitness_measure == fitness_goal)
        {
	fitness=500.0 - newenv->distance_to_target();
	if(fitness<0.1) fitness=0.1;
        //if (newenv->hero.collide)
	//	fitness+=50;
	}
        
        if(fitness_measure ==fitness_rnd)
        {
         //todo assign random fitness, this needs to get reassigned
         //often...
   	 fitness = randint(10,100);
        }
   
	//calculate fitness as meeting minimal criteria
	if(fitness_measure == fitness_drift)
	{
	if(newenv->reachgoal)
	{
		fitness=1.0;
                if(newenv->reachpoi)
			fitness=500.0;
	}
	else
	{
		fitness=SNUM/1000.0;
	}
	}

        if(fitness_measure == fitness_std)
        {
          fitness=SNUM;
          float mod = 500.0 - newenv->closest_to_target;
          if(mod<0) mod=0.0;
          float mod2 = 500.0 - newenv->closest_to_poi; 
          if(mod2<0) mod2=0.0;

	  fitness+=mod;
	  fitness+=mod2;
/*
          if(newenv->reachgoal)
		fitness+=250.0;
          else
                fitness+=mod;

           if(newenv->reachpoi)
		fitness+=250.0;
           else if (newenv->reachgoal)
                fitness+=mod2;
*/
        }

	//fitness as novelty studies
	//data.push_back(fitness);
	
	float x=newenv->hero.location.x;
	float y=newenv->hero.location.y;
	
	/* ENABLE FOR DISCRETIZATION STUDIES
	if(param>0)
	{
	 long bins=powerof2(param);
	 x=discretize(x,bins,0.0,200.0);
	 y=discretize(y,bins,0.0,200.0);
	}
	*/
	
	if(novelty_measure==novelty_sample || novelty_measure==novelty_sample_free)	
	if(false)
	{
	//novelty point is the ending location of the navigator
	data.push_back(x);
	data.push_back(y);
	}

	if(novelty_measure==novelty_accum)
	{
		accum->transform();	
		for(int x=0;x<accum->size;x++)
			data.push_back(accum->buffer[x]);
	}

	if(record!=NULL)
	{
		record->ToRec[0]+=fitness;
		record->ToRec[1]=newenv->hero.location.x;
		record->ToRec[2]=newenv->hero.location.y;
		record->ToRec[3]+=newenv->reachgoal;
		record->ToRec[4]+=newenv->reachpoi;
                record->ToRec[5]=newenv->hero.collide;		
	}

	if(novelty_measure==novelty_accum)
		delete accum;

	dc.push_back(data);
	
        delete newenv;
	return fitness;
}

bool contained(float x,float y)
{
Point p(x,y);
if (x>200)
 return false;
if (x<0)
 return false;
if (y>200)
 return false;
if (y<0)
 return false;

return true;
/*
if(envList[0]->end.distance(p) < 400)
{
 //cout <<"contained.." << endl;
 return true;
}
//cout <<"notcontained..." << endl;
return false;
*/
}
#define DISABLE_PROB 0.05
void mutate_genome(Genome* new_genome)
{
	double mut_power=NEAT::weight_mut_power;
			//NOTE:  A link CANNOT be added directly after a node was added because the phenotype
			//       will not be appropriately altered to reflect the change
			if(1) {
				//If we didn't do a structural mutation, we do the other kinds

				if (randfloat()<NEAT::mutate_random_trait_prob) {
					//cout<<"mutate random trait"<<endl;
					new_genome->mutate_random_trait();
				}
				if (randfloat()<NEAT::mutate_link_trait_prob) {
					//cout<<"mutate_link_trait"<<endl;
					new_genome->mutate_link_trait(1);
				}
				if (randfloat()<NEAT::mutate_node_trait_prob) {
					//cout<<"mutate_node_trait"<<endl;
					new_genome->mutate_node_trait(1);
					new_genome->mutate_node_parameters(NEAT::time_const_mut_power,NEAT::time_const_mut_prob,
					                                   NEAT::bias_mut_power,NEAT::bias_mut_prob);
				}
				if (randfloat()<NEAT::mutate_link_weights_prob) {
					//cout<<"mutate_link_weights"<<endl;
					new_genome->mutate_link_weights(mut_power,1.0,GAUSSIAN);
				}
				if (randfloat()<NEAT::mutate_toggle_enable_prob) {
					//cout<<"mutate toggle enable"<<endl;
					new_genome->mutate_toggle_enable(1);
					new_genome->struct_change=3; //JLADD

				}
                                if (randfloat()<DISABLE_PROB)
                                   new_genome->mutate_gene_disable();
				if (randfloat()<NEAT::mutate_gene_reenable_prob) {
					//cout<<"mutate gene reenable"<<endl;
					new_genome->mutate_gene_reenable();
					new_genome->struct_change=3; //JLADD
				}
			}
}
 
double evolvability(Organism* org,char* fn) {
 fstream file; 
 file.open(fn,ios::app|ios::out);
 cout <<"Evolvability..." << endl;
 file << "---" << endl;
  Organism *new_org= new Organism(*org);
 for(int i=0;i<200;i++) { 
  new_org->gnome = new Genome(*org->gnome);
  if(i!=0) //first copy is clean
  for(int j=0;j<3;j++) mutate_genome(new_org->gnome);
  new_org->net=new_org->gnome->genesis(0);
  noveltyitem* nov_item = maze_novelty_map(new_org);
  for(int k=0;k<nov_item->data[0].size();k++)
    file << nov_item->data[0][k] << " ";
  delete new_org->gnome;
  delete new_org->net;
  delete nov_item;
  file << endl;
 }  
 file.close();
 return 0.0;
}


//evaluates an individual and stores the novelty point
noveltyitem* maze_novelty_map(Organism *org,data_record* record)
{
  static int best = 0;
  noveltyitem *new_item = new noveltyitem;
  new_item->genotype=new Genome(*org->gnome);
  new_item->phenotype=new Network(*org->net);
  vector< vector<float> > gather;

  vector<float> constraint_vector;
  bool apply_constraints=constraint_switch; //false;
  bool remove_regular=constraint_switch; //false;
  int c1old=0,c2old=0;
  double fitness=0.0;
  static float highest_fitness=0.0;
  
new_item->viable=true;

  if(record!=NULL && minimal_criteria) 
  for(int x=0;x<mcList.size();x++)
  {
  record->ToRec[3]=0;
   mazesim(org->net,gather,record,mcList[x]); 
  if(!record->ToRec[3]) {
    new_item->viable=false;
   // cout << "not viable..." << endl;
    break;   
  }
  //cout << "viable..." << endl;
  }

 gather.clear();
 if(record!=NULL) {
 record->ToRec[0]=0;
 record->ToRec[3]=0;
 record->ToRec[4]=0;
 }

if(true) //new_item->viable)
	for(int x=0;x<envList.size();x++)
        {
         if(record!=NULL) {
           c1old = record->ToRec[3];
           c2old = record->ToRec[4];
 	 }

         fitness+=mazesim(org->net,gather,record,envList[x]);

         if(record!=NULL) {
	   constraint_vector.push_back(record->ToRec[3]-c1old);
           constraint_vector.push_back(record->ToRec[4]-c2old);
           c1old=record->ToRec[3];
           //if(record->ToRec[5]==1)
           //  new_item->viable=false;
         }
	else { constraint_vector.push_back(0);
	 constraint_vector.push_back(0); }
        } 
else 
{
for(int x=0;x<envList.size();x++) constraint_vector.push_back(0); 
}

         //minimal criteria must be met in *all* scenarios...
         
        if(record!=NULL)
         {
           if(area_of_interest)
           {
            if(!contained(record->ToRec[1],record->ToRec[2]))
 	    {
         	new_item->viable=false;
         	//cout << "not viable..." << endl;
            }
            //else cout << "viable..." << endl;
           }
            else
           if( record->ToRec[3]<envList.size())
           { new_item->viable=false;
            //cout << record->ToRec[3] << endl; 
           }
            else { 
             //cout << "viable... " << endl;
            }
          }
  	
        if(fitness>highest_fitness)
		highest_fitness=fitness;
	
	//keep track of highest fitness so far in record
	if(record!=NULL)
	{
        record->ToRec[5]=new_item->viable;
	if(record->ToRec[3]>best)
        {
         best=record->ToRec[3];
         cout << "best: " << best << endl;
 	}
	record->ToRec[RECSIZE-1]=highest_fitness;
	}

 	 //push back novelty characterization
         if(!remove_regular)
         for(int i=0;i<gather.size();i++)
	  new_item->data.push_back(gather[i]);
         
         if(apply_constraints)
	  new_item->data.push_back(constraint_vector);
 	 //set fitness (this is 'real' objective-based fitness, not novelty)
 	 new_item->fitness=fitness;
         org->fitness=fitness;
 
  return new_item;
}


//Perform evolution on single pole balacing, for gens generations
Population *maze_generational(char* outputdir,const char* mazefile,int param,const char *genes, int gens,bool novelty) 
{
    float archive_thresh=3.0;

    noveltyarchive archive(archive_thresh,*maze_novelty_metric,true,push_back_size,minimal_criteria,true);


    Population *pop;

    Genome *start_genome;
    char curword[20];
    int id;
    data_rec Record;

    ostringstream *fnamebuf;
    int gen;

    ifstream iFile(genes,ios::in);

    read_in_environments(mazefile,envList);	

    if(strlen(mc_mazes)>0)
	 read_in_environments(mc_mazes,mcList);
    
    if(outputdir!=NULL) strcpy(output_dir,outputdir);
    cout<<"START GENERATIONAL MAZE EVOLUTION"<<endl;

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    iFile.close();
  
      cout<<"Start Genome: "<<start_genome<<endl;
      
      //Spawn the Population
      cout<<"Spawning Population off Genome"<<endl;
      
      pop=new Population(start_genome,NEAT::pop_size);
      
      cout<<"Verifying Spawned Pop"<<endl;
      pop->verify();

      for (gen=0;gen<=gens;gen++) {
	cout<<"Generation "<<gen<<endl;
	bool win = maze_generational_epoch(pop,gen,Record,archive,novelty);

  //writing out stuff 
  if(gen%NEAT::print_every == 0 )
  {
  char filename[100];
  sprintf(filename,"%s_record.dat",output_dir);
  char fname[100];
  sprintf(fname,"%s_archive.dat",output_dir);
  archive.Serialize(fname);
  Record.serialize(filename);
  sprintf(fname,"%sgen_%d",output_dir,gen);
  pop->print_to_file_by_species(fname);
  }

  if (win)
  {
   char fname[100];
   sprintf(fname,"%s_wingen",output_dir);
    ofstream winfile(fname);
    winfile << gen << endl;
  sprintf(fname,"%s_archive.dat",output_dir);
  archive.Serialize(fname);
  sprintf(fname,"%s_record.dat",output_dir);
  Record.serialize(fname);
    //break;
  }
      
}


    return pop;

}

int maze_generational_epoch(Population *pop,int generation,data_rec& Record, noveltyarchive& archive, bool novelty) {
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
static double best_fitness =0.0;  
static vector<Organism*> measure_pop;
  //char cfilename[100];
  //strncpy( cfilename, filename.c_str(), 100 );

  //ofstream cfilename(filename.c_str());

  static bool win=false;
  static bool firstflag=false;
  int winnernum;
  int indiv_counter=0;
 
  int evolveupdate=100;
  if(generation%evolveupdate==0)
  {
   char fn[100];
   sprintf(fn,"%s_evolvability%d.dat",output_dir,generation/evolveupdate);
   for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
    evolvability(*curorg,fn);
     }
   }
  
  //Evaluate each organism on a test
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
        
        data_record* newrec=new data_record();
	newrec->indiv_number=indiv_counter;
	//evaluate individual, get novelty point    
	(*curorg)->noveltypoint = maze_novelty_map((*curorg),newrec);
        (*curorg)->noveltypoint->indiv_number = indiv_counter;
        (*curorg)->datarec = newrec;	

	if((newrec->ToRec[3]>=envList.size())) { // && newrec->ToRec[4]>0.0)) {
	        if((newrec->ToRec[4]>=envList.size())) {	
                (*curorg)->winner=true;
                win=true;	
                }
                if(fitness_measure==fitness_goal) {
                (*curorg)->winner=true;
                win=true;
                }
	}
  
        if((*curorg)->noveltypoint->fitness > best_fitness)
	{
		best_fitness = (*curorg)->noveltypoint->fitness;
		cout << "NEW BEST " << best_fitness << endl;
	}	
        //add record of new indivdual to storage
	Record.add_new(newrec);
	indiv_counter++;
/*
	if((*curorg)->noveltypoint->viable)
             cout << "viable..." << endl;
        else
             cout << "not viable..." << endl;
*/
	if( !(*curorg)->noveltypoint->viable && minimal_criteria)
	{
		(*curorg)->fitness = SNUM/1000.0;
                //new_org->novelty = 0.00000001;
                //reset behavioral characterization
                (*curorg)->noveltypoint->reset_behavior();
                //cout << "fail" << endl;
               // cout << " :( " << endl;
	}	
       	
	//update fittest list
	archive.update_fittest(*curorg);
	
        if(!novelty)
    	   (*curorg)->fitness = (*curorg)->noveltypoint->fitness;
  }
 
  if(novelty)
  {
  
	//NEED TO CHANGE THESE TO GENERATIONAL EQUIVALENTS...
	//assign fitness scores based on novelty
 	archive.evaluate_population(pop,true);
	///now add to the archive (maybe remove & only add randomly?)
	archive.evaluate_population(pop,false);
       
 
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
	if( !(*curorg)->noveltypoint->viable && minimal_criteria)
	{
		(*curorg)->fitness = SNUM/1000.0;
                //new_org->novelty = 0.00000001;
                //reset behavioral characterization
                //cout << "fail" << endl;
               // cout << " :( " << endl;
	}	
      } 
        cout << "ARCHIVE SIZE:" << archive.get_set_size() << endl;  
        cout << "THRESHOLD:" << archive.get_threshold() << endl;
	archive.end_of_gen_steady(pop);
	//adjust novelty of infeasible individuals
	/*
	if(!newrec->ToRec[3] && novelty_measure != novelty_sample_free)
	{
		(*curorg)->fitness = 0.00001;
	}
	*/	
  }

  char fn[100];
  sprintf(fn,"%sdist%d",output_dir,generation);
  //pop->print_distribution(fn);
  
  //Average and max their fitnesses for dumping to file and snapshot
  for(curspecies=(pop->species).begin();curspecies!=(pop->species).end();++curspecies) {

    //This experiment control routine issues commands to collect ave
    //and max fitness, as opposed to having the snapshot do it, 
    //because this allows flexibility in terms of what time
    //to observe fitnesses at

    (*curspecies)->compute_average_fitness();
    (*curspecies)->compute_max_fitness();
  }

  //Take a snapshot of the population, so that it can be
  //visualized later on
  //if ((generation%1)==0)
  //  pop->snapshot();

  //Only print to file every print_every generations
 // if  (win||
 //      ((generation%(NEAT::print_every))==0))
if(win && !firstflag)
{
    for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
      if ((*curorg)->winner) {
	winnernum=((*curorg)->gnome)->genome_id;
	cout<<"WINNER IS #"<<((*curorg)->gnome)->genome_id<<endl;
	char filename[100];
        sprintf(filename,"%s_winner", output_dir);	
       (*curorg)->print_to_file(filename);
      }
    }    
   firstflag = true;
}


 
 
 //Create the next generation
  pop->epoch(generation);


  return win;
}

