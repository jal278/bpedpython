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

using namespace std;

enum novelty_measure_type { novelty_sample, novelty_accum, novelty_sample_free };
static novelty_measure_type novelty_measure = novelty_sample;

enum fitness_measure_type { fitness_goal, fitness_drift, fitness_std };
static fitness_measure_type fitness_measure = fitness_std;

static int number_of_samples = 1;
static int simulated_timesteps = 400;

void set_samples(int s)
{
 number_of_samples=s;
}

void set_timesteps(int s)
{
 simulated_timesteps=s;
}

void set_fit_measure(string m)
{
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
static int param=-1;
static int push_back_size = 200;
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

//fitness simulation of maze navigation
Population *maze_fitness_realtime(char* outputdir,const char *mazefile,int par,const char* genes) {
    Population *pop;
    Genome *start_genome;
    char curword[20];
    int id;

	//create new maze environment
	env=new Environment(mazefile);
	if(outputdir!=NULL) strcpy(output_dir,outputdir);
	param=par;
	
	//starter gene file
    ifstream iFile(genes,ios::in);
	
    cout<<"START MAZE NAVIGATOR FITNESS REAL-TIME EVOLUTION VALIDATION"<<endl;

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    iFile.close();

    cout<<"Start Genome: "<<start_genome<<endl;

    //Spawn the Population from starter gene
    cout<<"Spawning Population off Genome"<<endl;
    pop=new Population(start_genome,NEAT::pop_size);
      
    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
      
    //Start the evolution loop using rtNEAT method calls 
    maze_fitness_realtime_loop(pop);

	//clean up
	delete env;
    return pop;
}

//actual rtNEAT loop for fitness run of maze navigation
int maze_fitness_realtime_loop(Population *pop) {
  bool firstflag=false; //indicates whether maze has been solved yet
  int indiv_counter=0;
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
  vector<Species*>::iterator curspec; //used in printing out debug info                                                         
  vector<Species*> sorted_species;  //Species sorted by max fit org in Species                                                  

  data_rec Record; //holds run information
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

  //Initially, we evaluate the whole population                                                                               
  //Evaluate each organism on a test                                                                                          
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {

    //shouldn't happen                                                                                                        
    if (((*curorg)->gnome)==0) {
      cout<<"ERROR EMPTY GEMOME!"<<endl;
      cin>>pause;
    }

	//map the novelty point to each individual (this runs the maze simulation)
	(*curorg)->noveltypoint = maze_novelty_map((*curorg));
	(*curorg)->noveltypoint->indiv_number = indiv_counter;
	indiv_counter++;
	(*curorg)->fitness = (*curorg)->noveltypoint->fitness;
    
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

  //run for 2000 generations (250*2000 = 500,000 evaluations)
  for 
(offspring_count=0;offspring_count<NEAT::pop_size*2001;offspring_count++) 
{
    
	if(offspring_count>=pop_size*1000 && firstflag)
		break;
    if(offspring_count % (NEAT::pop_size*NEAT::print_every) == 0 )
	{
			cout << offspring_count << endl;
			char fname[30];
			sprintf(fname,"%sgen_%d",output_dir,offspring_count/NEAT::pop_size);
			pop->print_to_file_by_species(fname);


    char filename[30];
    sprintf(filename,"%srecord.dat",output_dir);

    Record.serialize(filename);
	}
	
    //Every pop_size reproductions, adjust the compat_thresh to better match the num_species_targer
    //and reassign the population to new species                                              
    if (offspring_count % compat_adjust_frequency == 0) {
		count++;
		#ifdef DEBUG_OUTPUT
		cout << "Adjusting..." << endl;
		#endif
		
			
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

    #ifdef DEBUG_OUTPUT
	cout<<"Evaluating new baby: "<<endl;
	#endif
	
	//create record for new individual
	data_record *newrec = new data_record();
	newrec->indiv_number=indiv_counter;
	//evaluate new individual
	new_org->noveltypoint = maze_novelty_map(new_org,newrec);
	new_org->noveltypoint->indiv_number = indiv_counter;
	//grab novelty
	newrec->ToRec[RECSIZE-2]=new_org->noveltypoint->novelty;
	//set organism's fitness 
	new_org->fitness = new_org->noveltypoint->fitness;
	indiv_counter++;
	
	#ifdef DEBUG_OUTPUT
	cout << "Fitness: " << new_org->fitness << endl;
	cout << "Novelty: " << new_org->noveltypoint->novelty << endl;
	cout << "RFit: " << new_org->noveltypoint->fitness << endl;
	#endif
	
	//add record of new individual to storage
	Record.add_new(newrec);
	
    //Now we reestimate the baby's species' fitness
    new_org->species->estimate_average();

    //Remove the worst organism                                                                                               
    pop->remove_worst();
	
	//store first solution organism
	if(!firstflag && (newrec->ToRec[3]>0.0 && newrec->ToRec[4]>0.0)) {
		firstflag=true;
		char filename[30];
		sprintf(filename,"%sgen_first",output_dir);
		pop->print_to_file_by_species(filename);
		cout << "Maze solved by indiv# " << indiv_counter << endl;	
	}
  	
  }
  
  //finish up, write out the record and the final generation
	cout << "COMPLETED..." << endl;
    char filename[30];
    sprintf(filename,"%srecord.dat",output_dir);
    Record.serialize(filename);
	sprintf(filename,"%sfit_rtgen_final",output_dir);
	pop->print_to_file_by_species(filename);
  
   return 0;
}

//novelty maze navigation run
Population *maze_novelty_realtime(char* outputdir,const char* mazefile,int par,const char* genes) {
	
    Population *pop;
    Genome *start_genome;
    char curword[20];
    int id;

	//crgate new maze environment
	env=new Environment(mazefile);
	//param=par;
	push_back_size=par;
	if(outputdir!=NULL) strcpy(output_dir,outputdir);
		
	//starter genes file
    ifstream iFile(genes,ios::in);
	
    cout<<"START MAZE NAVIGATOR NOVELTY REAL-TIME EVOLUTION VALIDATION"<<endl;

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    iFile.close();

    cout<<"Start Genome: "<<start_genome<<endl;

    //Spawn the Population from starter gene
    cout<<"Spawning Population off Genome"<<endl;
    pop=new Population(start_genome,NEAT::pop_size);
      
    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
      
    //Start the evolution loop using rtNEAT method calls 
    maze_novelty_realtime_loop(pop);

	//clean up
	delete env;
    return pop;
}

//actual rtNEAT loop for novelty maze navigation runs
int maze_novelty_realtime_loop(Population *pop) {
	bool firstflag=false; //indicates whether the maze has been solved yet
	
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
  vector<Species*>::iterator curspec; //used in printing out debug info                                                         

  vector<Species*> sorted_species;  //Species sorted by max fit org in Species 

   float archive_thresh=1.0; //initial novelty threshold

  //archive of novel behaviors
  noveltyarchive archive(archive_thresh,*maze_novelty_metric,true,push_back_size);
	
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
	//make sure that we are resetting novelty metric of
	//individuals that fail to meet goal?
	indiv_counter++;
  }

  //assign fitness scores based on novelty
  archive.evaluate_population(pop,true);
  //add to archive
  archive.evaluate_population(pop,false);
 
  
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) (*curorg)->fitness = 0.00000001;
  
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
(offspring_count=0;offspring_count<NEAT::pop_size*2001;offspring_count++) 
{
	//only continue past generation 1000 if not yet solved
	if(offspring_count>=pop_size*1000 && firstflag)
		break;
	
	
	//end of generation
    if(offspring_count % (NEAT::pop_size*1) == 0)
	{
			archive.end_of_gen_steady(pop);
			//archive.add_randomly(pop);
			archive.evaluate_population(pop,false);
			cout << "ARCHIVE SIZE:" << 
			archive.get_set_size() << endl;
	}

	//write out current generation and fittest individuals
    if( offspring_count % (NEAT::pop_size*NEAT::print_every) == 0 )
	{
        	cout << offspring_count << endl;
			char fname[30];
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
	
	   //update fittest individual list		
	   archive.update_fittest(pop);
	   //refresh generation's novelty scores
	   archive.evaluate_population(pop,true);
			
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
	//calculate novelty of new individual
	archive.evaluate_individual(new_org,pop->organisms);
	newrec->ToRec[5] = archive.get_threshold();
	newrec->ToRec[6] = archive.get_set_size();
	newrec->ToRec[RECSIZE-2] = new_org->noveltypoint->novelty;

	if(!newrec->ToRec[3] && novelty_measure != novelty_sample_free)
	{
		new_org->fitness = 0.00001;
                //cout << "fail" << endl;
	}	
	else
	{
		//cout << "succeed" << endl;
		//cout << new_org->fitness << endl;
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

	//write out the first individual to solve maze
	if(!firstflag && (newrec->ToRec[3]>0.0 && newrec->ToRec[4]>0.0)) {
		firstflag=true;
		char filename[30];
		sprintf(filename,"%srtgen_first",output_dir);
		pop->print_to_file_by_species(filename);
		cout << "Maze solved by indiv# " << indiv_counter << endl;	
		//break;
	}

    //Remove the worst organism                                                                                               
    pop->remove_worst();
	

  }
  
  //write out run information, archive, and final generation
  cout << "COMPLETED...";
  char filename[30];
  sprintf(filename,"%srecord.dat",output_dir);
  char fname[30];
  sprintf(fname,"%srtarchive.dat",output_dir);
  archive.Serialize(fname);
  Record.serialize(filename);
  
  sprintf(fname,"%sfittest_final",output_dir);
  archive.serialize_fittest(fname);

  sprintf(fname,"%srtgen_final",output_dir);
  pop->print_to_file_by_species(fname);

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
		newenv->interpret_outputs(net->outputs[0]->activation,net->outputs[1]->activation);
	  	//update the environment
		newenv->Update();
	  newenv->distance_to_poi(); 
	  double dist = newenv->distance_to_target();
	  if(dist<=1) dist=1;
	  double fitness = 5.0/dist; //used for accumulated fitness (obselete)
	  
	  return fitness;
  }
double mazesim(Network* net, vector< vector<float> > &dc, data_record *record)
{
	vector<float> data;
	
	int timesteps=simulated_timesteps;
	int stepsize=10000;
	
	double fitness=0.0;
	Environment *newenv;
        position_accumulator *accum; 
	
	newenv=mazesimIni(env,net,dc);
	
	//data collection vector initialization
	dc.clear();
	dc.push_back(data);

	if(novelty_measure == novelty_sample || 
           novelty_measure ==novelty_sample_free)
		dc[0].reserve(timesteps/stepsize);
	if(novelty_measure == novelty_accum)
	{
		dc[0].reserve(100);	
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
				dc[0].push_back(newenv->hero.location.x);
				dc[0].push_back(newenv->hero.location.y);
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
	fitness=300.0 - newenv->distance_to_target();
	if(fitness<0.1) fitness=0.1;
	}

	//calculate fitness as meeting minimal criteria
	if(fitness_measure == fitness_drift)
	{
	if(newenv->reachgoal)
	{
		fitness=1.0;
	}
	else
	{
		fitness=0.00001;
	}
	}

        if(fitness_measure == fitness_std)
        {
          fitness=0.0;
          if(newenv->reachgoal)
		fitness+=250.0;
          else
                fitness+=250.0 - newenv->distance_to_target();

           if(newenv->reachpoi)
		fitness+=250.0;
           else
                fitness+=250.0 - newenv->closest_to_poi;
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
	dc[0].push_back(x);
	dc[0].push_back(y);
	}

	if(novelty_measure==novelty_accum)
	{
		accum->transform();	
		for(int x=0;x<accum->size;x++)
			dc[0].push_back(accum->buffer[x]);
	}

	if(record!=NULL)
	{
		record->ToRec[0]=fitness;
		record->ToRec[1]=newenv->hero.location.x;
		record->ToRec[2]=newenv->hero.location.y;
		record->ToRec[3]=newenv->reachgoal;
		record->ToRec[4]=newenv->reachpoi;		
	}

	if(novelty_measure==novelty_accum)
		delete accum;
	delete newenv;
	return fitness;
}

//evaluates an individual and stores the novelty point
noveltyitem* maze_novelty_map(Organism *org,data_record* record)
{

  noveltyitem *new_item = new noveltyitem;
  new_item->genotype=new Genome(*org->gnome);
  new_item->phenotype=new Network(*org->net);
  vector< vector<float> > gather;

  double fitness;
  static float highest_fitness=0.0;

	fitness=mazesim(org->net,gather,record);
  	if(fitness>highest_fitness)
		highest_fitness=fitness;
	
	//keep track of highest fitness so hard in record
	if(record!=NULL)
	{

	//TODO: Perhaps remove non-viable organisms from merged populations
	//or set their behavioral characterization to some null value..
	if(!record->ToRec[3])
		new_item->viable=false;

		/*
	record->ToRec[19]=org->gnome->parent1;
	record->ToRec[18]=org->gnome->parent2;
	record->ToRec[17]=org->gnome->struct_change;
		*/

	record->ToRec[RECSIZE-1]=highest_fitness;
	}

 	 //push back novelty characterization
	 new_item->data.push_back(gather[0]);
 	 //set fitness (this is 'real' objective-based fitness, not novelty)
 	 new_item->fitness=fitness;
  return new_item;
}


//Perform evolution on single pole balacing, for gens generations
Population *maze_generational(char* outputdir,const char* mazefile,int param,const char *genes, int gens,bool novelty) 
{
    float archive_thresh=0.1;
    noveltyarchive archive(archive_thresh,*maze_novelty_metric,true,push_back_size);

    Population *pop;

    Genome *start_genome;
    char curword[20];
    int id;
    data_rec Record;

    ostringstream *fnamebuf;
    int gen;

    ifstream iFile(genes,ios::in);

    env=new Environment(mazefile);
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
  if(gen%200==0)
  {
  char filename[30];
  sprintf(filename,"%s_record.dat",output_dir);
  char fname[30];
  sprintf(fname,"%s_archive.dat",output_dir);
  archive.Serialize(fname);
  Record.serialize(filename);
  sprintf(fname,"%sgen_%d",output_dir,gen);
  pop->print_to_file_by_species(fname);
  }
  if (win)
  {
   char fname[30];
   sprintf(fname,"%s_wingen",output_dir);
    ofstream winfile(fname);
    winfile << gen << endl;
    break;
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
 
  int switch_amount = 20;

  //Evaluate each organism on a test
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
        
        data_record* newrec=new data_record();
	newrec->indiv_number=indiv_counter;
	//evaluate individual, get novelty point    
	(*curorg)->noveltypoint = maze_novelty_map((*curorg),newrec);
        (*curorg)->noveltypoint->indiv_number = indiv_counter;
        (*curorg)->datarec = newrec;	

	if((newrec->ToRec[3]>0.0 && newrec->ToRec[4]>0.0)) {
		(*curorg)->winner=true;
                win=true;	
	}
  
        if((*curorg)->noveltypoint->fitness > best_fitness)
	{
		best_fitness = (*curorg)->noveltypoint->fitness;
		cout << "NEW BEST " << best_fitness << endl;
	}	
        //add record of new indivdual to storage
	Record.add_new(newrec);
	indiv_counter++;
	
	//update fittest list
	archive.update_fittest(*curorg);
	
        if(!novelty)
    	   (*curorg)->fitness = (*curorg)->noveltypoint->fitness;
  }


  //adjust target every so often
  if(novelty)
  {
  if((generation)%switch_amount==0 && generation>0)
  {
     //merge populations together...
     //then make the measure_pop equal to the current population 
     Population *new_pop = archive.merge_populations(pop,measure_pop);
     cout << "Populations merged..." << endl;
     pop = new_pop;
  }
  
  if((generation % switch_amount)==0)
  {
    cout << "CREATE MEASUREPOP " << endl;
    measure_pop.clear();
    for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg)    {  
       measure_pop.push_back(new Organism(*(*curorg))); 
       //make it a copy so it dont get deleted
    }
    cout << "MEASUREPOP SIZE " << measure_pop.size() << endl;
  }
  
  }
  if(novelty)
  {
  
	//NEED TO CHANGE THESE TO GENERATIONAL EQUIVALENTS...
	//assign fitness scores based on novelty
 	archive.evaluate_population(pop,measure_pop,true);
	///now add to the archive (maybe remove & only add randomly?)
	archive.evaluate_population(pop,measure_pop,false);
        
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

