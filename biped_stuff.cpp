#include <ode/odeconfig.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>

#include "neat.h"
#include "organism.h"
#include "noveltyset.h"
#include "datarec.h"
#include <ode/ode.h>
#include "biped.h"
#include "experiments.h"
#include "calc_evol.h"

extern bool evaluate_switch;
extern bool seed_mode;
extern bool minimal_criteria;
extern bool population_dirty;
extern char output_dir[30];
extern char seed_name[40];

int novelty_function = NF_COGSAMPSQ;
vector<dGeomID> geoms;
vector<Creature*> creatures;
//NEAT + NS stuff
inline float dist(float x1, float y1, float x2, float y2)
{
    float xd = x1-x2;
    float yd = y1-y2;
    return xd*xd+yd*yd;
}

static	void calculate_delta(dVector3 v1, dVector3 v2, dVector3 o)
{
    for (int x=0; x<3; x++)
	o[x]=v2[x]-v1[x];
}

static	void calculate_power(dVector3 v,int pow)
{
    for (int x=0; x<3; x++)
    {
	float temp=v[x];
	bool sign=false;
	if (temp<0.0)
	    sign=true;
	for (int k=1; k<pow; k++)
	    v[x]*=temp;
	if (sign)
	    v[x]=(-v[x]);
    }
}

//novelty metric for maze simulation
float walker_novelty_metric(noveltyitem* x,noveltyitem* y)
{
    float dist=0.0;

    int size = x->data[0].size();
    int size2 = y->data[0].size();
    if(size!=size2) {
     cout << size << " " << size2 << endl;
     cout << "fucking horseshit." << endl;
     exit(0);
    }

    for (int k=0; k<size; k++)
    {
	float delta = x->data[0][k]-y->data[0][k];
	dist+=delta*delta;
    }

    return dist;

}

static dWorldID world;
static dSpaceID space;
static dGeomID floorplane;

static vector<dBodyID> bodies;

static dJointGroupID contactgroup;

// this is called by dSpaceCollide when two objects in space are
// potentially colliding.
static void nearCallback (void *data, dGeomID o1, dGeomID o2)
{
    dBodyID b1,b2;
    dBodyID test;
    assert(o1);
    assert(o2);

    b1 = dGeomGetBody(o1);
    b2 = dGeomGetBody(o2);

    if (b1 && b2 && dAreConnected (b1,b2)) return;

    if (o1 == floorplane || o2 == floorplane)
    {
        if (o1==floorplane)
            test=b2;
        if (o2==floorplane)
            test=b1;
//test should equal the body that is colliding with floor

        for (int x=0; x<creatures.size(); x++)
        {
            int bsize=creatures[x]->bodies.size();
            for (int y=0; y<bsize; y++)
                if (test==creatures[x]->bodies[y])
                    creatures[x]->onground[y]=true;
        }
    }

    const int N = 32;
    dContact contact[N];
    int n = dCollide (o1,o2,N,&(contact[0].geom),sizeof(dContact));
    if (n > 0)
    {
        for (int i=0; i<n; i++)
        {
            contact[i].surface.mode = 0;
            contact[i].surface.mu = dInfinity; //50.0; // was: dInfinity
            dJointID c = dJointCreateContact (world,contactgroup,&contact[i]);
            dJointAttach (c, dGeomGetBody(contact[i].geom.g1), dGeomGetBody(contact[i].geom.g2));
        }
    }
}

//executes one step of simulation 
void simulationStep(bool bMoviePlay)
{
    double timestep=0.01;
    if (!bMoviePlay)
    {
        dSpaceCollide (space,0,&nearCallback);
        dWorldStep(world,timestep);
    }

    for (int x=0; x<creatures.size(); x++)
        creatures[x]->Update(timestep);

    if (!bMoviePlay)
        dJointGroupEmpty (contactgroup);

}

void create_world(Controller* controller,bool log,bool bMoviePlay)
{
// create world
    dRandSetSeed(10);
    dInitODE();
    creatures.clear();
    world = dWorldCreate();
    space = dHashSpaceCreate (0);
    contactgroup = dJointGroupCreate (0);
    dWorldSetGravity (world,0,0,-9.8);
    floorplane = dCreatePlane (space,0,0,1, 0.0);
    dWorldSetERP(world,0.1);
    dWorldSetCFM(world,1E-4);

    Biped* biped = new Biped(log,bMoviePlay);
    dVector3 pos={0.0,0.0,0.0};

    biped->Create(world,space,pos,controller);
    creatures.push_back(biped);
}

void destroy_world()
{
    dJointGroupEmpty (contactgroup);
    dJointGroupDestroy (contactgroup);

    for (int x=0; x<geoms.size(); x++)
        dGeomDestroy(geoms[x]);

    for (int x=0; x<creatures.size(); x++)
    {
        creatures[x]->Destroy();
        delete creatures[x];
    }
    creatures.clear();
    bodies.clear();
    geoms.clear();

    dSpaceDestroy (space);
    dWorldDestroy (world);
    dCloseODE();
}

void update_behavior(vector<float> &k, Creature* c,bool good=true)
{

    if (novelty_function==NF_COGSAMPSQ)
    {
        dVector3& o_com= ((Biped*)c)->orig_com;
        dVector3& c_com= ((Biped*)c)->curr_com;
        dVector3 com;
        dVector3 delta;
        c->CenterOfMass(com);
        calculate_delta(o_com,com,delta);
        calculate_power(delta,2);

        if (good)
        {
            k.push_back(delta[0]);
            k.push_back(delta[1]);
        }
        else
        {
            k.push_back(0.0);
            k.push_back(0.0);
        }
    }
}

dReal evaluate_controller(Controller* controller,noveltyitem* ni,data_record* record,bool log)
{
    vector<float> k;
    dReal fitness;
    int timestep=0;
    const int simtime=1500;
    create_world(controller,log);
    while (!creatures[0]->abort() && timestep<simtime)
    {
        simulationStep();
        timestep++;
        if (timestep%100 == 0 && novelty_function % 2 == 1)
        {
            update_behavior(k,creatures[0]);
        }
        if (log && timestep%100==0)
            cout << creatures[0]->fitness() << endl;
    }
    int time=timestep;
    //for (int x=timestep+1; x<=simtime; x++)
    //    if (x%100==0)
    while(k.size()< (simtime/100*2))
            update_behavior(k,creatures[0]); 

    fitness=creatures[0]->fitness();
    ((Biped*)creatures[0])->lft.push_back(timestep);
    ((Biped*)creatures[0])->rft.push_back(timestep);
    if (ni!=NULL)
    {
        //ni->time=time;
        ni->novelty_scale = 1.0; 
        ni->data.push_back(k);
    }

    if (record!=NULL)
    {
        dVector3 com;
        creatures[0]->CenterOfMass(com);
        record->ToRec[0]=fitness;
        record->ToRec[1]=com[0];
        record->ToRec[2]=com[1];
        record->ToRec[3]=com[2];
        record->ToRec[4]=timestep;
    }

    destroy_world();
    return fitness;
}

noveltyitem* biped_evaluate(NEAT::Organism *org,data_record* data)
{
    noveltyitem *new_item = new noveltyitem;
    new_item->genotype=new Genome(*org->gnome);
    new_item->phenotype=new Network(*org->net);

    CTRNNController* cont = new CTRNNController(org->net);
    new_item->fitness=evaluate_controller(cont,new_item,data);
    if (new_item->fitness < 2.5) new_item->viable=false;
    else new_item->viable=true;
    delete cont;

    return new_item;
}


//novelty maze navigation run
Population *biped_novelty_realtime(char* outputdir,const char* mazefile,int par,const char* genes,bool novelty) {
	
    Population *pop;
    Genome *start_genome;
    char curword[20];

    int id;

	if(outputdir!=NULL) strcpy(output_dir,outputdir);
		
    if(!seed_mode)
        strcpy(seed_name,genes);
	//starter genes file
    ifstream iFile(seed_name,ios::in);
	
    cout<<"START BIPED NAVIGATOR NOVELTY REAL-TIME EVOLUTION VALIDATION"<<endl;
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
    if(evaluate_switch) { 
      int dist=0;
      double evol=0.0;
      evolvability_biped(pop->organisms[0],"dummyfile",&dist,&evol);
      cout << endl << dist << " " << evol << endl;
      return 0;
    }

    }
    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
    pop->set_evaluator(&biped_evaluate);     
    //pop->set_compatibility(&behavioral_compatibility); 
    //Start the evolution loop using rtNEAT method calls 
    biped_novelty_realtime_loop(pop,novelty);

    //clean up
    return pop;
}

//actual rtNEAT loop for novelty maze navigation runs
int biped_novelty_realtime_loop(Population *pop,bool novelty) {
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
  vector<Species*>::iterator curspec; //used in printing out debug info                                                         

  vector<Species*> sorted_species;  //Species sorted by max fit org in Species 

//was 1.0*number_of_samples+1.0 for earlier results...
   float archive_thresh=(1.0);// * 20.0 * envList.size(); //initial novelty threshold
  //if(!minimal_criteria)
  //	archive_thresh*=20;
 //if(constraint_switch)
   //archive_thresh/=200.0;
  cout << "Archive threshold: " << archive_thresh << endl;
  //archive of novel behaviors
  noveltyarchive archive(archive_thresh,*walker_novelty_metric,true,30,minimal_criteria);
	
  data_rec Record; //stores run information
	
  int count=0;
  int pause;

  //Real-time evolution variables                                                                                             
  int offspring_count;
  Organism *new_org;

  //We try to keep the number of species constant at this number                                                    
  int num_species_target=NEAT::pop_size/20;
  
  //This is where we determine the frequency of compatibility threshold adjustment
  int compat_adjust_frequency = NEAT::pop_size/20;
  if (compat_adjust_frequency < 1)
    compat_adjust_frequency = 1;

  //Initially, we evaluate the whole population                                                                               
  //Evaluate each organism on a test                   
  int indiv_counter=0; 
 pop->evaluate_all();

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
(offspring_count=0;offspring_count<NEAT::pop_size*2001;offspring_count++) 
{
//fix compat_threshold, so no speciation...
//      NEAT::compat_threshold = 1000000.0;
	//only continue past generation 1000 if not yet solved
	//if(offspring_count>=pop_size*1000 && firstflag)
        // if(firstflag)
	// break;

int evolveupdate=50000;
if(NEAT::evolvabilitytest && offspring_count % evolveupdate ==0) {
   char fn[100];
   sprintf(fn,"%s_evolvability%d.dat",output_dir,offspring_count/evolveupdate);
   for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
    evolvability_biped(*curorg,fn);
   }
}
	
	//end of generation
    if(offspring_count % (NEAT::pop_size*1) == 0)
	{
          /*
          if((offspring_count/NEAT::pop_size)%change_extinction_length==0)
            change_extinction_point();
          if((offspring_count/NEAT::pop_size)%change_goal_length==0)
            change_goal_location();
          */
          if(population_dirty) {
           pop->evaluate_all();
           population_dirty=false;
          }
                 if(novelty) {
			archive.end_of_gen_steady(pop);
			//archive.add_randomly(pop);
			archive.evaluate_population(pop,false);
			cout << "ARCHIVE SIZE:" << 
			archive.get_set_size() << endl;
                 }
                 cout << "GEN" << offspring_count/NEAT::pop_size << endl;
	 char fn[100];
         sprintf(fn,"%sdist%d",output_dir,offspring_count/NEAT::pop_size);
         if(NEAT::printdist)
         pop->print_distribution(fn);
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
	
/*	data_record* newrec=new data_record();
	newrec->indiv_number=indiv_counter;
	//evaluate individual, get novelty point
	new_org->noveltypoint = maze_novelty_map(new_org,newrec);
	new_org->noveltypoint->indiv_number = indiv_counter;
	new_org->fitness=new_org->noveltypoint->fitness;
*/
        data_record* newrec=new_org->datarec;
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
              //  cout << " :( " << endl;
	}	
        else
        {
            // cout << ":)" << new_org->noveltypoint->indiv_number << endl;
        }	
        //add record of new indivdual to storage
	//Record.add_new(newrec);
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

    //Remove the worst organism                                                                                               
    //if(rand_repl || fitness_measure ==fitness_rnd)
    // pop->remove_random();
    //else    
     pop->remove_worst();

  }
  
  //write out run information, archive, and final generation
  cout << "COMPLETED...";
  char filename[100];
  sprintf(filename,"%srecord.dat",output_dir);
  char fname[100];
  sprintf(fname,"%srtarchive.dat",output_dir);
  archive.Serialize(fname);
  //Record.serialize(filename);
  
  sprintf(fname,"%sfittest_final",output_dir);
  archive.serialize_fittest(fname);

  sprintf(fname,"%srtgen_final",output_dir);
  pop->print_to_file_by_species(fname);
  delete pop;
  exit(0);
  return 0;
}
 
#define BIPEDMUTATIONS 200
void evolvability_biped(Organism* org,char* fn,int* di,double* ev) {
 fstream file; 
 file.open(fn,ios::app|ios::out);
 cout <<"Evolvability..." << endl;
// file << "---" <<  " " << org->winner << endl;
 double points[BIPEDMUTATIONS*2];
 float minx=-3.0,maxx=3.0,miny=-3.0,maxy=3.0;
 double ox,oy,fit;
 int nodes;
 int connections; 
 data_record rec;
 for(int i=0;i<BIPEDMUTATIONS;i++) {
  Genome *new_gene= new Genome(*org->gnome); 
  //new_org->gnome = new Genome(*org->gnome);
  if(i!=0) //first copy is clean
   for(int j=0;j<1;j++) mutate_genome(new_gene);
  Organism *new_org= new Organism(0.0,new_gene,0);
  
  noveltyitem* nov_item = biped_evaluate(new_org,&rec);
  if(i==0) {
   fit=nov_item->fitness;
   nodes=new_org->net->nodecount();
   connections=new_org->net->linkcount();
   ox=rec.ToRec[1];
   oy=rec.ToRec[2];
  }
   //for(int k=0;k<nov_item->data[0].size();k++)
  //  file << nov_item->data[0][k] << " ";
  points[i*2]=(rec.ToRec[1]-minx)/(maxx-minx);
  points[i*2+1]=(rec.ToRec[2]-miny)/(maxy-miny);
  delete new_org;
  delete nov_item;
  //file << endl;
 }  
 int dist = distinct(points,BIPEDMUTATIONS);
 if(di!=NULL) *di=dist;
 double evol = test_indiv(points,BIPEDMUTATIONS);
 if(ev!=NULL) *ev=evol;
 file << dist << " " << evol << " " << ox << " " << oy << " " << nodes << " " <<connections << " " << fit << endl; 
 file.close();
}
