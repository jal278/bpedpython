#include "experiments.h"
#include "noveltyset.h"

#include "datarec.h"
#include "maze.h"

#include "histogram.h"
#include "calc_evol.h"
#include "genome.h"
//#define DEBUG_OUTPUT 1
#include <algorithm>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>
static Environment* env;
static vector<Environment*> envList;
static vector<Environment*> mcList;
static ofstream *logfile;
void evolvability(Organism* org,char* fn,int *a=NULL,double* b=NULL,bool recall=false);
using namespace std;
enum novelty_measure_type { novelty_sample, novelty_accum, novelty_sample_free };
static novelty_measure_type novelty_measure = novelty_sample;

enum fitness_measure_type { fitness_goal, fitness_drift, fitness_std,fitness_rnd,fitness_spin,fitness_changegoal,fitness_collisions };
static fitness_measure_type fitness_measure = fitness_goal;

static bool mc_reach_onepoint=true;

bool population_dirty=false;

static bool extinction=true;
void set_extinction(bool _ext) {
    extinction=_ext;
}
static Point extinction_point(0.0,0.0);
static float extinction_radius=50.0;
static int change_extinction_length=5;

void change_extinction_point() {
    float minx,maxx,miny,maxy;
    envList[0]->get_range(minx,miny,maxx,maxy);

    extinction_point.x = randfloat()*(maxx-minx)+minx;
    extinction_point.y = randfloat()*(maxy-miny)+miny;
    if (extinction)
        population_dirty=true;
}
bool extinct(Point k) {
    if (k.distance(extinction_point)<extinction_radius)
        return true;
    return false;
}


static Point changing_goal(0.0,0.0);
static int change_goal_length=5;


void change_goal_location() {
    static int count=0;
    static float vx=0.0;
    static float vy=0.0;
    static float newx;
    static float newy;
    float minx,maxx,miny,maxy;
    envList[0]->get_range(minx,miny,maxx,maxy);

    if (count==0) {
        newx = randfloat()*(maxx-minx)+minx;
        newy = randfloat()*(maxy-miny)+miny;
        changing_goal.x=newx;
        changing_goal.y=newy;
    }

    if (count%10==0) {
        newx = randfloat()*(maxx-minx)+minx;
        newy = randfloat()*(maxy-miny)+miny;
        vx=(newx-changing_goal.x)/10.0;
        vy=(newy-changing_goal.y)/10.0;
    }

    count++;

    changing_goal.x+=vx;
    changing_goal.y+=vy;
    if (changing_goal.x<minx) changing_goal.x=minx;
    else if (changing_goal.x>maxx) changing_goal.x=maxx;
    if (changing_goal.y<miny) changing_goal.y=miny;
    else if (changing_goal.y>maxy) changing_goal.y=maxy;


    if (fitness_measure==fitness_changegoal)
    {
        cout << "New goal: " << changing_goal.x << " " << changing_goal.y << endl;
        population_dirty=true;
    }

}




static int number_of_samples = 1;
static int simulated_timesteps = 400;
bool seed_mode = false;
char seed_name[100]="";
static char mc_mazes[40]="";
bool minimal_criteria=false;
bool evaluate_switch=false;
static bool goal_attract=true;

static bool activity_stats=false;

static bool constraint_switch=false;
static bool area_of_interest=false;
static bool rand_repl=false;
void set_evaluate(bool val) {
    evaluate_switch=val;
}
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
    if (strlen(seed_name)>0)
        seed_mode=true;
}

void set_fit_measure(string m)
{
    if (m=="rnd")
        fitness_measure=fitness_rnd;
    if (m=="std")
        fitness_measure=fitness_std;
    if (m=="drift")
        fitness_measure=fitness_drift;
    if (m=="goal")
        fitness_measure=fitness_goal;
    if (m=="spin")
        fitness_measure=fitness_spin;
    if (m=="changegoal")
        fitness_measure=fitness_changegoal;
    if (m=="collisions")
        fitness_measure=fitness_collisions;
    cout << "Fitness measure " << fitness_measure << endl;
}

void set_nov_measure(string m)
{
    if (m=="std" || m=="sample")
        novelty_measure=novelty_sample;
    if (m=="accum")
        novelty_measure=novelty_accum;
    if (m=="sample_free")
        novelty_measure=novelty_sample_free;
    cout << "Novelty measure " << novelty_measure << endl;
}

char output_dir[30]="";


static int param=-1;
static int push_back_size = 20;

//used for discretization studies
double discretize(double x,long bins,double low, double high)
{
    double norm = x-low;
    double binsize = (high-low)/bins;
    int bin = (int)(norm/binsize);
    if (bin==bins)
        bin--;
    double result = (double)binsize*bin+binsize/2.0+low;
    return result;
}

long powerof2(int num)
{
    long x=1;
    if (num==0) return 1;
    for (int y=0; y<num; y++)
        x*=2;
    return x;
}

//novelty metric for maze simulation
float maze_novelty_metric(noveltyitem* x,noveltyitem* y)
{
    float diff = 0.0;
    for (int k=0; k<(int)x->data.size(); k++)
    {
        diff+=hist_diff(x->data[k],y->data[k]);
    }
    return diff;
}

static void read_in_environments(const char* mazefile, vector<Environment*>& envLst)
{
    ifstream listfile(mazefile);

    while (!listfile.eof())
    {
        string filename;
        getline(listfile,filename);
        if (filename.length() == 0)
            break;
        cout << "Reading maze: " << filename << endl;
        Environment* new_env = new Environment(filename.c_str());
        envLst.push_back(new_env);
    }

}

#define MAX_NICHES 300000
class passive_niche {
	public:
        
	//niches, whether niches have been explored
	vector<Organism*> niches[MAX_NICHES];
	bool explored[MAX_NICHES];
	int order[MAX_NICHES];
	//params
	int niche_size;
        int evals;
	int density;
	int nc;
	bool firstsolved;
	bool random;
        float minx,miny,maxx,maxy;

        void calc_evolvability(char*fn) {
	    for(int i=0;i<nc;i++) {
		cout << "evolvability niche " << i << endl;
		for(int j=0;j<5;j++) {	
			int ns=niches[order[i]].size();
		    Organism *org = niches[order[i]][randint(0,ns-1)];
	            evolvability(org,fn);
		}
	    }
 	}	

	passive_niche(bool _r=false) {
         random=_r;
	 density=30;
 	 niche_size=10;
         evals=100001;
	 for(int i=0;i<MAX_NICHES;i++) explored[i]=false;
	 for(int i=0;i<MAX_NICHES;i++) order[i]=false;
         nc=0;
	 firstsolved=true;
       	}


        void print_niches() {
	for(int x=0;x<density;x++) 
	 {
	  for(int y=0;y<density;y++)
	   {
		cout << explored[x*density+y];
           }
		cout <<endl;
	 }
	}

	int scale(int d,float val, float min,float max) {
		return (int)(d*(val-min)/((max+0.01f)-min));
	}

	int map_into_niche(Organism* o) {
 	 float x= o->noveltypoint->data[0][0];
	 float y= o->noveltypoint->data[0][1];
         if (!random)
          return ((int)x)/10*300+((int)y)/10; //to match other simulations
         else
          return rand()%400;
	 //return scale(density,x,minx,maxx)*density+scale(density,y,miny,maxy);
        }

	void insert_population(Population* pop) {
	 vector<Organism*>::iterator curorg;
         for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
	  insert_org((*curorg)); 
	 }
	}

	void insert_org(Organism* org) {
	 int target_niche = map_into_niche(org);
	 if(target_niche<0)
		return; 
	int sz=niches[target_niche].size();
	 if(sz==0) {
		explored[target_niche]=true;
		order[nc]=target_niche;
		nc++;
	 }
	 if(sz>=niche_size) {
		//return;
		remove_one_from_niche(target_niche);
	 }
	 niches[target_niche].push_back(org);
	}

        void remove_one_from_niche(int n) {
		//cout << "deleting from " << n << endl;

                int to_rem = randint(0,niches[n].size()-1);
		Organism* o = niches[n][to_rem];
		niches[n].erase(niches[n].begin()+to_rem);
		delete o;
	}
    int exploredcount() {
	int c=0;
	for(int i=0;i<MAX_NICHES;i++) if (explored[i]) c++;
	return  c;
	}
    void run_niche(Population* initpop) {
    //num evals
    envList[0]->get_range(minx,miny,maxx,maxy);
	 insert_population(initpop);
    int e=0;
    int num_niches=density*density;
    while(e<evals) {
	cout << "evals " << e <<endl;
	cout << "explored " << exploredcount() << endl;
	vector<Organism*> children;
	//print_niches();
	int conn=0;
	int nodes=0;
	int count=0;
	for(int i=0;i<num_niches;i++) {
		if(niches[i].size()==0)
			continue;
		Genome *new_gene= new Genome(*niches[i][randint(0,niches[i].size()-1)]->gnome);
		//for(vector<Organism*>::iterator it=niches[i].begin();it!=niches[i].end();it++)
		//{
                //Genome *new_gene= new Genome(*(*it)->gnome);
		mutate_genome(new_gene,true);
		Organism* new_org= new Organism(0.0,new_gene,0);
		initpop->evaluate_organism(new_org);
		if(new_org->datarec->ToRec[3] > 0 && firstsolved) {
			cout << "solved " << e << endl;
			firstsolved=false;
		}
		nodes+=new_org->net->nodecount();
		conn+=new_org->net->linkcount();
		count++;
		children.push_back(new_org);
		e++;
		int upcnt=10000;
		if(e%upcnt==0) {
			char fn[100];
        sprintf(fn,"%s_evolvability%d.dat",output_dir,e/upcnt);
			calc_evolvability(fn);
		}
        	//}
          }
	   cout << "avgnodes" << ((float)nodes)/count << endl;
	   cout << "avgconns" << ((float)conn)/count << endl;
	
		for(vector<Organism*>::iterator it=children.begin();it!=children.end();it++)
		insert_org(*it);

     }
  }

};

void enumerate_behaviors(const char* mazefile, long long par,const char* outfile,int count) {
read_in_environments(mazefile,envList);

ofstream ofile(outfile);
ofile << par << endl;

for(int x=0;x<count;x++) {
Genome *g = new Genome(3,2,2,2);

long long partemp=par;
for(int i=17;i>=0;i--) {
long long val = (partemp % 3) - 1;
g->genes[i]->lnk->weight = (double)val;
partemp /= 3; 
}

Organism* new_org= new Organism(0.0,g,0);
noveltyitem* nov_item = maze_novelty_map(new_org);
ofile << nov_item->data[0][0] << " " << nov_item->data[0][1] << endl;
delete nov_item;
delete new_org;
par++;
}
}

//passive algorithm
Population *maze_passive(char* outputdir,const char* mazefile,int par,const char* genes,bool novelty) {

    Population *pop;
    Genome *start_genome;
    char curword[20];

    int id;

    read_in_environments(mazefile,envList);

    if (strlen(mc_mazes)>0)
        read_in_environments(mc_mazes,mcList);

    push_back_size=par;
    if (outputdir!=NULL) strcpy(output_dir,outputdir);

    if (!seed_mode)
        strcpy(seed_name,genes);
    //starter genes file
    ifstream iFile(seed_name,ios::in);

    cout<<"START MAZE NAVIGATOR NOVELTY REAL-TIME EVOLUTION VALIDATION"<<endl;
    if (!seed_mode)
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
    if (!seed_mode)
        pop=new Population(start_genome,NEAT::pop_size);
    else
    {
        pop=new Population(seed_name);//start_genome,NEAT::pop_size,0.0);
        if (evaluate_switch) {
            int dist=0;
            double evol=0.0;
            evolvability(pop->organisms[0],"dummyfile",&dist,&evol,true);
            cout << endl << dist << " " << evol << endl;
            return 0;
        }

    }
    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
    pop->set_evaluator(&maze_novelty_map);
    pop->evaluate_all();
 
    passive_niche pn(novelty);
    pn.run_niche(pop);

    //pop->set_compatibility(&behavioral_compatibility);
    //Start the evolution loop using rtNEAT method calls




    //clean up
    delete env;
    return pop;
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

    if (strlen(mc_mazes)>0)
        read_in_environments(mc_mazes,mcList);
    //param=par;
    push_back_size=par;
    if (outputdir!=NULL) strcpy(output_dir,outputdir);

    if (!seed_mode)
        strcpy(seed_name,genes);
    //starter genes file
    ifstream iFile(seed_name,ios::in);

    cout<<"START MAZE NAVIGATOR NOVELTY REAL-TIME EVOLUTION VALIDATION"<<endl;
    if (!seed_mode)
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
    if (!seed_mode)
        pop=new Population(start_genome,NEAT::pop_size);
    else
    {
        pop=new Population(seed_name);//start_genome,NEAT::pop_size,0.0);
        if (evaluate_switch) {
            int dist=0;
            double evol=0.0;
            evolvability(pop->organisms[0],"dummyfile",&dist,&evol,true);
            cout << endl << dist << " " << evol << endl;
            return 0;
        }

    }
    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
    pop->set_evaluator(&maze_novelty_map);
    //pop->set_compatibility(&behavioral_compatibility);
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
    int num_species_target=NEAT::pop_size/20;

    //This is where we determine the frequency of compatibility threshold adjustment
    int compat_adjust_frequency = NEAT::pop_size/20;
    if (compat_adjust_frequency < 1)
        compat_adjust_frequency = 1;

    char sol_evo_fn[100];
    sprintf(sol_evo_fn,"%s_solution_evolvability.dat",output_dir);

//activity stat log file
    char asfn[100];
    sprintf(asfn,"%s_activitystats.dat",output_dir);
    ofstream activity_stat_file(asfn);
    if (activity_stats)
        reset_activity();
    //Initially, we evaluate the whole population
    //Evaluate each organism on a test
    int indiv_counter=0;
    pop->evaluate_all();

    if (novelty) {
        //assign fitness scores based on novelty
        archive.evaluate_population(pop,true);
        //add to archive
        archive.evaluate_population(pop,false);
    }

    if (novelty && minimal_criteria)
        for (curorg=(pop->organisms).begin(); curorg!=(pop->organisms).end(); ++curorg)
        {
            (*curorg)->fitness = SNUM/1000.0;
        }
//Get ready for real-time loop
    //Rank all the organisms from best to worst in each species
    pop->rank_within_species();

    //This average must be kept up-to-date by rtNEAT in order to select species probabailistically for reproduction
    pop->estimate_all_averages();

    cout <<"Entering real time loop..." << endl;

    //Now create offspring one at a time, testing each offspring,
    // and replacing the worst with the new offspring if its better
    for
    (offspring_count=0; offspring_count<NEAT::pop_size*2001; offspring_count++)
    {
//fix compat_threshold, so no speciation...
//      NEAT::compat_threshold = 1000000.0;
        //only continue past generation 1000 if not yet solved
        //if(offspring_count>=pop_size*1000 && firstflag)
        // if(firstflag)
        // break;

        int evolveupdate=6250;
        if (NEAT::evolvabilitytest && offspring_count % evolveupdate ==0) {
            char fn[100];
            sprintf(fn,"%s_evolvability%d.dat",output_dir,offspring_count/evolveupdate);
            for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
                evolvability(*curorg,fn);
            }
        }

        if (activity_stats&& offspring_count % 10000 == 0)
        {
            pop->update_statistics();
            activity_stat_file << offspring_count << " "  << calculate_diversity()
            << " " << calculate_cumulative_activity() << " " <<
            calculate_average_activity() << endl;
        }

        //end of generation
        if (offspring_count % (NEAT::pop_size*1) == 0)
        {
            if ((offspring_count/NEAT::pop_size)%change_extinction_length==0)
                change_extinction_point();
            if ((offspring_count/NEAT::pop_size)%change_goal_length==0)
                change_goal_location();

            if (population_dirty) {
                pop->evaluate_all();
                population_dirty=false;
            }
            if (novelty) {
                archive.end_of_gen_steady(pop);
                //archive.add_randomly(pop);
                archive.evaluate_population(pop,false);
                cout << "ARCHIVE SIZE:" <<
                     archive.get_set_size() << endl;
            }
        double mx=0.0;  
double tot=0.0;    
    Organism* b;  
	for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
tot+=(*curorg)->noveltypoint->fitness;
if( (*curorg)->noveltypoint->fitness > mx) {
mx=(*curorg)->noveltypoint->fitness; b=(*curorg); }
} 
           cout << "GEN" << offspring_count/NEAT::pop_size << " " << tot << " " << mx <<  endl;
            char fn[100];
            sprintf(fn,"%sdist%d",output_dir,offspring_count/NEAT::pop_size);
            if (NEAT::printdist)
                pop->print_distribution(fn);
        }

        //write out current generation and fittest individuals
        if ( offspring_count % (NEAT::pop_size*NEAT::print_every) == 0 )
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
            if (novelty) {
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
        for (curspec=(pop->species).begin(); curspec!=(pop->species).end(); curspec++) {
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
        if (novelty) {
            archive.evaluate_individual(new_org,pop->organisms);
            //newrec->ToRec[5] = archive.get_threshold();
            newrec->ToRec[6] = archive.get_set_size();
            newrec->ToRec[RECSIZE-2] = new_org->noveltypoint->novelty;
        }
        if ( !new_org->noveltypoint->viable && minimal_criteria)
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

        //if solution, do evolvability
        
        if (newrec->ToRec[3]>=envList.size() && newrec->ToRec[4]>=envList.size()) {
            if (false) //NEAT::evolvabilitytest)
            {
                cout << "solution found" << endl;
                evolvability(new_org,sol_evo_fn);
            }
        }
        

        if (!weakfirst && (newrec->ToRec[3]>=envList.size())) {
            weakfirst=true;
            //NEAT::evolvabilitytest=true; //TODO REMOVE LATER
            char filename[100];
            sprintf(filename,"%srtgen_weakfirst",output_dir);
            pop->print_to_file_by_species(filename);
            cout << "Maze weakly solved by indiv# " << indiv_counter << endl;
//disable quit for now
            if (fitness_measure == fitness_goal && false)
                firstflag=true;
        }
        //write out the first individual to solve maze
        if (!firstflag && (newrec->ToRec[3]>=envList.size() && newrec->ToRec[4]>=envList.size())) {
            firstflag=true;
            char filename[100];
            sprintf(filename,"%srtgen_first",output_dir);
            pop->print_to_file_by_species(filename);
            cout << "Maze solved by indiv# " << indiv_counter << endl;
            //break;
        }

        //Remove the worst organism
        if (rand_repl || fitness_measure ==fitness_rnd)
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
    //Record.serialize(filename);

    sprintf(fname,"%sfittest_final",output_dir);
    archive.serialize_fittest(fname);

    sprintf(fname,"%srtgen_final",output_dir);
    pop->print_to_file_by_species(fname);
    delete pop;
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
    for (int i=0; i<10; i++)
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
    if (dist<=1) dist=1;
    double fitness = 5.0/dist; //used for accumulated fitness (obselete)

    return fitness;
}
double mazesim(Network* net, vector< vector<float> > &dc, data_record *record,Environment* the_env,Organism* o=NULL)
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

    if (novelty_measure == novelty_sample ||
            novelty_measure ==novelty_sample_free)
        data.reserve(timesteps/stepsize);
    if (novelty_measure == novelty_accum)
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
    if (number_of_samples>0)
        stepsize=timesteps/number_of_samples;

    for (int i=0; i<timesteps; i++)
    {
        fitness+=mazesimStep(newenv,net,dc);
        //if taking additional samples, collect during run
        if (novelty_measure==novelty_sample ||
                novelty_measure==novelty_sample_free)
            if ((timesteps-i-1)%stepsize==0)
            {
               if(!newenv->hero.collide) {
                data.push_back(newenv->hero.location.x);
                data.push_back(newenv->hero.location.y);
               }
               else {
                data.push_back(-10.0);
                data.push_back(-10.0);
               }
            }

        float loc[2]={newenv->hero.location.x,newenv->hero.location.y};
        if (novelty_measure==novelty_accum)
        {
            accum->add_point(loc);
        }
    }
    if (extinction) {
        if (extinct(newenv->hero.location)) {
            if (o!=NULL) o->eliminate=true;
        }
    }
    //calculate fitness of individual as closeness to target
    if (fitness_measure == fitness_goal)
    {
        fitness=300.0 - newenv->distance_to_target(); //was 500 for MCNS
        if (fitness<0.1) fitness=0.1;
        //if (newenv->hero.collide)
        //	fitness+=50;
    }
    if (fitness_measure == fitness_collisions)
    {
        fitness= (-newenv->hero.collisions);
    }
    if (fitness_measure == fitness_spin)
    {
        fitness=log(newenv->hero.total_spin+0.1);
        if (fitness>7) fitness=7.0;
        if (fitness<0) fitness=0.0;
        fitness=7.01-fitness;
    }

    if (fitness_measure == fitness_changegoal)
    {
        fitness=500-changing_goal.distance(newenv->hero.location);
        if (fitness<0.1) fitness=0.1;
    }

    if (fitness_measure ==fitness_rnd)
    {
        //todo assign random fitness, this needs to get reassigned
        //often...
        fitness = randint(10,100);
    }

    //calculate fitness as meeting minimal criteria
    if (fitness_measure == fitness_drift)
    {
        if (newenv->reachgoal)
        {
            fitness=1.0;
            if (newenv->reachpoi)
                fitness=500.0;
        }
        else
        {
            fitness=SNUM/1000.0;
        }
    }

    if (fitness_measure == fitness_std)
    {
        fitness=SNUM;
        float mod = 500.0 - newenv->closest_to_target;
        if (mod<0) mod=0.0;
        float mod2 = 500.0 - newenv->closest_to_poi;
        if (mod2<0) mod2=0.0;

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

    if (novelty_measure==novelty_sample || novelty_measure==novelty_sample_free)
        if (false)
        {
            //novelty point is the ending location of the navigator
            data.push_back(x);
            data.push_back(y);
        }

    if (novelty_measure==novelty_accum)
    {
        accum->transform();
        for (int x=0; x<accum->size; x++)
            data.push_back(accum->buffer[x]);
    }

    if (record!=NULL)
    {
        record->ToRec[0]+=fitness;
        record->ToRec[1]=newenv->hero.location.x;
        record->ToRec[2]=newenv->hero.location.y;
        record->ToRec[3]+=newenv->reachgoal;
        record->ToRec[4]+=newenv->reachpoi;
        record->ToRec[5]= (-newenv->hero.collisions);
    }

    if (novelty_measure==novelty_accum)
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

void mutate_genome(Genome* new_genome,bool traits)
{
    Network* net_analogue;
    double mut_power=NEAT::weight_mut_power;
				double inno=0;
				int id=0;
    new_genome->mutate_link_weights(mut_power,1.0,GAUSSIAN);
    if(traits) {
	vector<Innovation*> innos;
    if (randfloat()<NEAT::mutate_add_node_prob) 
				new_genome->mutate_add_node(innos,id,inno);
    else if (randfloat()<NEAT::mutate_add_link_prob) {
				//cout<<"mutate add link"<<endl;
				net_analogue=new_genome->genesis(0);
				new_genome->mutate_add_link(innos,inno,NEAT::newlink_tries);
				delete net_analogue;
    }


	if(randfloat()<0.05)
	new_genome->mutate_random_trait();
	if(randfloat()<0.05)
	new_genome->mutate_link_trait(1);
    }

    return;
}

#define MUTATIONS 200
void evolvability(Organism* org,char* fn,int* di,double* ev,bool recall) {
    bool solution=false;
    fstream file;
    file.open(fn,ios::app|ios::out);
    cout <<"Evolvability..." << endl;
// file << "---" <<  " " << org->winner << endl;
    double points[MUTATIONS*2];
    float minx,maxx,miny,maxy;
    envList[0]->get_range(minx,miny,maxx,maxy);
    double ox,oy,fit;
    int nodes;
    int connections;
    for (int i=0; i<MUTATIONS; i++) {
        Genome *new_gene= new Genome(*org->gnome);
        //new_org->gnome = new Genome(*org->gnome);
        if (i!=0) //first copy is clean
            for (int j=0; j<1; j++) mutate_genome(new_gene);
        Organism *new_org= new Organism(0.0,new_gene,0);

        noveltyitem* nov_item = maze_novelty_map(new_org);
        if (i==0) {
            fit=nov_item->fitness;
            nodes=new_org->net->nodecount();
            connections=new_org->net->linkcount();
            ox=nov_item->data[0][0];
            oy=nov_item->data[0][1];
        }
        if (nov_item->fitness>340) solution=true;
        if(recall) {
         for(int k=0;k<nov_item->data[0].size();k++)
          file << nov_item->data[0][k] << " ";
         file << endl;
        }
        points[i*2]=nov_item->data[0][0];
        points[i*2+1]=nov_item->data[0][1];
        //HOW IT WAS:
        //points[i*2]=(nov_item->data[0][0]-minx)/(maxx-minx);
        //points[i*2+1]=(nov_item->data[0][1]-miny)/(maxy-miny);
        delete new_org;
        delete nov_item;
        //file << endl;
    }
    int dist = distinct(points,MUTATIONS,2);
    if (di!=NULL) *di=dist;
    double evol = 0.0; // test_indiv(points,MUTATIONS);
    if (ev!=NULL) *ev=evol;
    if(!recall)
    file << dist << " " << evol << " " << ox << " " << oy << " " << nodes << " " <<connections << " " << fit << " " << solution << endl;
    file.close();
    return;
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

    if (record!=NULL && minimal_criteria)
        for (int x=0; x<mcList.size(); x++)
        {
            record->ToRec[3]=0;
            mazesim(org->net,gather,record,mcList[x]);
            if (!record->ToRec[3]) {
                new_item->viable=false;
                // cout << "not viable..." << endl;
                break;
            }
            //cout << "viable..." << endl;
        }

    gather.clear();
    if (record!=NULL) {
        record->ToRec[0]=0;
        record->ToRec[3]=0;
        record->ToRec[4]=0;
    }

    if (true) //new_item->viable)
        for (int x=0; x<envList.size(); x++)
        {
            if (record!=NULL) {
                c1old = record->ToRec[3];
                c2old = record->ToRec[4];
            }

            org->eliminate=false;
            fitness+=mazesim(org->net,gather,record,envList[x],org);
             if (org->eliminate) {
                new_item->viable=false;
                org->eliminate=false;
            }

            if (record!=NULL) {
                constraint_vector.push_back(record->ToRec[3]-c1old);
                constraint_vector.push_back(record->ToRec[4]-c2old);
                c1old=record->ToRec[3];
                new_item->secondary=record->ToRec[5];
                //if(record->ToRec[5]==1)
                //  new_item->viable=false;
            }
            else {
                constraint_vector.push_back(0);
                constraint_vector.push_back(0);
            }
        }
    else
    {
        for (int x=0; x<envList.size(); x++) constraint_vector.push_back(0);
    }

    //minimal criteria must be met in *all* scenarios...

    if (record!=NULL)
    {
        if (area_of_interest)
        {
            if (!contained(record->ToRec[1],record->ToRec[2]))
            {
                new_item->viable=false;
                //cout << "not viable..." << endl;
            }
            //else cout << "viable..." << endl;
        }
        else if (mc_reach_onepoint)
            if ( record->ToRec[3]<envList.size())
            {
                new_item->viable=false;
                //cout << record->ToRec[3] << endl;
            }
            else {
                //cout << "viable... " << endl;
            }
    }

    if (fitness>highest_fitness)
        highest_fitness=fitness;

    //keep track of highest fitness so far in record
    if (record!=NULL)
    {
        record->ToRec[5]=new_item->viable;
        if (record->ToRec[3]>best)
        {
            best=record->ToRec[3];
            cout << "best: " << best << endl;
        }
        record->ToRec[RECSIZE-1]=highest_fitness;
    }

    //push back novelty characterization
    if (!remove_regular)
        for (int i=0; i<gather.size(); i++)
            new_item->data.push_back(gather[i]);

    if (apply_constraints)
        new_item->data.push_back(constraint_vector);
    //set fitness (this is 'real' objective-based fitness, not novelty)
    new_item->fitness=fitness;
    org->fitness=fitness;

    return new_item;
}






static int maxgens;
//Perform evolution on single pole balacing, for gens generations
Population *maze_generational(char* outputdir,const char* mazefile,int param,const char *genes, int gens,bool novelty)
{

    char logname[100];
    sprintf(logname,"%s_log.txt",outputdir);
    logfile=new ofstream(logname);

    maxgens=gens;
    float archive_thresh=3.0;

    noveltyarchive archive(archive_thresh,*maze_novelty_metric,true,push_back_size,minimal_criteria,true);

//if doing multiobjective, turn off speciation, TODO:maybe turn off elitism
    if (NEAT::multiobjective) NEAT::speciation=false;

    Population *pop;

    Genome *start_genome;
    char curword[20];
    int id;
    data_rec Record;

    ostringstream *fnamebuf;
    int gen;

    if (!seed_mode)
        strcpy(seed_name,genes);
    if(seed_mode)
	cout << "READING IN SEED" << endl;
    ifstream iFile(seed_name,ios::in);

    read_in_environments(mazefile,envList);

    if (strlen(mc_mazes)>0)
        read_in_environments(mc_mazes,mcList);

    if (outputdir!=NULL) strcpy(output_dir,outputdir);
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
    cout << "Start genomes node: " << start_genome->nodes.size() << endl;
    pop=new Population(start_genome,NEAT::pop_size);

    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
   
//set evaluator
    pop->set_evaluator(&maze_novelty_map);
//pop->set_compatibility(&behavioral_compatibility);    
for (gen=0; gen<=maxgens; gen++)  { //WAS 1000
        cout<<"Generation "<<gen<<endl;
        if (gen%change_extinction_length==0)
            change_extinction_point();
        if (gen%change_goal_length==0)
            change_goal_location();
        bool win = maze_generational_epoch(&pop,gen,Record,archive,novelty);


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

    delete logfile;
    delete pop;
    return pop;

}

int maze_generational_epoch(Population **pop2,int generation,data_rec& Record, noveltyarchive& archive, bool novelty) {
    Population* pop= *pop2;
    vector<Organism*>::iterator curorg;
    vector<Species*>::iterator curspecies;
    static double best_fitness =0.0;
    static double best_secondary =  -100000.0;
    static vector<Organism*> measure_pop;
//char cfilename[100];
//strncpy( cfilename, filename.c_str(), 100 );

//ofstream cfilename(filename.c_str());

    static bool win=false;
    static bool firstflag=false;
    static bool weakfirst=false;
    int winnernum;
    int indiv_counter=0;

    int evolveupdate=100;
   if (generation==0) pop->evaluate_all();

    if (NEAT::multiobjective) {  //merge and filter population
        for (curorg=(pop->organisms).begin(); curorg!=(pop->organisms).end(); ++curorg) {
            measure_pop.push_back(new Organism(*(*curorg),true)); //TODO:maybe make a copy?
        }

//evaluate this 'super-population'
        archive.rank(measure_pop);
        /*
        for(int i=0;i<measure_pop.size();i++) {
        cout << measure_pop[i]->noveltypoint->competition << " " << measure_pop[i]->noveltypoint->novelty << " " <<  measure_pop[i]->noveltypoint->rank << endl;
        }
        */
        if (generation!=0) {
//chop population down by half (maybe delete orgs that aren't used)
            int start=measure_pop.size()/2;
            vector<Organism*>::iterator it;
            for (it=measure_pop.begin()+start; it!=measure_pop.end(); it++)
                delete (*it);
            measure_pop.erase(measure_pop.begin()+(measure_pop.size()/2),measure_pop.end());
        }
//delete old pop, create new pop
        Genome* sg=pop->start_genome;
	delete pop;
        pop=new Population(measure_pop);
        pop->start_genome=sg;
        pop->set_evaluator(&maze_novelty_map);
        *pop2= pop;
    }

    if (NEAT::evolvabilitytest && generation%evolveupdate==0)
    {
        char fn[100];
        sprintf(fn,"%s_evolvability%d.dat",output_dir,generation/evolveupdate);
        for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
            evolvability(*curorg,fn);
        }
    }

//Evaluate each organism on a test
    float divtotal=0.0;
    for (curorg=(pop->organisms).begin(); curorg!=(pop->organisms).end(); ++curorg) {

//newrec->indiv_number=indiv_counter;
//data_record* newrec=new data_record();
//evaluate individual, get novelty point
//(*curorg)->noveltypoint = maze_novelty_map((*curorg),newrec);
//(*curorg)->noveltypoint->indiv_number = indiv_counter;
//(*curorg)->datarec = newrec;
        data_record* newrec = (*curorg)->datarec;
        if (!weakfirst && (newrec->ToRec[3]>=envList.size())) {
            weakfirst=true;
            char filename[100];
            cout << "Maze weakly solved by indiv# " << indiv_counter << endl;
//disable quit for now
        }
        //write out the first individual to solve maze
        if (!firstflag && (newrec->ToRec[3]>=envList.size() && newrec->ToRec[4]>=envList.size())) {
            firstflag=true;
            char filename[100];
            cout << "Maze solved by indiv# " << indiv_counter << endl;
            //break;
        }

        if ((newrec->ToRec[3]>=envList.size())) { // && newrec->ToRec[4]>0.0)) {
            if ((newrec->ToRec[4]>=envList.size())) {
                (*curorg)->winner=true;
                win=true;
            }
            if (fitness_measure==fitness_goal) {
                (*curorg)->winner=true;
                win=true;
                if ((*curorg)->noveltypoint->secondary >best_secondary) {
                    best_secondary=(*curorg)->noveltypoint->secondary;
                    cout << "NEW BEST SEC " << best_secondary << endl;

                }
            }
        }
        divtotal+=(*curorg)->noveltypoint->genodiv;
        if ((*curorg)->noveltypoint->fitness > best_fitness)
        {
            best_fitness = (*curorg)->noveltypoint->fitness;
            cout << "NEW BEST " << best_fitness << endl;
        }


//add record of new indivdual to storage
//TODO: PUT BACK IN (to fix record.dat...)
//Record.add_new(newrec);
        indiv_counter++;
        /*
        if((*curorg)->noveltypoint->viable)
             cout << "viable..." << endl;
        else
             cout << "not viable..." << endl;
        */
        if ( !(*curorg)->noveltypoint->viable && minimal_criteria)
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
        if (!novelty)
            (*curorg)->fitness = (*curorg)->noveltypoint->fitness;
    }
    cout << "DIVTOTAL:" << divtotal << endl;

    (*logfile) << generation << " " << best_fitness << " " << best_secondary << endl;
    
    if (novelty)
    {

//NEED TO CHANGE THESE TO GENERATIONAL EQUIVALENTS...
//assign fitness scores based on novelty
        archive.evaluate_population(pop,true);
///now add to the archive (maybe remove & only add randomly?)
        archive.evaluate_population(pop,false);

        if (NEAT::multiobjective)
            archive.rank(pop->organisms);

        for (curorg=(pop->organisms).begin(); curorg!=(pop->organisms).end(); ++curorg) {
            if ( !(*curorg)->noveltypoint->viable && minimal_criteria)
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
    if (NEAT::printdist)
        pop->print_distribution(fn);
//Average and max their fitnesses for dumping to file and snapshot
    for (curspecies=(pop->species).begin(); curspecies!=(pop->species).end(); ++curspecies) {

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
    if (win && !firstflag)
    {
        for (curorg=(pop->organisms).begin(); curorg!=(pop->organisms).end(); ++curorg) {
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

//writing out stuff
    if (generation%NEAT::print_every == 0 )
    {
        char filename[100];
        sprintf(filename,"%s_record.dat",output_dir);
        char fname[100];
        sprintf(fname,"%s_archive.dat",output_dir);
        archive.Serialize(fname);
//Record.serialize(filename);
        sprintf(fname,"%sgen_%d",output_dir,generation);
        pop->print_to_file_by_species(fname);
    }

    if (NEAT::multiobjective) {
        for (curorg=measure_pop.begin(); curorg!=measure_pop.end(); curorg++) delete (*curorg);
//clear the old population
        measure_pop.clear();
        if (generation!=maxgens)
            for (curorg=(pop->organisms).begin(); curorg!=(pop->organisms).end(); ++curorg) {
                measure_pop.push_back(new Organism(*(*curorg),true));
            }
    }
//Create the next generation
    pop->epoch(generation);


    return win;
}

