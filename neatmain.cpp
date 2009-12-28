#include <tclap/CmdLine.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <unistd.h>	
#include "neat.h"
#include "population.h"
#include "experiments.h"

using namespace std;
using namespace TCLAP;

int main(int argc, char **argv) {

  CmdLine cmd("Maze evolution", ' ', "0.1");
   
  ValueArg<string> maze("m","maze","Maze file",false,"maze.txt","string");
  cmd.add(maze);
  
  ValueArg<string> mcmaze("","mcmaze","MC Maze file",false,"mcmaze.txt","string");
  cmd.add(mcmaze);

  ValueArg<string> genes("z","sg","Starter genes",false,"mazestartgenes","string");
  cmd.add(genes);

  ValueArg<string> settings("s","settings","Settings file",false,"maze.ne","string");
  cmd.add(settings);
 
  ValueArg<string> output("o","output","Output directory",false,"./","string");
  cmd.add(output);

  ValueArg<string> seed_genome("c","seed","Seed Genome",false,"","string");
  cmd.add(seed_genome);

  SwitchArg goal_attract("","goalnotsticky","Goal is not attractor",false);
  cmd.add(goal_attract);

  SwitchArg area_of_interest("","aoi","Enforce pruning of behavior space",false);
  cmd.add(area_of_interest);

  SwitchArg noveltySwitch("n","novelty","Enable novelty search",false);
  cmd.add(noveltySwitch);

  SwitchArg constraintSwitch("","constraint","Enable constraint-based NS",false);
  cmd.add(constraintSwitch);

  SwitchArg generationalSwitch("","gen","Enable generational search",false);
  cmd.add(generationalSwitch);

  SwitchArg mcSwitch("","mc","Enable minimal criteria",false);
  cmd.add(mcSwitch);

  ValueArg<string> nov_measure("","nm","Novelty Measure",false,"std","string");
  cmd.add(nov_measure);

  ValueArg<string> fit_measure("f","fm","Fitness Measure",false,"std","string");
  cmd.add(fit_measure);

  ValueArg<int> extra_param("p","parameter","Extra Parameter",false,0,"int");
  cmd.add(extra_param);

  ValueArg<int> num_samples("a","samples","Num Samples",false,1,"int");
  cmd.add(num_samples);

  ValueArg<int> time_steps("b","timesteps","Num Timesteps",false,400,"int");
  cmd.add(time_steps);

  ValueArg<int> rng_seed("r","random_seed","Random Seed",false,-1,"int");
  cmd.add(rng_seed);

  cmd.parse(argc,argv);

  char mazename[100]="hard_maze.txt";
  char filename[100]="./";
  char settingsname[100]="maze.ne";
  char startgenes[100]="mazestartgenes";
  int param;
  NEAT::Population *p;

  //***********RANDOM SETUP***************//
  /* Seed the random-number generator with current time so that
      the numbers will be different every time we run.    */
  srand( (unsigned)time( NULL )  + getpid());
 
  if(rng_seed.getValue()!=-1)
	srand((unsigned)rng_seed.getValue());

  strcpy(settingsname,settings.getValue().c_str());
  strcpy(mazename,maze.getValue().c_str());
  strcpy(filename,output.getValue().c_str());
  strcpy(startgenes,genes.getValue().c_str());

  NEAT::load_neat_params(settingsname,true);
  param = extra_param.getValue();
  cout<<"loaded"<<endl;

  cout << "Maze: " << mazename << endl;
  cout << "Start genes: " << startgenes << endl;

  set_mcmaze(mcmaze.getValue());
  set_fit_measure(fit_measure.getValue());
  set_nov_measure(nov_measure.getValue());
  set_aoi(area_of_interest.getValue());

  cout << "Timesteps: " << time_steps.getValue() << endl;
  set_timesteps(time_steps.getValue());

  cout << "Num Samples: " << num_samples.getValue() << endl;
  set_samples(num_samples.getValue());

  set_seed(seed_genome.getValue()); 

  cout << "Goal not sticky? "  << goal_attract.getValue() << endl;
  set_goal_attract(!goal_attract.getValue());

  cout << "Minimal criteria engaged? " << mcSwitch.getValue() << endl;
  set_minimal_criteria(mcSwitch.getValue());
  
  set_constraint_switch(constraintSwitch.getValue());
  if(!generationalSwitch.getValue())
{
if(!noveltySwitch.getValue())
      p = maze_novelty_realtime(filename,mazename,param,startgenes,false);
  else
      p = maze_novelty_realtime(filename,mazename,param,startgenes,true);
}
else
{
 p = maze_generational(filename,mazename,param,startgenes,2000,noveltySwitch.getValue());
}
  return(0);
 
}
