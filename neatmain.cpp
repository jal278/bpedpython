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

  ValueArg<string> genes("z","sg","Starter genes",false,"mazestartgenes","string");
  cmd.add(genes);

  ValueArg<string> settings("s","settings","Settings file",false,"maze.ne","string");
  cmd.add(settings);
 
  ValueArg<string> output("o","output","Output directory",false,"./","string");
  cmd.add(output);

  SwitchArg noveltySwitch("n","novelty","Enable novelty search",false);
  cmd.add(noveltySwitch);

  ValueArg<string> nov_measure("g","nm","Novelty Measure",false,"std","string");
  cmd.add(nov_measure);

  ValueArg<string> fit_measure("f","fm","Fitness Measure",false,"std","string");
  cmd.add(fit_measure);

  ValueArg<int> extra_param("p","parameter","Extra Parameter",false,0,"int");
  cmd.add(extra_param);

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
   
if(!noveltySwitch.getValue())
      p = maze_fitness_realtime(filename,mazename,param,fit_measure.getValue(),startgenes);
  else
      p = maze_novelty_realtime(filename,mazename,param,nov_measure.getValue(),startgenes);


  return(0);
 
}
