#include "noveltyset.h"
#include <string.h>
#include "population.h"
#include "organism.h"

//for sorting by novelty
bool cmp(const noveltyitem *a, const noveltyitem* b)
{
return a->novelty < b->novelty;
}

//for sorting by fitness
bool cmp_fit(const noveltyitem *a, const noveltyitem *b)
{
return a->fitness < b->fitness;
}

noveltyitem::noveltyitem(const noveltyitem& item)
{
	added=item.added;
        //TODO: this might cause memory leak in
        //merge_population?
	genotype=new Genome(*(item.genotype));
	phenotype=new Network(*(item.phenotype));
	age=item.age;
	fitness=item.fitness;
	novelty=item.novelty;
	generation=item.generation;
	indiv_number=item.indiv_number;
	for(int i=0;i<(int)item.data.size();i++)
	{
		vector<float> temp;
		for(int j=0;j<(int)item.data[i].size();j++)
			temp.push_back(item.data[i][j]);
		data.push_back(temp);		
	}
}

//merge two populations together according to novelty
Population* noveltyarchive::merge_populations(Population* p1, vector<Organism*> p2)
{

vector<Organism*> total_orgs;
vector<Organism*> merged_orgs;
vector<Organism*>::iterator it;
vector<noveltyitem*>::iterator novit;

//compile the organisms together
for(it = p1->organisms.begin(); it!= p1->organisms.end();it++)
{
	total_orgs.push_back(*it);
	(*it)->blacklist=false;
}

for(it = p2.begin(); it!= p2.end(); it++)
{
	total_orgs.push_back(*it);
	(*it)->blacklist=false;
}

//throw in the archive as well
for(novit = novel_items.begin();novit != novel_items.end(); novit++)
{
	//TODO: just creating these organisms will be a mem leak
	//eventually refactor...
	Organism* arch_org = new Organism(0.1,(*novit)->genotype,0,NULL);
	arch_org->noveltypoint = (*novit);
	total_orgs.push_back(arch_org);
	//or at least delete...?
}

int size = total_orgs.size(); //remove one since we are adding 1st
cout << size << " " << p1->organisms.size() << " " << p2.size() << endl;

//randomly add first member to merged organisms
Organism* last_added = total_orgs[rand()%size];
last_added->blacklist=true;
merged_orgs.push_back(last_added);

//find the closest archive point to each individual
for(it = total_orgs.begin(); it!=total_orgs.end(); it++)
{
	double closest = 100000000.0;
	(*it)->closest = closest;
}

//now greedily add point furthest from archive + merged pop so far
//for(int x=0;x<(size/2)-1;x++)
for(int x=0;x<(p1->organisms.size()-1);x++)
{
	Organism* best=NULL;
	double best_dist= -1000.0;
	for(it = total_orgs.begin(); it!=total_orgs.end(); it++)
	{
		if ((*it)->blacklist)
			continue;

                double new_dist = (*novelty_metric)((*it)->noveltypoint,
					last_added->noveltypoint);

		if (new_dist < (*it)->closest)
			(*it)->closest = new_dist;

		if ((*it)->closest > best_dist)
		{
			best_dist = ((*it)->closest);
			best = *it;
		}
	}
        best->blacklist=true;
	merged_orgs.push_back(best);
}
return new Population(merged_orgs);


}

//evaluate the novelty of the whole population
void noveltyarchive::evaluate_population(Population* p1,vector<Organism*> p2,bool fitness)
{
	vector<Organism*>::iterator it;
	for(it=p1->organisms.begin();it<p1->organisms.end();it++)
		evaluate_individual((*it),p2,fitness);
}

//evaluate the novelty of the whole population
void noveltyarchive::evaluate_population(Population* pop,bool fitness)
{
	Population *p = (Population*)pop;
	vector<Organism*>::iterator it;
	for(it=p->organisms.begin();it<p->organisms.end();it++)
		evaluate_individual((*it),pop->organisms,fitness);
}

//evaluate the novelty of a single individual
void noveltyarchive::evaluate_individual(Organism* ind,vector<Organism*> pop,bool fitness)
{
	float result;
	if(fitness)  //assign fitness according to average novelty
	{
		if(!histogram)
			result = novelty_avg_nn(ind->noveltypoint,-1,false,&pop);
		else
		{		
			result = novelty_histogram(ind->noveltypoint);
		}
		ind->fitness = result;
	} 
	else  //consider adding a point to archive based on dist to nearest neighbor
	{
		if(!histogram)
		{
		result = novelty_avg_nn(ind->noveltypoint,1,false);
		ind->noveltypoint->novelty=result;

                if(!minimal_criteria)
                   ind->noveltypoint->viable=true;

		if(ind->noveltypoint->viable && add_to_novelty_archive(result))
				add_novel_item(ind->noveltypoint);
		}
	}
}
