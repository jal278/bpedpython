#include "modularity.hpp"

using namespace boost;

int main(int argc, char **argv)
{
    std::cout<<"created"<<std::endl;
typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
typedef property_map<Graph, vertex_index_t>::type IndexMap;
typedef std::pair<int,int> Edge;
Edge edge_array[] = {Edge(0,2),Edge(2,0),Edge(1,3),Edge(3,1)};
Graph g(4);
IndexMap index= get(vertex_index,g);
for(int i=0;i<4;i++)
 add_edge(i,i,g); //edge_array[i].first,edge_array[i].second,g);


    std::cout<<"loaded"<<std::endl;

    std::vector<std::set<Graph::vertex_descriptor> > mods;
    float mod1 = mod::h_modules(g, mods);
    std::cout<<"modularity:"<<mod1<<std::endl;
    BOOST_FOREACH(std::set<Graph::vertex_descriptor> s, mods)
      {
	std::cout<<"---"<<std::endl;
	BOOST_FOREACH(Graph::vertex_descriptor v, s)
	  std::cout<<index[v]<<std::endl;
      }
}

