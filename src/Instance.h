#ifndef __INSTANCE__H__
#define __INSTANCE__H__

#include "Tools.h"
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <unordered_set>

using namespace std;

class Instance
{

public:

	struct Edge
	{
		u_int v1, v2;
		int weight;
	};

	// input graph
	string graph_name;

	// number of nodes and edges
	u_int n_nodes, n_edges;

	// array of edges
	vector<Edge> edges;

	// incident edge indices for each node
	vector<unordered_set<u_int> > incidentEdges;

	// constructor
	Instance();

	// import from file
	void import( string file );

	// create randomly and write it to file
	void create( string out_file, u_int _n_nodes, u_int _n_edges );

};

#endif //__INSTANCE__H__
