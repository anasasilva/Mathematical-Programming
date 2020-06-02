#include "Instance.h"

Instance::Instance() :
	graph_name( "" ), n_nodes( 0 ), n_edges( 0 )
{
}

void Instance::import( string file )
{
	ifstream ifs( file.c_str() );
	if( ifs.fail() ) {
		cerr << "could not open input file " << file << "\n";
		exit( -1 );
	}

	graph_name = file.substr( 5, 3 );

	cout << "Reading instance from file " << file << "\n";

	ifs >> n_nodes >> n_edges;
	cout << "Number of nodes: " << n_nodes << "\n";
	cout << "Number of edges: " << n_edges << "\n";

	edges.resize( n_edges );
	incidentEdges.resize( n_nodes );

	u_int id;
	while( ifs >> id ) {
		ifs >> edges[id].v1 >> edges[id].v2 >> edges[id].weight;
		incidentEdges[edges[id].v1].insert( id );
		incidentEdges[edges[id].v2].insert( id );
	}
	ifs.close();
}

void Instance::create( string _graph_name, u_int _n_nodes, u_int _n_edges )
{
	graph_name = _graph_name;
	n_nodes = _n_nodes;
	n_edges = _n_edges;
	edges.resize( n_edges );
	incidentEdges.resize( n_nodes );

	// generate random graph (in a dumb way ...)
	srand( (unsigned) time( NULL ) );
	vector<unordered_set<u_int> > adjacent_nodes( n_nodes );
	// add all root edges
	u_int m = 0;
	for( u_int i = 1; i < n_nodes; i++ ) {
		edges[m].v1 = 0;
		edges[m].v2 = i;
		edges[m].weight = 0;
		incidentEdges[0].insert( m );
		incidentEdges[i].insert( m );
		m++;
	}
	// add further random edges
	while( m < n_edges ) {
		int i = (rand() % (n_nodes - 1)) + 1;
		int j = (rand() % (n_nodes - 2)) + 1;
		if( i <= j ) j++;
		else swap( i, j );

		if( adjacent_nodes[i].count( j ) == 0 ) {
			edges[m].v1 = i;
			edges[m].v2 = j;
			edges[m].weight = rand() % 1000;
			incidentEdges[i].insert( m );
			incidentEdges[j].insert( m );
			adjacent_nodes[i].insert( j );
			m++;
		}
	}

	// write to file
	stringstream file;
	file << "data/" << graph_name << ".dat";
	ofstream ofs( file.str().c_str() );
	if( ofs.fail() ) {
		cerr << "could not open output file " << file.str() << "\n";
		exit( -1 );
	}

	ofs << n_nodes << endl;
	ofs << n_edges << endl;
	for( u_int e = 0; e < n_edges; e++ )
		ofs << e << " " << edges[e].v1 << " " << edges[e].v2 << " " << edges[e].weight << endl;
	ofs.close();
}
