#ifndef __MAIN__CPP__
#define __MAIN__CPP__

#include <iostream>
#include "Tools.h"
#include "Instance.h"
#include "kMST_ILP.h"

using namespace std;

void usage()
{
	cout << "USAGE:\t<program> -f filename -m model [-k <nodes to connect>]\n";
	cout << "EXAMPLE:\t" << "./kmst -f data/g01.dat -m dcc -k 5\n\n";
	exit( 1 );
} // usage

int main( int argc, char *argv[] )
{
	// read parameters
	int opt;
	// default values
	string file( "data/g01.dat" );
	string model_type( "dcc" );
	int k = 5;
	while( (opt = getopt( argc, argv, "f:m:k:" )) != EOF ) {
		switch( opt ) {
			case 'f': // instance file
				file = optarg;
				break;
			case 'm': // algorithm to use
				model_type = optarg;
				break;
			case 'k': // nodes to connect
				k = atoi( optarg );
				break;
			default:
				usage();
				break;
		}
	}

	// read instance
	Instance inst;
	inst.import( file );
//	inst.create( "g10", 2001, 40000 );

	// solve instance
	kMST_ILP ilp( inst, model_type, k );
	ilp.solve();

	return 0;
} // main

#endif // __MAIN__CPP__
