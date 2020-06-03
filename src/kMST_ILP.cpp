#include "kMST_ILP.h"
#include <sstream>
#include <fstream>

kMST_ILP::kMST_ILP( Instance &_instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k ), epsInt( 0.0 ), epsOpt( 0.0 )
{
	if( k == 0 ) k = instance.n_nodes;
}

void kMST_ILP::solve()
{
	try {

		// initialize CPLEX solver
		initCPLEX();

		// add common constraints
		modelCommon();
		// add model-specific constraints
		if( model_type == "scf" ) modelSCF();
		else if( model_type == "mcf" ) modelMCF();
		else if( model_type == "mtz" ) modelMTZ();
		else if( model_type == "dcc" ) cout << "DC-CUT model: no additional constraints\n";
		else if( model_type == "cec" ) cout << "CE-CUT model: no additional constraints\n";
		else {
			cerr << "No existing model chosen\n";
			exit( -1 );
		}

		// build model
		cplex = IloCplex( model );
		cplex.exportModel( "model.lp" );

		// set parameters
		cplex.setParam( IloCplex::Param::Threads, 1 ); // only use a single thread
		cplex.setParam( IloCplex::Param::TimeLimit, 3600 ); // set time limit to 1 hour
		cplex.setParam( IloCplex::Param::WorkMem, 8192 ); // set memory limit to 8 GB

		epsInt = cplex.getParam( IloCplex::Param::MIP::Tolerances::Integrality );
		epsOpt = cplex.getParam( IloCplex::Param::Simplex::Tolerances::Optimality );

		// set cut-callback for cycle-elimination cuts ("cec") or directed connection cuts ("dcc")
		// both for integer (mandatory!) and fractional (optional) solutions
		CutCallback cb( model_type, epsOpt, instance, x, z );
		if( model_type == "dcc" || model_type == "cec" ) {
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate | IloCplex::Callback::Context::Id::Relaxation;
			cplex.use( &cb, contextmask );
		}

		// solve model
		cout << "Calling CPLEX solve ...\n";
		cplex.solve();
		cout << "CPLEX finished.\n\n";
		cout << "CPLEX status: " << cplex.getStatus() << "\n";
		cout << "Branch-and-Bound nodes: " << cplex.getNnodes() << "\n";
		cout << "Objective value: " << cplex.getObjValue() << "\n";
		cout << "CPU time: " << Tools::CPUtime() << "\n\n";

	}
	catch( IloException &e ) {
		cerr << "kMST_ILP: exception " << e.getMessage();
		exit( -1 );
	}
	catch( ... ) {
		cerr << "kMST_ILP: unknown exception.\n";
		exit( -1 );
	}
}

// ----- private methods -----------------------------------------------

void kMST_ILP::initCPLEX()
{
	cout << "initialize CPLEX ... ";
	try {
		env = IloEnv();
		model = IloModel( env );
	}
	catch( IloException &e ) {
		cerr << "kMST_ILP: exception " << e.getMessage();
	}
	catch( ... ) {
		cerr << "kMST_ILP: unknown exception.\n";
	}
	cout << "done.\n";
}

void kMST_ILP::modelCommon()
{
	try {

		// +++++++++++++++++++++++++++++++++++++++++++++++
		// TODO create variables, build common constraints
		// +++++++++++++++++++++++++++++++++++++++++++++++
		
		int nr_nodes = instance.n_nodes;
		int nr_edges = instance.n_edges;

		// Initialization of x - edge or arc selection variables
		// number of x variables = number of arcs (edges in both directions)
		x = IloBoolVarArray(env, 2 * nr_edges);

		for(u_int e = 0; e < 2 * nr_edges; e++)
		{
			stringstream variableX;
			variableX << "x" << e;
			x[e] = IloBoolVar(env, variableX.str().c_str());
		}

		// Objective Function
		IloExpr objective(env);
		for(u_int e = nr_nodes - 1; e < nr_edges; e++)
		{
			// sums the weights in both directions
			// one of them will be 0
			objective += instance.edges.at(e).weight * x[e];
			objective += instance.edges.at(e).weight * x[e + nr_edges];
		}
		model.add(IloMinimize(env, objective));
		objective.end();
	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelCommon: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelCommon: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelSCF()
{
	try {

		// ++++++++++++++++++++++++++++++++++++++++++
		// TODO build single commodity flow model
		// ++++++++++++++++++++++++++++++++++++++++++

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelSCF: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelSCF: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelMCF()
{
	try {

		// ++++++++++++++++++++++++++++++++++++++++++
		// TODO build multi commodity flow model
		// ++++++++++++++++++++++++++++++++++++++++++

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelMCF: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelMCF: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelMTZ()
{
	try {

		int nr_nodes = instance.n_nodes;
		int nr_edges = instance.n_edges;

		// Initialization of u - order node variable
		// Vary from 0 to k because MTZ orders the nodes
		IloIntVarArray u (env, nr_nodes);
		for(u_int v = 0; v < nr_nodes; v++)
		{
			stringstream variableV;
			variableV << "v" << v;
			u[v] = IloIntVar(env, 0, k, variableV.str().c_str());
		}

		// Initialization of z - node selection variables
		// Makes sure only k variables are chosen
		z = IloBoolVarArray(env, nr_nodes);
		for(u_int v = 0; v < nr_nodes; v++)
		{
			stringstream variableZ;
			variableZ << "z" << v;
			z[v] = IloBoolVar(env, variableZ.str().c_str());
		}

		// The ordered nodes in u have to been between 0 and k
		for(u_int v = 1; v < nr_nodes; v++)
		{
			IloNumExpr limitOrder(env);
			limitOrder += u[v];
			model.add(limitOrder <= k);
			limitOrder.end();
		}

		// The route node (the artificial one) is the one with the first id, u[0] = 0
		// And has to have only 1 outgoing edge
		model.add(u[0] == 0);

		IloNumExpr routeOut(env);
		for(u_int e=0; e < nr_nodes - 1; e++)
		{
			routeOut += x[e];
		}
		model.add(routeOut == 1);
		routeOut.end();

		// The number of edges must be k-1 
		IloNumExpr nrEdges(env);
		for(u_int e = nr_nodes - 1; e < nr_edges; e++)
		{
			// both arcs are considered
			// at least one of them will be 0
			nrEdges += x[e];
			nrEdges += x[e + nr_edges];
		}
		model.add(nrEdges == k-1);
		nrEdges.end();

		// // ????????????????
		// IloNumExpr Expr7(env);
		// for(u_int v = 0; v < nr_nodes; v++)
		// {
		// 	Expr7 += u[v];
		// }
		// model.add(Expr7 == (k * (k+1) )/2);
		// Expr7.end();

		// The node z must be 1 if the node is selected
		// The node must have an order in u, so the u value must not be 0
		for(u_int v = 0; v < nr_nodes; v++)
		{
			model.add(u[v] <= k * z[v]);
		}

		// There must be k nodes selected 
		// So there must be k values in z different from 0
		IloNumExpr nrNodesK(env);
		for(u_int v = 0; v < nr_nodes; v++)
		{
			nrNodesK += z[v];
		}
		model.add(nrNodesK == k);
		nrNodesK.end();

		// To iterate the list of incidents edges of each vertex
		unordered_set<u_int>::iterator it;

		// If the edges x are belong to the solution
		// Both of the incident vertexes must belong to it too
		for (u_int v = 1; v < nr_nodes; v++)
		{
			for(it=instance.incidentEdges.at(v).begin(); it != instance.incidentEdges.at(v).end(); it++)
			{
				// the route vertex (artificial one) doesn't belong to the solution
				if(instance.edges.at(*it).v1 != 0) {
					model.add(x[(*it)] <= u[instance.edges.at(*it).v1]);
					model.add(x[(*it) + nr_edges] <= u[instance.edges.at(*it).v1]);
				}

				model.add(x[(*it)] <= u[instance.edges.at(*it).v2]);
				model.add(x[(*it)+ nr_edges] <= u[instance.edges.at(*it).v2]);
			}
		}

		// MTZ
		for(u_int e=0; e < nr_edges; e++)
		{
			IloNumExpr Expr1(env);
			IloNumExpr Expr2(env);

			Expr1 += u[instance.edges.at(e).v1];
			Expr1 += x[e];

			Expr2 += u[instance.edges.at(e).v2];
			Expr2 += k * (1 - x[e]);

			model.add(Expr1 <= Expr2);
			Expr1.end();
			Expr2.end();

			// For the arcs in the order direction
			IloNumExpr Expr3(env);
			IloNumExpr Expr4(env);

			Expr3 += u[instance.edges.at(e).v2];
			Expr3 += x[e + nr_edges];

			Expr4 += u[instance.edges.at(e).v1];
			Expr4 += k * (1 - x[(e)+ nr_edges]);

			model.add(Expr3 <= Expr4 );
			Expr3.end();
			Expr4.end();
		}

		// Make sure that a node is in the solution when the edge connecting it 
		// to another node and the other node are in the solution
		for(u_int e = nr_nodes - 1; e < nr_edges; e++)
		{
			IloNumExpr Expr1(env);
			IloNumExpr Expr2(env);

			Expr1 += z[instance.edges.at(e).v1];
			// at most only 1 will equal 1
			Expr1 += x[e];
			Expr1 += x[e + nr_edges];

			Expr2 += z[instance.edges.at(e).v2] + 1;

			model.add(Expr1 <= Expr2);
			Expr1.end();
			Expr2.end();

			// For the arcs in the other direction
			IloNumExpr Expr3(env);
			IloNumExpr Expr4(env);

			Expr3 += z[instance.edges.at(e).v2];
			Expr3 += x[e];
			Expr3 += x[e + nr_edges];

			Expr4 += z[instance.edges.at(e).v2] + 1;

			model.add(Expr3 <= Expr4);
			Expr3.end();
			Expr4.end();
		}

		// The incoming arc of each selected node must be at most 1
		for(u_int v = 1; v < nr_nodes; v++)
		{
			IloNumExpr nodeFlow(env);

			for(it=instance.incidentEdges.at(v).begin(); it != instance.incidentEdges.at(v).end(); it++)
			{
				if(instance.edges.at(*it).v2 == v)
				{	// incoming arc
					nodeFlow += x[*it];
				}
				else
				{	// outgoing arc
					nodeFlow += x[(*it)+  nr_edges];
				}
			}
			nodeFlow.end();
		}
	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelMTZ: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelMTZ: unknown exception.\n";
		exit( -1 );
	}
}

kMST_ILP::~kMST_ILP()
{
	// free CPLEX resources
//	x.end();
//	z.end();
//	if( model_type == "scf" || model_type == "mcf" ) f.end();
//	else if( model_type == "mtz" ) d.end();
	cplex.end();
	model.end();
	env.end();
}
