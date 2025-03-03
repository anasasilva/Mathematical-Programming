#include "kMST_ILP.h"
#include <sstream>
#include <fstream>

kMST_ILP::kMST_ILP( Instance &_instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k ), epsInt( 0.0 ), epsOpt( 0.0 )
{
	if( k == 0 ) k = instance.n_nodes; // excludes the artificial od
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
		cplex.setParam( IloCplex::Param::TimeLimit, 1800 ); // set time limit to half an hour
		cplex.setParam( IloCplex::Param::WorkMem, 4096 ); // set memory limit to 4 GB

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
		
		int nr_nodes = instance.n_nodes;
		int nr_edges = instance.n_edges;
		stringstream variableName;

		// Initialization of x - edge or arc selection variables
		// number of x variables = number of arcs (edges in both directions)
		x = IloBoolVarArray(env, 2 * nr_edges);

		for(u_int e = 0; e < 2 * nr_edges; e++)
		{
			variableName << "x" << e;
			x[e] = IloBoolVar(env, variableName.str().c_str());
		}

		// Initialization of z - node selection variables
		// Makes sure only k nodes are chosen
		z = IloBoolVarArray(env, nr_nodes);
		for(u_int v = 0; v < nr_nodes; v++)
		{
			variableName << "z" << v;
			z[v] = IloBoolVar(env, variableName.str().c_str());
		}

		// The route node (the artificial one) has to have only 1 outgoing edge
		IloNumExpr routeOut(env);
		for(u_int e = 0; e < nr_nodes - 1; e++)
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


		// Make sure that a node is in the solution when the edge connecting it 
		// to another node and the other node are in the solution
		for(u_int e = nr_nodes - 1; e < nr_edges; e++)
		{
			IloNumExpr Expr1(env);
			IloNumExpr Expr2(env);

			// at most only 1 will equal 1
			Expr1 += x[e];
			Expr1 += x[e + nr_edges];

			Expr2 += z[instance.edges.at(e).v2];

			model.add(Expr1 <= Expr2);
			Expr1.end();
			Expr2.end();

			// For the arcs in the other direction
			IloNumExpr Expr3(env);
			IloNumExpr Expr4(env);

			Expr3 += x[e];
			Expr3 += x[e + nr_edges];

			Expr4 += z[instance.edges.at(e).v2];

			model.add(Expr3 <= Expr4);
			Expr3.end();
			Expr4.end();
		}

		// There must be k nodes selected 
		// So there must be k + 1 values in z that are 1 (includes artificial node)
		IloNumExpr nrNodesK(env);
		for (u_int v = 0; v < nr_nodes; v++)
		{
			nrNodesK += z[v];
		}
		model.add(nrNodesK == k + 1);
		nrNodesK.end();


		// Make sure the nodes are in the solution when the edge that connects them is in the solution
		for (u_int e = 0; e < 2 * nr_edges; e++)
		{
			model.add(x[e] <= z[instance.edges.at(e % nr_edges).v1]);
			model.add(x[e] <= z[instance.edges.at(e % nr_edges).v2]);
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

		int nr_nodes = instance.n_nodes;
		int nr_edges = instance.n_edges;
		stringstream varName;
		unordered_set<u_int>::iterator it;

		// Initialization of flow variable
		// Nr of flow varibles = nr of arcs = 2 * nr of edges
		IloIntVarArray f(env, 2 * nr_edges);
		for (u_int e = 0; e < 2 * nr_edges; e++)
		{
			varName << "x" << e;
			f[e] = IloIntVar(env, 0, k, varName.str().c_str());
		}

		// Flow expression of the artificial node
		// Node 0 must have k commodities going out
		IloNumExpr flow0(env);
		for (it = instance.incidentEdges.at(0).begin(); it != instance.incidentEdges.at(0).end(); it++)
		{
			// outgoing edge from 0
			if (instance.edges.at(*it).v1 == 0)
			{	// considers both arcs
				flow0 += f[*it];
				flow0 -= f[(*it) + nr_edges];
			}
			// incoming edge from 0
			else
			{
				flow0 -= f[*it];
				flow0 += f[(*it) + nr_edges];
			}

		}
		model.add(flow0 == k);
		flow0.end();


		// Flow conservation constraint: the difference between what comes in and out of a node is 1
		// or 0 (if the initial flow is 0 which means the node doesn't belong to the solution)
		for(u_int v = 1; v < nr_nodes; v++)
		{
			IloNumExpr flowConstraint(env);

			for(it=instance.incidentEdges.at(v).begin(); it != instance.incidentEdges.at(v).end(); it++)
			{
				if(instance.edges.at(*it).v1 == v)
				{	// outgoing edge
					flowConstraint -= f[*it];
					flowConstraint += f[(*it) + nr_edges];
				}
				else
				{	// incoming edge
					flowConstraint += f[*it];
					flowConstraint -= f[(*it) + nr_edges];
				}
			}
			model.add(flowConstraint == z[v]);

			flowConstraint.end();
		}

		// Make sure that if an edge as flow, then the edge belongs to the solution
		for (u_int e = 0; e < 2 * nr_edges; e++)
		{
			model.add(f[e] <= (k * x[e]));
		}
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

		int nr_nodes = instance.n_nodes;
		int nr_edges = instance.n_edges;
		stringstream varName;
		unordered_set<u_int>::iterator it;

		// Initialization of f ... size: |V| * |E|
		IloNumVarArray f (env, 2 * nr_edges * (nr_nodes - 1));
		for (u_int w = 0; w < nr_nodes - 1; w++)
		{
			for (u_int e = 0; e < 2 * nr_edges; e++)
			{
				varName << "f" << e << "_" << w;
				f[(2 * nr_edges * w) + e] = IloNumVar(env, 0, 1, varName.str().c_str());
			}
		}

		// For each commodity, we send out exactly 1 unit from the route/artificial node
		for (u_int w = 0; w < nr_nodes - 1; w++)
		{
			IloNumExpr routeFlow(env);
			for (it = instance.incidentEdges.at(0).begin(); it != instance.incidentEdges.at(0).end(); it++)
			{
				if(instance.edges.at(*it).v1 == 0)
				{	// outgoing edge
					routeFlow += f[(2 * w * nr_edges ) + (*it)];
					routeFlow -= f[(2 * w * nr_edges) + ((*it) + nr_edges)];
				}
			}
			model.add(routeFlow <= 1);
			routeFlow.end();
		}

		// In total the route node must send exactly k commodities
		IloNumExpr totalCommoditiesRoute(env);
		for (u_int w = 0; w < nr_nodes - 1; w++)
		{
			for (it = instance.incidentEdges.at(0).begin(); it != instance.incidentEdges.at(0).end(); it++)
			{
				if(instance.edges.at(*it).v1 == 0)
				{	// outgoing edge
					totalCommoditiesRoute += f[(2 * nr_edges * w) + (*it)];
					totalCommoditiesRoute -= f[(2 * nr_edges * w) + ((*it) + nr_edges)];
				}
			}
		}
		model.add(totalCommoditiesRoute == k);
		totalCommoditiesRoute.end();

		// A node can receive at most one commodity
		for (u_int w = 0; w < nr_nodes - 1; w++)
		{
			IloNumExpr sendCommodities(env);
			for (it = instance.incidentEdges.at(w + 1).begin(); it != instance.incidentEdges.at(w + 1).end(); it++)
			{
				if (instance.edges.at(*it).v2 == w + 1)
				{	// incoming edge normal
					sendCommodities += f[(2 * nr_edges * w) + (*it)];
				}
				else
				{	// incoming artificial edge
					sendCommodities += f[(2 * nr_edges * w) + ((*it) + nr_edges)];
				}
			}
			model.add(sendCommodities <= 1);
			sendCommodities.end();
		}

		// A node must receive in total k commodities
		IloNumExpr totalCommodities(env);
		for (u_int w = 0; w < nr_nodes - 1; w++)
		{
			for (it = instance.incidentEdges.at(w + 1).begin(); it != instance.incidentEdges.at(w + 1).end(); it++)
			{
				if (instance.edges.at(*it).v2 == w + 1)
				{	// incoming edge normal
					totalCommodities += f[(2 * nr_edges * w) + (*it)];
				}
				else
				{	// incoming artificial edge
					totalCommodities += f[(2 * nr_edges * w) + ((*it) + nr_edges)];
				}
			}
		}
		model.add(totalCommodities == k);
		totalCommodities.end();

		// If the node is not the target node of the comodity, the node can't consume it
		// It has to push it forward
		for (u_int w = 0; w < nr_nodes - 1; w++)
		{
			for (u_int j = 1; j < nr_nodes; j++)
			{
				if(j == w + 1)
					continue;

				IloNumExpr Expr1(env);

				for (it = instance.incidentEdges.at(j).begin(); it != instance.incidentEdges.at(j).end(); it++)
				{
					if (instance.edges.at(*it).v2 == j)
					{

						Expr1 += f[(2 * nr_edges * w) + (*it)];
						Expr1 -= f[(2 * nr_edges * w) + (*it + nr_edges)];

					}
					else if(instance.edges.at(*it).v1 == j)
					{

						Expr1 -= f[(2 * nr_edges * w) + (*it)];
						Expr1 += f[(2 * nr_edges * w) + (*it + nr_edges)];

					}
				}

				model.add(Expr1 == 0);
				Expr1.end();
			}
		}

		// Make sure that if an edge has flow, then the edge belongs to the solution
		for (u_int w = 0; w < nr_nodes - 1; w++)
		{
			for (u_int e = 0; e < 2 * nr_edges; e++)
			{
				model.add(f[(2 * w * nr_edges) + e] <= (k * x[e]));
			}
		}
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
		stringstream variableName;
		unordered_set<u_int>::iterator it;

		// Initialization of u - order node variable
		// Vary from 0 to k because MTZ orders the nodes
		IloIntVarArray u (env, nr_nodes);
		for(u_int v = 0; v < nr_nodes; v++)
		{
			variableName << "v" << v;
			u[v] = IloIntVar(env, 0, k, variableName.str().c_str());
		}

		// The ordered nodes in u have to been between 0 and k
		for(u_int v = 1; v < nr_nodes; v++)
		{
			IloNumExpr limitOrder(env);
			limitOrder += u[v];
			model.add(limitOrder <= k);
			limitOrder.end();
		}

		// The node z must be 1 if the node is selected
		// The node must have an order in u, so the u value must not be 0
		for(u_int v = 0; v < nr_nodes; v++)
		{
			model.add(u[v] <= k * z[v]);
		}

		// The route node (the artificial one) is the one with the first id, u[0] = 0
		model.add(u[0] == 0);

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

		// MTZ constraint
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

		// The incoming arc of each selected node must be at most 1
		for(u_int v = 1; v < nr_nodes; v++)
		{
			IloNumExpr nrArcs(env);

			for(it=instance.incidentEdges.at(v).begin(); it != instance.incidentEdges.at(v).end(); it++)
			{
				if(instance.edges.at(*it).v2 == v)
				{	// incoming edge
					nrArcs += x[*it];
				}
				else
				{	// outgoing edge
					nrArcs += x[(*it)+  nr_edges];
				}
			}
			model.add(nrArcs <= 1);
			nrArcs.end();
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
