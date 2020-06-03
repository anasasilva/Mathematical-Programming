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
		// modelCommon();
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

		char VarName[24];
		for(u_int e=0; e < 2* instance.n_edges; e++){
			sprintf( VarName, "x_%d", e);
			cout << VarName << ": " << cplex.getValue(x[e]) << endl;
		}

		for(u_int e=0; e < instance.n_nodes; e++)
		{
			sprintf( VarName, "u_%d", e);
			cout << VarName << ": " << cplex.getValue(u[e]) << endl;
		}

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
		
		int sq_nodes = instance.n_nodes * instance.n_nodes;
		int nr_nodes = instance.n_nodes;

		// Variables
		x = IloBoolVarArray (env, sq_nodes); // edge or arc selection variables
		f = IloIntVarArray (env, sq_nodes); // flow variables
		d = IloIntVarArray (env, sq_nodes); // node distance variables	
		z = IloBoolVarArray (env, nr_nodes); // node selection variables

		for (u_int i = 0; i < nr_nodes; i++ ) {

			for (u_int j = 0; j < nr_nodes; j++) {
				u_int index = i * nr_nodes + j;

				stringstream xname;
				xname << "x" << i << j;
				x[index] = IloBoolVar(env, xname.str().c_str());

				stringstream fname;
				fname << "f" << i << j;
				f[index] = IloBoolVar(env, 0, 1, fname.str().c_str());

				stringstream dname;
				dname << "d" << i << j;
				d[index] = IloBoolVar(env, dname.str().c_str());
			}

			stringstream zname;
			zname << "z" << i;
			z[i] = IloBoolVar(env, zname.str().c_str());
		}



		// Constraints
		for (u_int i = 0; i < nr_nodes; i++ ) {
			

			for (u_int j = 0; j < nr_nodes; j++) {
				u_int index1 = i * nr_nodes + j;
				u_int index2 = j * nr_nodes + i;

				model.add( f[index1] + f[index2] = x[index1] );



			}
		}
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

		unordered_set<u_int>::iterator it;

		u_int n = instance.n_nodes;
		u_int m = instance.n_edges;

		// (42) // Initialization of x (array of used edges in the solution)
		x = IloBoolVarArray(env, 2*m);
		char VarName[24];
		for(u_int e=0; e < 2*m; e++)
		{
			sprintf( VarName, "x_%d",e);
			x[e] = IloBoolVar(env,VarName);
		}

		// Initialization of u (array of used nodes in the solution)
		u = IloIntVarArray(env, n);
		for(u_int e=0; e < n; e++)
		{
			sprintf( VarName, "u_%d",e);
			u[e] = IloIntVar(env, 0, k,VarName);
		}

		// (43) // Initialization of y (support array to include only k nodes)
		IloBoolVarArray y(env, n);
		for(u_int e=0; e < n; e++)
		{
			sprintf( VarName, "y_%d",e);
			y[e] = IloBoolVar(env,VarName);
		}

		// (29 // objective function
		IloExpr expr(env);
		for(u_int e=n-1; e < m; e++)
		{
			expr += instance.edges.at(e).weight * x[e];
			expr += instance.edges.at(e).weight * x[e+m];
		}
		model.add(IloMinimize(env, expr));
		expr.end();

		// (31) // u[v] has to be between 0 and k
		for(u_int v=1; v < n; v++)
		{
			IloNumExpr Expr2(env);
			Expr2 += u[v];
			model.add(Expr2 <= k);
			Expr2.end();
		}

		// (32) // u[0] must be 0
		IloNumExpr Expr5(env);
		Expr5 += u[0];
		model.add(Expr5 == 0);
		Expr5.end();

		// (33) // exactly 1 outgoing edge from the root
		IloNumExpr Expr4(env);
		for(u_int e=0; e < n-1; e++)
		{
			Expr4 += x[e];
		}
		model.add(Expr4 == 1);
		Expr4.end();

		// (34) // k-1 edges are allowed in the solution
		IloNumExpr Expr6(env);
		for(u_int e=n-1; e < m; e++)
		{
			Expr6 += x[e];
			Expr6 += x[e+m];
		}
		model.add(Expr6 == k-1);
		Expr6.end();

		// (35) // sum of the u value
		IloNumExpr Expr7(env);
		for(u_int v=0; v < n; v++)
		{
			Expr7 += u[v];
		}
		model.add(Expr7 == (k*(k+1))/2);
		Expr7.end();

		// (36) // force y to be 1 if u is in the solution
		for(u_int v=0; v < n; v++)
		{
			IloNumExpr Expr33(env);
			Expr33 += u[v];
			model.add(Expr33 <= k*y[v]);
			Expr33.end();
		}

		// (37) // only k different nodes allowed
		IloNumExpr Expr34(env);
		for(u_int v=0; v < n; v++)
		{
			Expr34 += y[v];
		}
		model.add(Expr34 == k);
		Expr34.end();

		// (39) (40) // if f_ij in the solution -> u[i] and u[j] must be greater than 0
		for(u_int v=1; v < n; v++)
		{
			for(it=instance.incidentEdges.at(v).begin(); it != instance.incidentEdges.at(v).end(); it++)
			{
				if(instance.edges.at(*it).v1 == 0)
				{
					IloNumExpr Expr6(env);
					IloNumExpr Expr7(env);
					Expr6 += x[(*it)];
					Expr7 += x[(*it)+m];
					model.add(Expr6 <= u[instance.edges.at(*it).v2]);
					model.add(Expr7 <= u[instance.edges.at(*it).v2]);
					Expr6.end();
					Expr7.end();
				}
				else
				{
					IloNumExpr Expr6(env);
					IloNumExpr Expr7(env);
					Expr6 += x[(*it)];
					Expr7 += x[(*it)];
					model.add(Expr6 <= u[instance.edges.at(*it).v1]);
					model.add(Expr7 <= u[instance.edges.at(*it).v2]);
					Expr6.end();
					Expr7.end();

					IloNumExpr Expr40(env);
					IloNumExpr Expr41(env);
					Expr40 += x[(*it)+m];
					Expr41 += x[(*it)+m];
					model.add(Expr40 <= u[instance.edges.at(*it).v1]);
					model.add(Expr41 <= u[instance.edges.at(*it).v2]);
					Expr40.end();
					Expr41.end();
				}
			}
		}

		// (30) // miller tucker zemlin formulation
		for(u_int e=0; e < m; e++)
		{
			IloNumExpr Expr4(env);
			Expr4 += u[instance.edges.at(e).v1];
			Expr4 += x[e];

			model.add(Expr4 <= u[instance.edges.at(e).v2] + k*(1-x[e]));
			Expr4.end();

			IloNumExpr Expr12(env);
			Expr12 += u[instance.edges.at(e).v2];
			Expr12 += x[(e)+m];

			model.add(Expr12 <= u[instance.edges.at(e).v1] + k*(1-x[(e)+m]));
			Expr12.end();
		}

		// (41) //
		for(u_int e=n-1; e < m; e++)
		{
			IloNumExpr Expr50(env);
			Expr50 += y[instance.edges.at(e).v1];
			Expr50 += x[e];
			Expr50 += x[e+m];

			model.add(Expr50 <= y[instance.edges.at(e).v2] + 1);
			Expr50.end();

			IloNumExpr Expr51(env);
			Expr51 += y[instance.edges.at(e).v2];
			Expr51 += x[e];
			Expr51 += x[e+m];

			model.add(Expr51 <= y[instance.edges.at(e).v1] + 1);
			Expr51.end();
		}

		// (38) // every node \ 0 must have an incoming flow of <= 1
		for(u_int v=1; v < n; v++)
		{
			IloNumExpr Expr11(env);

			for(it=instance.incidentEdges.at(v).begin(); it != instance.incidentEdges.at(v).end(); it++)
			{
				if(instance.edges.at(*it).v2 == v)
				{	// incoming edge
					Expr11 += x[*it];
				}
				else
				{	// outgoing edge
					Expr11 += x[(*it)+m];
				}
			}
			model.add(Expr11 <= 1);

			Expr11.end();
		}

		// build model
		// cplex = IloCplex( model );
		// // export model to a text file
		// //cplex.exportModel( "model.lp" );
		// // set parameters
		// setCPLEXParameters();

		// // solve model
		// cout << "Calling CPLEX solve ..." << endl;
		// cplex.solve();

	//	for(u_int i=0; i<n; i++)
	//	{
	//		cout << "u[" << i << "] = " << cplex.getValue(u[i]) << endl;
	//		//printf("x_%d;\n",x[i]);
	//	}
	//	for(u_int i=0; i<n; i++)
	//	{
	//		cout << "y[" << i << "] = " << cplex.getValue(y[i]) << endl;
	//		//printf("x_%d;\n",x[i]);
	//	}
	//	for(u_int i=0; i<2*m; i++)
	//	{
	//		cout << "x[" << i << "] = " << cplex.getValue(x[i]) << endl;
	//		//printf("x_%d;\n",x[i]);
	//	}

		// cout << "CPLEX finished." << endl << endl;
		// cout << "CPLEX status: " << cplex.getStatus() << endl;
		// cout << "Branch-and-Bound nodes: " << cplex.getNnodes() << endl;
		// cout << "Objective value: " << cplex.getObjValue() << endl;
		// cout << "CPU time: " << Tools::CPUtime() << endl;
		// cout << "Epsilon: " << cplex.getParam(IloCplex::EpInt) << endl << endl;

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
