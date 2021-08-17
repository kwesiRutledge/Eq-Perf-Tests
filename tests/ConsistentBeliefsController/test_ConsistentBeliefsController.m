function tests = test_ConsistentBeliefsController
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function [sys_out,L1,P_u,P_x0] = get_simple_1dsys()
	%Description:

	L1 = Language([1,2,1,2],[1,1,1,1],[3,1,1,1]);
	T = length(L1.words{1});

	A1 = 0.9;
	B1 = 1;
	C1 = 1; %[1,0];
	
	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	eta_v = 0.3;
	Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

	eta_u = 1.0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w;
	f2 = -eta_w;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

	sys_out = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

function test1_constructor(testCase)
	%Description:
	%	Tests the ability to construct the standard ConsistentBeliefController object.

	%% Constant %%

	[system1,L1,U,X0] = get_simple_1dsys(); 
	TempSequences = [ repmat(L1,3,1) , [ L1 ; L1 ; Language(L1.words{1}) ] ];

	K_set = { 1 , 2 };
	k_set = { 1 , 1 };

	%% Algorithm %%

	cbc1 = ConsistentBeliefsController( system1 , TempSequences , K_set , k_set );

function [sys_out,L1,P_u,P_x0] = get_simpler_1dsys_v2()
	%Description:

	L1 = Language([1,1],[2,2]);
	T = length(L1.words{1});

	A1 = 1;
	B1 = 1;
	C1 = 1; %[1,0];
	
	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	% Pw2 = -2*eta_w + Pw1;

	eta_u = eta_w; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w;
	f2 = -eta_w;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1) ];

	sys_out = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 , 'U' , P_u );

function test1_compute_control(testCase)
	%test1_apply_control
	%Description:
	%	Verifies that the controller will output the correct values
	%	when the history is defined appropriately.

	%% Constants %%

	[system1,L1,P_u,P_x0] = get_simpler_1dsys_v2();
	KnowledgeSequences = [ [ L1 ; Language(L1.words{1}) ] , [ L1 ; Language(L1.words{2}) ] ];

	[ n_x , n_u , n_y , n_w , n_v ] = system1.Dimensions();
	T = length(L1.words{1});
	K = { zeros(T*n_x,T*n_w) , zeros(T*n_x,T*n_w) };
	k = { [1;3] , [2;4] };

	%% Algorithm %%

	cbc1 = ConsistentBeliefsController( system1 , KnowledgeSequences , K , k );
	cbc1.x_hist = [0.1;0.4];
	cbc1.u_hist = 0;

	assert( cbc1.compute_control() == 3 )

function test1_simulation(testCase)
	%test1_external_behavior_explained_by
	%Description:
	%

	[BG1,contr1,P_target,lcsas,P_x0,P_u] = get_reach_contr1(); 

	[x, u , y , ~ ] = contr1.simulate_1run(lcsas, P_x0 );

	assert(P_target.contains(x(end)) )

function test2_simulation(testCase)
	%test1_external_behavior_explained_by
	%Description:
	%

	%% Constants %%

	[sys_out,L1,P_u,P_x0] = get_simple_1dsys();

	[BG2,contr2,P_target] = get_reach_contr2( sys_out , L1 , P_u , P_x0 );

	if ~isa(contr2,'POB_Feedback')
		error(['contr2 is not of class POB_Feedback. Instead it is of class ' class(contr2)] )
	end

	[x, u , y , ~ ] = contr2.simulate_1run(sys_out, P_x0 );

	assert(P_target.contains(x(end)) )

function test1_pob_feedback(testCase)
	%Description:
	%	Catch a bad POB_Feedback object constructor
	%	

	%% Constants
	[sys1,L1,P_u,P_x0] = get_simple_1dsys();

	%% Algorithm

	try 
		POB_Feedback(1,sys1,1)
		assert(false)
	catch e
		% disp(e.message)
		assert(strcmp(e.message,['Unexpected type for first input: double. Choose BeliefGraph or LCSAS object as first input.']))
	end

function test2_pob_feedback(testCase)
	%Description:
	%	Create using a simple lcsas example.

	%% Constants
	[sys1,L1,P_u,P_x0] = get_simple_1dsys();

	%% Algorithm

	pob1 = POB_Feedback(sys1,eye(2),ones(2,1))

	assert(all(all(pob1.F_set == eye(2))) && all(pob1.u0_set == ones(2,1)))

function test1_GetReachableSetAt(testCase)
	%Description:
	%	Testing the error catching capabilities of the function.

	%% Constants
	[sys1,L1,P_u,P_x0] = get_simple_1dsys();
	bad_mode_val = 4;

	%% Algorithm

	pob1 = POB_Feedback(sys1,eye(2),ones(2,1));

	try
		pob1.GetReachableSetAt(2,bad_mode_val)
	catch e
		disp(e.message)
		assert(strcmp(e.message,['There is not a word ' num2str(bad_mode_val) ' in the language L. It has only ' num2str(sys1.L.cardinality()) ' words.'  ]))
	end

function test2_GetReachableSetAt(testCase)
	%Description:
	%	Testing the correct behavior of the function for a simple static system.

	%% Constants
	n_x = 2; n_y = 2;

	A1 = eye(2); B1 = [0;1]; C1 = eye(2);
	eta_v = 0; eta_w = 0;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));

	eta_x0 = 0.3;
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	d1 = Aff_Dyn(A1,B1,zeros(2,1),C1,Pw1,Pv1);

	sys1 = LCSAS( [d1] , Language([1,1,1,1]) , 'X0' , P_x0 );


	%% Algorithm

	pob1 = POB_Feedback(sys1,{zeros(2)},{zeros(2,1)});

	reach1 = pob1.GetReachableSetAt(1,1)

	assert( (reach1 <= P_x0) && (P_x0 <= reach1) )
	