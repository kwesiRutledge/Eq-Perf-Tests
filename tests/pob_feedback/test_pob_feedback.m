function tests = test_pob_feedback
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

function [BG_out,controller_out,P_target,lcsas,P_x0,P_u] = get_reach_contr1()
	%get_reach_contr1
	%Description:
	%	This function synthesizes a controller that steers the simple one dimensional system towards the target 

	%% Constant %%

	[lcsas,~,P_u,P_x0] = get_simple_1dsys();
	
	P_target = Polyhedron('A',-1,'b',-2);

	%% Algorithms %%

	[BG_out,controller_out] = lcsas.synth_robust_reach_contr( P_x0 , P_u , 'P_target', P_target , 'debug' , 0 );

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
	Pw2 = 1.5*eta_w + Pw1;

	eta_u = eta_w; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w;
	f2 = -eta_w;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1) ];

	sys_out = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

function [BG_out,controller_out,P_target] = get_reach_contr2( lcsas_in, language_in, P_u, P_x0 )
	%get_reach_contr2
	%Description:
	%	This function synthesizes a controller that steers the simple one dimensional system towards the target
	%	the belief graph that it uses DOES NOT use projection.

	%% Constant %%
	
	P_target = Polyhedron('A',-1,'b',-2);

	%% Algorithms %%

	[BG_out,controller_out] = lcsas_in.synth_robust_reach_contr( P_x0 , P_u , 'P_target', P_target , 'debug' , 0 , 'UseProjection' , false );

function [BG_out] = get_unconnected_bg1()
	%get_unconnected_bg1
	%Description:
	%	Creates a very particular BeliefGraph to test the BeliefLanguage construction.

	%% Constants %%
	T = 4;
	[lcsas,~,P_u,P_x0] = get_simple_1dsys();

	%% Create Graph %%
	bn1 = BeliefNode(Language([1,2]),0);
	bn3 = BeliefNode(Language([1,2]),1);
	bn5 = BeliefNode(Language([1,2]),2);
	bn7 = BeliefNode(Language([1,2]),3);

	bn2 = BeliefNode(Language([2]),0);
	bn4 = BeliefNode(Language([2]),1);
	bn6 = BeliefNode(Language([2]),2);
	bn8 = BeliefNode(Language([2]),3);

	edges = [ 	1,3;
				3,5;
				5,7;
				2,4;
				4,6;
				6,8 ];

	BG_out = BeliefGraph(lcsas,P_u,P_x0,'return_empty',true);
	BG_out.N = [bn1,bn2,bn3,bn4,bn5,bn6,bn7,bn8];
	BG_out.E = edges;

function test1_belief_lang(testCase)
	%test1_belief_lang
	%Description:
	%	Verifies that the belief language of get_unconnected_bg1 is as I think it is
	%	2 words ([1,3,5,7].[2,4,6,8])

	%% Constants %%

	BG1 = get_unconnected_bg1();

	%% Algorithm %%

	BG1.BeliefLanguage = BG1.get_belief_language();

	assert( (BG1.BeliefLanguage.cardinality() == 2) & BG1.BeliefLanguage.contains([1,3,5,7]) & BG1.BeliefLanguage.contains([2,4,6,8]) )

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
	