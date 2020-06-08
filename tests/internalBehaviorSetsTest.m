function tests = internalBehaviorSetsTest
	disp(localfunctions)
	tests = functiontests(localfunctions);

function testInternalBehaviorSets1(testCase)
	%testInternalBehaviorSets1.m
	%Description:
	%	This test script is meant to evaluat the function get_all_consistent_behavior_sets.m which is a member function of
	%	the belief graph object.

	L1 = Language([1,2,1,2],[2,1,2,1],[1,1,1,1]);
	T = 4;

	%Create a simple Language Constrainted Switching System
	A1 = eye(2);
	B1 = [0;1];
	C1 = [1,0];

	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = zeros(2,1);

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	%Create BeliefGraph

	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , false );

	%Create initial node set
	debug_flag = 0;
	N0 = empty_bg.get_initial_beliefnodes('verbosity',debug_flag);

	n0 = N0(1);

	fb_method = 'output';
	[ ~ , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												n0.t+1, n0.subL, P_u,P_x0, ...
												'fb_method',fb_method,'debug_flag',debug_flag, ...
												'use_proj',empty_bg.UsedProjection );

	extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );

	[ ~ , empty_set_flags ] = empty_bg.find_empty_observation_polyhedra( extended_internal_behavior_sets );

	assert( ~any( empty_set_flags ) )

function testInternalBehaviorSets2(testCase)
	%testInternalBehaviorSets2.m
	%Description:
	%	This test is meant to evaluate the function get_all_consistent_behavior_sets.m which is a member function of
	%	the belief graph object. In this case, we choose a language that leads to some non-overlapping external behaviors.

	L1 = Language([1,1,1,1],[2,2,2,2],[3,3,3,3]);
	T = 4;

	%Create a simple Language Constrainted Switching System
	A1 = eye(2);
	B1 = [0;1];
	C1 = [1,0];

	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.0;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));

	eta_u = 0.0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = [0;1];
	f2 = [0;-1];

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	%Create BeliefGraph

	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , false );

	%Create initial node set
	debug_flag = 0;
	N0 = empty_bg.get_initial_beliefnodes('verbosity',debug_flag);

	n0 = N0(1);

	fb_method = 'output';
	[ ~ , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												n0.t+1, n0.subL, P_u,P_x0, ...
												'fb_method',fb_method,'debug_flag',debug_flag, ...
												'use_proj',empty_bg.UsedProjection );

	extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );

	[ ~ , empty_set_flags ] = empty_bg.find_empty_observation_polyhedra( extended_internal_behavior_sets );

	assert( all( empty_set_flags == false(length(empty_set_flags),1) ) )