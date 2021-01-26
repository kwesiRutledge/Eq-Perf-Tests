function tests = behavior_set2timeTest
	% disp(localfunctions)
	tests = functiontests(localfunctions);

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	%% Algorithm

	%Include Yalmip
	%include_fcns2({'mosek','gurobi','MPT3'})

	%Add Local Functions to Path
	addpath(genpath('../functions'));

    
function test_behavior_set2time1(testCase)
	%test_behavior_set2time1
	%Description:
	%	This test is meant to evaluate the function LCSAS.behavior_set2timeTest for:
	%		- internal behavior sets
	%		- constructed with method 1
	%		- output feedback

    % Include libraries
    include_relevant_libraries();
    
    % Algorithm
    
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
	fb_method = 'output';
	debug_flag = 0;

	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , false , 'fb_method' , fb_method , 'ConsistencySetVersion' , 1 );

	%Create initial node set
	
	t0 = 2;
	[ ~ , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												t0, L1, P_u,P_x0, ...
												'fb_method',fb_method,'debug_flag',debug_flag, ...
												'use_proj',empty_bg.UsedProjection );

	%extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );
	check_ib = false(1,length(initial_internal_behavior_sets));
	for ib_idx = 1:length(initial_internal_behavior_sets)
		check_ib(ib_idx) = ( t0 == lcsas0.behavior_set2time( initial_internal_behavior_sets(ib_idx) , 'InternalBehaviorSet_1word' , ...
													'ConsistencySetVersion' , empty_bg.ConsistencySetVersion , ...
													'FeedbackMethod' , empty_bg.FeedbackMethod ) );
	end

	assert( all(check_ib) )

function test_behavior_set2time2(testCase)
	%test_behavior_set2time2
	%Description:
	%	This test is meant to evaluate the function LCSAS.behavior_set2timeTest for:
	%		- internal behavior sets
	%		- constructed with method 1
	%		- state feedback

    % Include libraries
    include_relevant_libraries();
    
    % Algorithm
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
	fb_method = 'state';
	debug_flag = 0;

	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , false , 'fb_method' , fb_method , 'ConsistencySetVersion' , 1 );

	%Create initial node set
	
	t0 = 2;
	[ ~ , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												t0, L1, P_u,P_x0, ...
												'fb_method',fb_method,'debug_flag',debug_flag, ...
												'use_proj',empty_bg.UsedProjection );

	%extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );
	check_ib = false(1,length(initial_internal_behavior_sets));
	for ib_idx = 1:length(initial_internal_behavior_sets)
		check_ib(ib_idx) = ( t0 == lcsas0.behavior_set2time( initial_internal_behavior_sets(ib_idx) , 'InternalBehaviorSet_1word' , ...
													'ConsistencySetVersion' , empty_bg.ConsistencySetVersion , ...
													'FeedbackMethod' , empty_bg.FeedbackMethod ) );
	end

	assert( all(check_ib) )

function test_behavior_set2time3(testCase)
	%test_behavior_set2time3
	%Description:
	%	This test is meant to evaluate the function LCSAS.behavior_set2timeTest for:
	%		- internal behavior sets
	%		- constructed with method 2
	%		- state feedback

    % Include libraries
    include_relevant_libraries();
    
    % Algorithm
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
	fb_method = 'output';
	debug_flag = 0;
	cs_version = 2;

	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , false , 'fb_method' , fb_method , 'ConsistencySetVersion' , cs_version );

	%Create initial node set
	
	t0 = 2;
	[ ~ , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												t0, L1, P_u,P_x0, ...
												'fb_method',fb_method,'debug_flag',debug_flag, ...
												'use_proj',empty_bg.UsedProjection, ...
												'ConsistencySetVersion',cs_version );

	%extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );
	check_ib = false(1,length(initial_internal_behavior_sets));
	for ib_idx = 1:length(initial_internal_behavior_sets)
		check_ib(ib_idx) = ( t0 == lcsas0.behavior_set2time( initial_internal_behavior_sets(ib_idx) , 'InternalBehaviorSet_1word' , ...
													'ConsistencySetVersion' , empty_bg.ConsistencySetVersion , ...
													'FeedbackMethod' , empty_bg.FeedbackMethod ) );
	end

	assert( all(check_ib) )