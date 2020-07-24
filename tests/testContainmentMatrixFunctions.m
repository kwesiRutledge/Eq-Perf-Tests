function tests = testContainmentMatrixFunctions
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function testContainmentMatrixFunctions1(testCase)
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

	debug_flag = 0;
	N0 = empty_bg.get_initial_beliefnodes('verbosity',debug_flag);
	n0 = N0(1);

	[ initial_external_behavior_sets , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												n0.t+1, n0.subL, P_u,P_x0, ...
												'fb_method','output','debug_flag',debug_flag, ...
												'use_proj',true );

	%Extend the number of consistency sets by performing intersections.
	subL = n0.subL;
	[subL_powerset, subL_index_powerset] = subL.powerset();
	external_behavior_sets = initial_external_behavior_sets;
	for powerset_idx = (subL.cardinality()+1):length(subL_powerset)
		powerset_idcs_elt = subL_index_powerset{powerset_idx};
		external_behavior_sets(powerset_idx) = initial_external_behavior_sets(powerset_idcs_elt(1));
		
		for elt_idx = 2:length(powerset_idcs_elt)
			external_behavior_sets(powerset_idx) = external_behavior_sets(powerset_idx).intersect( initial_external_behavior_sets(powerset_idcs_elt(elt_idx)) );
		end

		%Update list tracking emptiness
		% consistency_set_is_empty(powerset_idx) = consistency_sets(powerset_idx).isEmptySet;
	end

	extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );

	[ ~ , empty_set_flags_eb ] = empty_bg.find_empty_observation_polyhedra( external_behavior_sets );
	[ ~ , empty_set_flags_ib ] = empty_bg.find_empty_observation_polyhedra( extended_internal_behavior_sets );

	assert( all( empty_set_flags_ib == empty_set_flags_eb ) )

function testContainmentMatrixFunctions2(testCase)
	%testInternalBehaviorSets2.m
	%Description:
	%	This test script is meant to evaluate the function get_all_consistent_behavior_sets.m which is a member function of
	%	the belief graph object.
	%	In this example, there are different dynamics for each mode, but because one state is unobservable the external behavior sets should
	%	still all look identical.

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

	f1 = [0;1];

	aff_dyn_list = [	Aff_Dyn(A1,B1, f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,-f1,C1,Pw1,Pv1) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	%Create BeliefGraph
	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , true );

	debug_flag = 0;
	N0 = empty_bg.get_initial_beliefnodes('verbosity',debug_flag);
	n0 = N0(1);

	[ initial_external_behavior_sets , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												n0.t+1, n0.subL, P_u, P_x0, ...
												'fb_method','output', 'debug_flag',debug_flag, ...
												'use_proj', true , 'reduce_representation' , true );

	% all(all(initial_external_behavior_sets(1).A == initial_external_behavior_sets(2).A))
	% [initial_external_behavior_sets(1).A,initial_external_behavior_sets(1).b]
	% [initial_internal_behavior_sets(1).A,initial_internal_behavior_sets(1).b]
	% [initial_internal_behavior_sets(1).Ae,initial_internal_behavior_sets(1).be]
	% all(all(initial_external_behavior_sets(1).b == initial_external_behavior_sets(2).b))
	% all(all(initial_external_behavior_sets(3).A == initial_external_behavior_sets(2).A))
	% all(all(initial_external_behavior_sets(3).b == initial_external_behavior_sets(2).b))

	%Extend the number of consistency sets by performing intersections.
	subL = n0.subL;
	[subL_powerset, subL_index_powerset] = subL.powerset();
	external_behavior_sets = initial_external_behavior_sets;
	for powerset_idx = (subL.cardinality()+1):length(subL_powerset)
		powerset_idcs_elt = subL_index_powerset{powerset_idx};
		external_behavior_sets(powerset_idx) = initial_external_behavior_sets(powerset_idcs_elt(1));
		
		for elt_idx = 2:length(powerset_idcs_elt)
			external_behavior_sets(powerset_idx) = external_behavior_sets(powerset_idx).intersect( initial_external_behavior_sets(powerset_idcs_elt(elt_idx)) );
		end
	end

	extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );

	[ ~ , empty_set_flags_eb ] = empty_bg.find_empty_observation_polyhedra( external_behavior_sets );
	[ ~ , empty_set_flags_ib ] = empty_bg.find_empty_observation_polyhedra( extended_internal_behavior_sets );

	assert( all( empty_set_flags_ib == false(length(subL_index_powerset),1) ) )

function testContainmentMatrixFunctions3(testCase)
	%testInternalBehaviorSets3.m
	%Description:
	%	This test script is meant to evaluate the function get_all_consistent_behavior_sets.m which is a member function of
	%	the belief graph object.
	%	In this example, there are different dynamics for each mode, and two of the modes are identical but one should be noticeably different.

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

	f1 = [1;0];

	aff_dyn_list = [	Aff_Dyn(A1,B1, f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,-f1,C1,Pw1,Pv1) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	%Create BeliefGraph
	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , true );

	N0 = empty_bg.get_initial_beliefnodes();
	n0 = N0(1);

	debug_flag = 0;
	[ initial_external_behavior_sets , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												n0.t+1, n0.subL, P_u, P_x0, ...
												'fb_method','output', 'debug_flag',debug_flag, ...
												'use_proj', true , 'reduce_representation' , true );

	% all(all(initial_external_behavior_sets(1).A == initial_external_behavior_sets(2).A))
	% [initial_external_behavior_sets(1).A,initial_external_behavior_sets(1).b]
	% [initial_internal_behavior_sets(1).A,initial_internal_behavior_sets(1).b]
	% [initial_internal_behavior_sets(1).Ae,initial_internal_behavior_sets(1).be]
	% all(all(initial_external_behavior_sets(1).b == initial_external_behavior_sets(2).b))
	% all(all(initial_external_behavior_sets(3).A == initial_external_behavior_sets(2).A))
	% all(all(initial_external_behavior_sets(3).b == initial_external_behavior_sets(2).b))

	%Extend the number of consistency sets by performing intersections.
	subL = n0.subL;
	[subL_powerset, subL_index_powerset] = subL.powerset();
	external_behavior_sets = initial_external_behavior_sets;
	for powerset_idx = (subL.cardinality()+1):length(subL_powerset)
		powerset_idcs_elt = subL_index_powerset{powerset_idx};
		external_behavior_sets(powerset_idx) = initial_external_behavior_sets(powerset_idcs_elt(1));
		
		for elt_idx = 2:length(powerset_idcs_elt)
			external_behavior_sets(powerset_idx) = external_behavior_sets(powerset_idx).intersect( initial_external_behavior_sets(powerset_idcs_elt(elt_idx)) );
		end
	end

	extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );

	[ ~ , empty_set_flags_eb ] = empty_bg.find_empty_observation_polyhedra( external_behavior_sets );
	[ ~ , empty_set_flags_ib ] = empty_bg.find_empty_observation_polyhedra( extended_internal_behavior_sets );

	assert( all( empty_set_flags_ib == [false(3,1);true;false;true;true] ) )

function testContainmentMatrixFunctions4(testCase)
	%testInternalBehaviorSets4.m
	%Description:
	%	This test script is meant to evaluate the function get_all_consistent_behavior_sets.m which is a member function of
	%	the belief graph object.
	%	In this example, there are different dynamics for each mode, and two of the modes are identical but one should be noticeably different.

	L1 = Language([1,2,1,2],[2,1,2,1],[3,1,1,1]);
	T = 4;

	%Create a simple Language Constrainted Switching System
	A1 = eye(2); A3 = 2*A1;
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

	f1 = [1;0];

	aff_dyn_list = [	Aff_Dyn(A1,B1, f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,-f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A3,B1, f1,C1,Pw1,Pv1) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

	%Create BeliefGraph
	empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true , 'use_proj_flag' , true );

	N0 = empty_bg.get_initial_beliefnodes();
	n0 = N0(1);

	debug_flag = 0;
	[ initial_external_behavior_sets , initial_internal_behavior_sets ] = lcsas0.get_consistency_sets_for_language( ...
												n0.t+1, n0.subL, P_u, P_x0, ...
												'fb_method','output', 'debug_flag',debug_flag, ...
												'use_proj', true , 'reduce_representation' , true );

	%Extend the number of consistency sets by performing intersections.
	subL = n0.subL;
	[subL_powerset, subL_index_powerset] = subL.powerset();
	external_behavior_sets = initial_external_behavior_sets;
	for powerset_idx = (subL.cardinality()+1):length(subL_powerset)
		powerset_idcs_elt = subL_index_powerset{powerset_idx};
		external_behavior_sets(powerset_idx) = initial_external_behavior_sets(powerset_idcs_elt(1));
		
		for elt_idx = 2:length(powerset_idcs_elt)
			external_behavior_sets(powerset_idx) = external_behavior_sets(powerset_idx).intersect( initial_external_behavior_sets(powerset_idcs_elt(elt_idx)) );
		end
	end

	extended_internal_behavior_sets = empty_bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets );

	%Create the containment matrix using the External Behavior Sets + Exact Containment
	containment_matrix_eb = false(length(subL_powerset));
	for x_idx = 1:size(containment_matrix_eb,1)
		containment_matrix_eb(x_idx,x_idx) = true;
	end


	for x_idx = 1:size(containment_matrix_eb,1)
		potential_y_idcs = [1:size(containment_matrix_eb,2)];
		potential_y_idcs = potential_y_idcs( potential_y_idcs ~= x_idx );
		for y_idx = potential_y_idcs
			%Consider the temporary combination
			%disp(['x_idx = ' num2str(x_idx) ', y_idx = ' num2str(y_idx) ])

			% Observe if Y_Set of X is contained by the Y_Set of Y
			ObservationSetX = external_behavior_sets(x_idx);
			ObservationSetY = external_behavior_sets(y_idx);
			
			containment_matrix_eb(x_idx,y_idx) = (ObservationSetX <= ObservationSetY);
		end
	end
	containment_matrix_eb;

	%Create the continment matrix using Internal Behavior Sets + Sufficient Condition
	containment_matrix_ib = empty_bg.internal_behavior_sets2containment_mat( extended_internal_behavior_sets );

	assert( all( containment_matrix_eb( containment_matrix_ib ) ) )
