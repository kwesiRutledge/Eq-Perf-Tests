function [results] = observer_comparison78( varargin )
	%observer_comparison78.m
	%Description:
	%	Testing the projection method.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	% if nargin >= 1
	% 	load_data_flag = varargin{1};
	% end

	% if nargin >= 3
	% 	c_sq.dim_x = varargin{2};
	% 	c_sq.dim_y = varargin{3};
	% end

	% if nargin >= 4
	% 	verbosity = varargin{4};
	% end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison78';
	save_file_name = [ 'results/' test_name '_results.mat'];

	n_x = 3;

	unit_box = Polyhedron('lb',-ones(1,n_x),'ub',ones(1,n_x));

	%%%%%%%%%%%%%
	%% Testing %%
	%%%%%%%%%%%%%

	%  ==========================================
	%% Create Simple Polyhedron based on unit_box

	unbounded_line_seqment = Polyhedron('Ae',[1,-2,-0.5],'be',0);
	X = unit_box.intersect(unbounded_line_seqment);

	results.X = X;

	figure;
	subplot(2,1,1)
	plot(X)
	subplot(2,1,2)
	plot(X.projection([1,2]))

	%  =================================================
	%% Create Consistency Sets for Systems Which I Know!

	L1 = Language([1,2,1,2],[2,1,2,1],[1,1,1,1]);
	T = 4;

	%Create a simple Language Constrainted Switching System
	A1 = 1;
	B1 = 1;
	C1 = 1;

	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = 1;

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

	all(all(initial_external_behavior_sets(1).A == initial_external_behavior_sets(2).A))
	[initial_external_behavior_sets(1).A,initial_external_behavior_sets(1).b]
	[initial_internal_behavior_sets(1).A,initial_internal_behavior_sets(1).b]
	[initial_internal_behavior_sets(1).Ae,initial_internal_behavior_sets(1).be]
	all(all(initial_external_behavior_sets(1).b == initial_external_behavior_sets(2).b))
	all(all(initial_external_behavior_sets(3).A == initial_external_behavior_sets(2).A))
	all(all(initial_external_behavior_sets(3).b == initial_external_behavior_sets(2).b))

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

	[ ~ , empty_set_flags_eb ] = empty_bg.find_empty_observation_polyhedra( external_behavior_sets )
	[ ~ , empty_set_flags_ib ] = empty_bg.find_empty_observation_polyhedra( extended_internal_behavior_sets )

	%assert( all( empty_set_flags_ib == empty_set_flags_eb ) )

	figure;
	for init_eb_idx = 1:n0.subL.cardinality()
		subplot(3,1,init_eb_idx)
		plot( initial_external_behavior_sets(init_eb_idx).projection([1,2]) )
	end

	temp_set = initial_external_behavior_sets(1).intersect( ...
				initial_external_behavior_sets(2));
	figure;
	plot(temp_set.projection([1,2]))


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%save([ save_file_name '.mat'])

end