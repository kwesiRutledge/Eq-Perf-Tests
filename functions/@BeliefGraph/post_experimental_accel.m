function ancest_nodes = post_experimental_accel(varargin)
	%Description:
	%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
	%	given in lcsas.
	%
	%Usage:
	%	BG.post_experimental_accel(BN,P_u,P_x0)
	%	BG.post_experimental_accel(BN,P_u,P_x0,'debug',debug_flag)
	%	BG.post_experimental_accel(BN,P_u,P_x0,'debug',debug_flag, 'fb_method' , fb_method )
	%
	%Inputs:
	%	BG 			- A Belief Graph object.
	%	lcsas 		- An array of Aff_Dyn() objects.
	%	P_u 		- A Polyhedron() object that defines the input constraints.
	%				  The input at each time must lie within the polyhedron.
	%	P_x0 		- A Polyhedron() object that defines the set of states from which the initial state is contained.
	%	fb_method 	- A string representing which feedback method is being used for this example.
	%				  Options are: 'state' or 'output'
	%	accel_flag	- A boolean which is used to decide whethor or not to use the accelerated or simple version of the algorithm.
	%				  Must be true or false.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 4
		error('Not enough input arguments.')
	end

	BG 		= varargin{1};
	BN 		= varargin{2};
	P_u 	= varargin{3};
	P_x0 	= varargin{4};

	if ~isa(BN,'BeliefNode')
		error('Please make sure that the second input to post is a BeliefNode object.')
	end

	if nargin > 4
		arg_idx = 5;
		while arg_idx <= nargin
			switch varargin{arg_idx}
				case 'debug'
					debug_flag = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				case 'fb_method'
					fb_method = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				case 'PerformBoundingBoxCheck'
					PerformBoundingBoxCheck = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				otherwise
					error(['Unexpected input: ' varargin{arg_idx}])
			end
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	if ~exist('debug_flag')
		debug_flag = 0;
	end

	if ~exist('fb_method')
		fb_method = 'output';
	end

	if ~exist('PerformBoundingBoxCheck')
		PerformBoundingBoxCheck = false;
	end

	lcsas = BG.lcsas;
	subL = BN.subL;

	%Get All Combinations of the node's subset
	%node_p_set = BN.idx_powerset_of_subL();
	subL_p_set = BN.subL.powerset();

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);
	n_y = size(lcsas.Dyn(1).C,1);

	t_p = BN.t+1;
	external_beh_dim = n_y*(t_p+1)+m*t_p;

	cg = constr_gen(0);
	ops0 = sdpsettings(	'verbose',debug_flag,'cachesolvers',1,...
						'solver','gurobi', 'gurobi.BarIterLimit', 15000);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[~,initial_ib_sets] = lcsas.get_consistency_sets_for_language(t_p,subL,P_u,P_x0, ...
																	'use_proj', false, ...
																	'fb_method','output');
	ib_sets = BG.get_all_consistent_internal_behavior_sets( initial_ib_sets );
	[ ~ , empty_set_flags ] = BG.find_empty_observation_polyhedra( ib_sets );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	containment_matrix = false(length(subL_p_set));
	for ut_idx = 1:size(containment_matrix,1)
		containment_matrix(ut_idx,ut_idx) = true;
	end

	%For each possible transition, see if its transition set is completely contained by another transition set
	for cm_idx1 = 1:size(containment_matrix,1)
		% L_ut = subL_p_set(ut_idx);
		if debug_flag >= 1
			disp(['  * cm_idx1 = ' num2str(cm_idx1)])
		end

		%Only check the indices of the powerset that have NON-empty
		%consistency sets.
		if empty_set_flags(cm_idx1)
			continue;
		end
		
		for cm_idx2 = [cm_idx1+1:size(containment_matrix,2)]

			if debug_flag >= 1
				disp(['  * cm_idx2 = ' num2str(cm_idx2)])
			end

			%if cm_idx2 is empty, then don't consider it.
			if empty_set_flags(cm_idx2)
				continue;
			end

			%Set up optimization
			ib_set_x = ib_sets(cm_idx1);
			ib_set_y = ib_sets(cm_idx2);

			%if the bounding box of ib_set_x is not included by the bounding box of ib_set_y then skip this computation.
			if PerformBoundingBoxCheck
				ib_set_x.outerApprox(); ib_set_y.outerApprox();
				if any(ib_set_x.Internal.ub([1:external_beh_dim]) < ib_set_y.Internal.lb([1:external_beh_dim])) | ...
				   any(ib_set_x.Internal.lb([1:external_beh_dim]) > ib_set_y.Internal.ub([1:external_beh_dim]))
					continue;
				end
			end

			dim_ibX = ib_set_x.Dim;
			dim_ibY = ib_set_y.Dim;

			Rt_X = [ eye(external_beh_dim), zeros(external_beh_dim, dim_ibX - external_beh_dim) ];
			Rt_Y = [ eye(external_beh_dim), zeros(external_beh_dim, dim_ibY - external_beh_dim) ];

			%Perform Containment Check via Optimization
			[dual_vars, constr] = cg.create_sadraddini_AH_inclusion_constr( ...
									zeros(external_beh_dim,1) , Rt_X , [ib_set_x.A;ib_set_x.Ae;-ib_set_x.Ae] , [ib_set_x.b;ib_set_x.be;-ib_set_x.be] , ...
									zeros(external_beh_dim,1) , Rt_Y , [ib_set_y.A;ib_set_y.Ae;-ib_set_y.Ae] , [ib_set_y.b;ib_set_y.be;-ib_set_y.be] );

			diagnostics = optimize(constr,[],ops0);

			containment_matrix(cm_idx1,cm_idx2) = (diagnostics.problem == 0);

            if containment_matrix(cm_idx1,cm_idx2)
                if debug_flag >= 2
                    disp(['temp_diff.isEmptySet = for:'])
                    disp(['- Phi1(' num2str([subL_p_set(cm_idx1).words{:}]) ')' ])
                    disp(['- Phi2(' num2str([subL_p_set(cm_idx2).words{:}]) ')' ])
                    disp(' ')
                end
                
                break;
            end 
		end
	end

	[~,observable_flags] = BG.containment_mat2observable_combos(containment_matrix);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Process visible_transitions matrix %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	visible_transitions = [1:length(ib_sets)];
	visible_transitions = visible_transitions( (~empty_set_flags) & observable_flags );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Ancestor Nodes %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ancest_nodes = [];
	for trans_idx = 1:length(visible_transitions)
		temp_L = subL_p_set(visible_transitions(trans_idx));
		%temp_consist_set = Consist_sets{visible_transitions(trans_idx)};
		temp_ib_set = ib_sets(trans_idx);

		c_node = BeliefNode(temp_L,BN.t+1,'FullTrajectorySet',temp_ib_set);
		ancest_nodes = [ancest_nodes,c_node];
	end

end