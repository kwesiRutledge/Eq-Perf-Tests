function ancest_nodes = post_proj(varargin)
	%Description:
	%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
	%	given in lcsas.
	%
	%Usage:
	%	BG.post_proj(BN,P_u,P_x0)
	%	BG.post_proj(BN,P_u,P_x0,'debug',debug_flag)
	%	BG.post_proj(BN,P_u,P_x0,'debug',debug_flag, 'fb_method' , fb_method )
	%	BG.post_proj(BN,P_u,P_x0,'debug',debug_flag, 'fb_method' , fb_method , 'accel_flag' , accel_flag )
	%	BG.post_proj(BN,P_u,P_x0,'use_unobservability_checks',false)
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
		error('Please make sure that the second input to post_proj is a BeliefNode object.')
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
				case 'accel_flag'
					accel_flag = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				case 'use_unobservability_checks'
					use_unobs_checks = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				otherwise
					error(['Unexpected input: ' varargin{arg_idx}])
			end
		end
	end

	%%%%%%%%%%%%%%
	%% Defaults %%
	%%%%%%%%%%%%%%

	if ~exist('debug_flag')
		debug_flag = 0;
	end

	if ~exist('fb_method')
		fb_method = 'output';
	end

	if ~exist('accel_flag')
		accel_flag = false;
	end

	if ~exist('use_unobs_checks')
		use_unobs_checks = true;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	lcsas = BG.lcsas;
	subL = BN.subL;

	%Get All Combinations of the node's subset
	%node_p_set = BN.idx_powerset_of_subL();
	[subL_powerset, subL_index_powerset] = subL.powerset();

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform Accelerated Computations If That Is Desired %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if accel_flag
		ancest_nodes = BG.post_accel(BN,P_u,P_x0,'debug',debug_flag, 'fb_method' , fb_method );
		return
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if debug_flag > 0
		disp('- Creating Consistency Sets.')
	end

	[ consistency_sets , initial_PhiSets ] = lcsas.get_consistency_sets_for_language( ...
										BN.t+1,subL,P_u,P_x0, ...
										'fb_method',fb_method,'debug_flag',debug_flag);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Extend the number of consistency sets by performing intersections.
	for powerset_idx = (subL.cardinality()+1):length(subL_powerset)
		powerset_idcs_elt = subL_index_powerset{powerset_idx};
		consistency_sets(powerset_idx) = consistency_sets(powerset_idcs_elt(1));
		
		for elt_idx = 2:length(powerset_idcs_elt)
			consistency_sets(powerset_idx) = consistency_sets(powerset_idx).intersect( consistency_sets(powerset_idcs_elt(elt_idx)) );
		end
	end

	PhiSets = BG.get_all_consistent_internal_behavior_sets( initial_PhiSets );

	%Check to see if the set is visible. i.e. Each set is
	%	- somehow unique compared to its 'siblings', and
	%	- not empty.

	if debug_flag > 0
		disp('- Detecting whether or not the Consistency Sets are independent enough to be found.')
	end


	%Check for emptiness of the consistency set.
	[ ~ , empty_set_flags ] = BG.find_empty_observation_polyhedra( consistency_sets );

	if use_unobs_checks

		containment_matrix = false(length(subL_powerset));
		for x_idx = 1:size(containment_matrix,1)
			containment_matrix(x_idx,x_idx) = true;
		end


		for x_idx = 1:size(containment_matrix,1)
			potential_y_idcs = [1:size(containment_matrix,2)];
			potential_y_idcs = potential_y_idcs( potential_y_idcs ~= x_idx );
			for y_idx = potential_y_idcs
				%Consider the temporary combination
				%disp(['x_idx = ' num2str(x_idx) ', y_idx = ' num2str(y_idx) ])

				% Observe if Y_Set of X is contained by the Y_Set of Y
				ObservationSetX = consistency_sets(x_idx);
				ObservationSetY = consistency_sets(y_idx);
				
				containment_matrix(x_idx,y_idx) = (ObservationSetX <= ObservationSetY);
			end
		end

		[ ~ , observation_set_is_observable ] = BG.containment_mat2observable_combos( containment_matrix );

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Ancestor Nodes %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	visible_transitions = [1:length(empty_set_flags)];
	if use_unobs_checks
		visible_transitions = visible_transitions( (~empty_set_flags) & observation_set_is_observable );
	else
		visible_transitions = visible_transitions( ~empty_set_flags );	
	end

	if debug_flag > 0
		disp('- Creating Belief Nodes.')
	end

	ancest_nodes = [];
	for trans_idx = 1:length(visible_transitions)
		temp_L = subL_powerset(visible_transitions(trans_idx));
		temp_consist_set = consistency_sets(visible_transitions(trans_idx));
		temp_full_set = PhiSets(visible_transitions(trans_idx));

		c_node = BeliefNode(temp_L,BN.t+1,'ConsistencySet',temp_consist_set,'FullTrajectorySet',temp_full_set);
		ancest_nodes = [ancest_nodes,c_node];
	end

end