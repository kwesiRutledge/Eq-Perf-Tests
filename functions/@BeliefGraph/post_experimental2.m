function [ ancest_nodes ] = post_experimental2(varargin)
	%Description:
	%
	%Usage:
	%	BG.post_experimental2(BN,P_u,P_x0)
	%	BG.post_experimental2(BN,P_u,P_x0,'debug',debug_flag)
	%	BG.post_experimental2(BN,P_u,P_x0,'debug',debug_flag)
	%	
	%Assumption:
	%	This function assumes that the BeliefGraph function contains the following member variables
	%		- UsedProjection: This flag indicates whether or not the method will use projection or not.
	%		- FeedbackMethod: Indicates if the system is using state or output feedback.
	%		- UsedAcceleratedAlgorithms: Indicates if the accelerated "post" algorithms should be used.
	%		- UsedUnobservabilityChecks: Indicates that the "post" algorithm checks for unobservable edges/transitions.
	%
	%Inputs:
	%	BG 			- A Belief Graph object.
	%	lcsas 		- An array of Aff_Dyn() objects.
	%	P_u 		- A Polyhedron() object that defines the input constraints.
	%				  The input at each time must lie within the polyhedron.
	%	P_x0 		- A Polyhedron() object that defines the set of states from which the initial state is contained.
	%	debug_flag  - A nonnegative integer indicating the level of debug information that the user wants.

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

    fb_method = 'output';
    use_unobs_checks = true;
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if debug_flag > 0
		disp('  + Creating Consistency Sets.')
	end

	[ consistency_sets , initial_ib_sets ] = lcsas.get_consistency_sets_for_language( ...
												BN.t+1,subL,P_u,P_x0, ...
												'fb_method',fb_method,'debug_flag',debug_flag, ...
												'use_proj',BG.UsedProjection, ...
												'ConsistencySetVersion', 2 );

	%Extend the number of consistency sets by performing intersections.
	consistency_set_is_empty = false(subL.cardinality(),1);
	for powerset_idx = (subL.cardinality()+1):length(subL_powerset)
		powerset_idcs_elt = subL_index_powerset{powerset_idx};
		consistency_sets(powerset_idx) = consistency_sets(powerset_idcs_elt(1));
		
		for elt_idx = 2:length(powerset_idcs_elt)
			consistency_sets(powerset_idx) = consistency_sets(powerset_idx).intersect( consistency_sets(powerset_idcs_elt(elt_idx)) );
		end

		%Update list tracking emptiness
		% consistency_set_is_empty(powerset_idx) = consistency_sets(powerset_idx).isEmptySet;
	end
	%extended_internal_behavior_sets = BG.get_all_consistent_internal_behavior_sets( initial_ib_sets , 'verbosity' , debug_flag );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Check to see if the set is visible. i.e. Each set is
	%	- somehow unique compared to its 'siblings', and
	%	- not empty.

	if debug_flag > 0
		disp('  + Detecting whether or not the Consistency Sets are independent enough to be found.')
	end

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

	[ ~ , empty_set_flags ] = BG.find_empty_observation_polyhedra( consistency_sets );

	% [ ~ , observation_set_is_observable ] = BG.containment_mat2observable_combos( containment_matrix_eb );

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
		%temp_full_set = extended_internal_behavior_sets(visible_transitions(trans_idx));

		c_node = BeliefNode(temp_L,BN.t+1,'ConsistencySet',temp_consist_set);
		ancest_nodes = [ancest_nodes,c_node];
	end


end