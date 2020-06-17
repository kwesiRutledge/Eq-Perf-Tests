function [ ancest_nodes ] = post_experimental(varargin)
	%Description:
	%
	%Usage:
	%	BG.post_experimental(BN,P_u,P_x0)
	%	BG.post_experimental(BN,P_u,P_x0,'debug',debug_flag)
	%	BG.post_experimental(BN,P_u,P_x0,'debug',debug_flag,'')
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

	[ ~ , phi_sets ] = lcsas.get_consistency_sets_for_language( ...
										BN.t+1,subL,P_u,P_x0, ...
										'fb_method',fb_method,'debug_flag',debug_flag, ...
										'use_proj',BG.UsedProjection );

	%Extend the number of consistency sets by performing intersections.
	extended_internal_behavior_sets = BG.get_all_consistent_internal_behavior_sets( phi_sets );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Check to see if the set is visible. i.e. Each set is
	%	- somehow unique compared to its 'siblings', and
	%	- not empty.

	if debug_flag > 0
		disp('  + Detecting whether or not the Consistency Sets are independent enough to be found.')
	end

	[ ~ , empty_set_flags ] = BG.find_empty_observation_polyhedra( extended_internal_behavior_sets );

	containment_matrix_ib = BG.internal_behavior_sets2containment_mat( extended_internal_behavior_sets , 'verbosity' , debug_flag > 1 );

	[ ~ , observation_set_is_observable ] = BG.containment_mat2observable_combos( containment_matrix_ib );

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
		%temp_consist_set = consistency_sets(visible_transitions(trans_idx));
		temp_full_set = extended_internal_behavior_sets(visible_transitions(trans_idx));

		c_node = BeliefNode(temp_L,BN.t+1,'FullTrajectorySet',temp_full_set);
		ancest_nodes = [ancest_nodes,c_node];
	end


end