function ancest_nodes = post(varargin)
	%Description:
	%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
	%	given in lcsas.
	%
	%Usage:
	%	ancest_nodes = BG.post(BN)
	%	ancest_nodes = BG.post(BN,'debug',debug_flag)
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
	%	debug_flag  - A nonnegative integer indicating the level of debug information that the user wants.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ BG , BN , post_settings ] = ip_post(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	lcsas = BG.lcsas;
	subL = BN.subL;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Select which Post function to call %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if BG.UsedProjection
		ancest_nodes = post_proj( BG , BN , post_settings );
		return
	else
		ancest_nodes = post_no_proj( BG , BN , post_settings );
		return
	end

end

function [ BG , BN , post_settings ] = ip_post(varargin)
	%Description:
	%	Input processing for the post function.

	% Checking That enough Arguments were given.

	if nargin < 4
		error(['Not enough input arguments were given to post; Expected 4, received ' num2str(nargin) '.'  ])
	end

	BG 		= varargin{1};
	BN 		= varargin{2};

	if ~isa(BN,'BeliefNode')
		error('Please make sure that the second input to post is a BeliefNode object.')
	end

	% Create Defaults

	post_settings = struct( ...
		'debug_flag', 0 , ...
		'fb_method', BG.FeedbackMethod, ...
		'accel_flag', BG.UsedAcceleratedAlgorithms , ...
		'use_unobs_checks', BG.UsedUnobservabilityChecks ...
		);

	if nargin > 4
		arg_idx = 5;
		while arg_idx <= nargin
			switch varargin{arg_idx}
				case 'debug'
					post_settings.debug_flag = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				case 'fb_method'
					post_settings.fb_method = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				case 'accel_flag'
					post_settings.accel_flag = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				case 'use_unobs_checks'
					post_settings.use_unobs_checks = varargin{arg_idx+1};
					arg_idx = arg_idx + 2;
				otherwise
					error(['Unexpected input: ' varargin{arg_idx}])
			end
		end
	end

end

function [ ancest_nodes ] = post_no_proj(BG,BN,post_settings)
	%Description:
	%
	%Usage:
	%	[ ancest_nodes ] = post_experimental(BG,BN,post_settings)
	%	
	%Assumption:
	%	This function assumes that the BeliefGraph function contains the following member variables
	%		- UsedProjection: This flag indicates whether or not the method will use projection or not.
	%		- FeedbackMethod: Indicates if the system is using state or output feedback.
	%		- UsedAcceleratedAlgorithms: Indicates if the accelerated "post" algorithms should be used.
	%		- UsedUnobservabilityChecks: Indicates that the "post" algorithm checks for unobservable edges/transitions.
	%
	%Inputs:
	%	BG 				- A Belief Graph object.
	%	lcsas 			- An array of Aff_Dyn() objects.
	%	post_settings	- A struct containing all of the special settings passed to post.

	% No input processing necessary, because this was already done by another algorithm.

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
    % use_unobs_checks = true;
    
    t_ancestor = BN.t+1; %The time of the "ancestor" of BN

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if post_settings.debug_flag > 0
		disp('  + Creating Consistency Sets.')
	end

	all_ibs = [];
	for powerset_idx = 1:length(subL_powerset)
		% Get An Element from the Powerset of subL
		powerset_elt = subL_powerset(powerset_idx);
		temp_sequence = repmat(powerset_elt,t_ancestor+1,1);

		% Get InternalBehaviorSet for this Path
		ibs_i = InternalBehaviorSet(lcsas,temp_sequence);

		% Prune / ignore certain sets based on if they are

		all_ibs = [all_ibs;ibs_i];

	end

	% %Extend the number of consistency sets by performing intersections.
	% extended_internal_behavior_sets = BG.get_all_consistent_internal_behavior_sets( all_ibs );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Check to see if the set is visible. i.e. Each set is
	%	- somehow unique compared to its 'siblings', and
	%	- not empty.

	if post_settings.debug_flag > 0
		disp('  + Detecting whether or not the Consistency Sets are independent enough to be found.')
	end

	empty_set_flags = all_ibs.IsEmpty();

	containment_matrix_ib = BG.internal_behavior_sets2containment_mat( all_ibs , ...
																		'verbosity' , post_settings.debug_flag > 1 );

	[ ~ , observation_set_is_observable ] = BG.containment_mat2observable_combos( containment_matrix_ib );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Ancestor Nodes %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	visible_transitions = [1:length(empty_set_flags)];
	if post_settings.use_unobs_checks
		visible_transitions = visible_transitions( (~empty_set_flags) & observation_set_is_observable );
	else
		visible_transitions = visible_transitions( ~empty_set_flags );	
	end

	if post_settings.debug_flag > 0
		disp('- Creating Belief Nodes.')
	end

	ancest_nodes = [];
	for trans_idx = 1:length(visible_transitions)
		temp_L = subL_powerset(visible_transitions(trans_idx));
		%temp_consist_set = consistency_sets(visible_transitions(trans_idx));
		temp_full_set = all_ibs(visible_transitions(trans_idx));

		c_node = BeliefNode(temp_L,BN.t+1,'FullTrajectorySet',temp_full_set);
		ancest_nodes = [ancest_nodes,c_node];
	end


end

function ancest_nodes = post_proj(BG,BN,post_settings)
	%Description:
	%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
	%	given in lcsas.
	%
	%Usage:
	%	ancest_nodes = post_proj(BG,BN,post_settings)
	%
	%Inputs:
	%	BG 				- A Belief Graph object.
	%	lcsas 			- An array of Aff_Dyn() objects.
	%	post_settings	- A struct containing all of the special settings passed to post.

	% No input processing should be needed here.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	lcsas = BG.lcsas;
	subL = BN.subL;

	U = lcsas.U;
	X0 = lcsas.X0;

	%Get All Combinations of the node's subset
	%node_p_set = BN.idx_powerset_of_subL();
	[subL_powerset, subL_index_powerset] = subL.powerset();

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);

	t_ancestor = BN.t+1; %The time of the "ancestor" of BN

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform Accelerated Computations If That Is Desired %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if post_settings.accel_flag
		ancest_nodes = BG.post_accel(BN,U,X0, ...
									'debug',post_settings.debug_flag, ...
									'fb_method' , post_settings.fb_method );
		return
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if post_settings.debug_flag > 0
		disp('- Creating Consistency Sets.')
	end

	[ consistency_sets , initial_PhiSets ] = lcsas.get_consistency_sets_for_language( ...
										t_ancestor,subL,U,X0, ...
										'fb_method',post_settings.fb_method, ...
										'debug_flag',post_settings.debug_flag);

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

	if post_settings.debug_flag > 0
		disp('- Detecting whether or not the Consistency Sets are independent enough to be found.')
	end


	%Check for emptiness of the consistency set.
	[ ~ , empty_set_flags ] = BG.find_empty_observation_polyhedra( consistency_sets );

	if post_settings.use_unobs_checks

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
	if post_settings.use_unobs_checks
		visible_transitions = visible_transitions( (~empty_set_flags) & observation_set_is_observable );
	else
		visible_transitions = visible_transitions( ~empty_set_flags );	
	end

	if post_settings.debug_flag > 0
		disp('- Creating Belief Nodes.')
	end

	ancest_nodes = [];
	for trans_idx = 1:length(visible_transitions)
		temp_L = subL_powerset(visible_transitions(trans_idx));
		temp_consist_set = consistency_sets(visible_transitions(trans_idx));
		temp_full_set = PhiSets(visible_transitions(trans_idx));

		c_node = BeliefNode(temp_L,t_ancestor,'ConsistencySet',temp_consist_set,'FullTrajectorySet',temp_full_set);
		ancest_nodes = [ancest_nodes,c_node];
	end

end