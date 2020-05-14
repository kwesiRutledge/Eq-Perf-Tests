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
	subL_p_set = BN.subL.powerset();

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

	Phi_sets = {};
	Consist_sets = {};
	visible_transitions = [1:length(subL_p_set)];
	for p_set_idx = 1:length(subL_p_set)
		if debug_flag > 1
			disp(['  + p_set_idx = ' num2str(p_set_idx)]);
		end
		
		temp_lang = subL_p_set(p_set_idx);
		temp_BN = BeliefNode(temp_lang,BN.t+1);

		[ Consist_sets{p_set_idx} , Phi_sets{p_set_idx} ] = lcsas.consistent_set(BN.t+1,subL_p_set(p_set_idx),P_u,P_x0,'fb_method',fb_method);

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Check to see if the set is visible. i.e. Each set is
	%	- somehow unique compared to its 'parents', and
	%	- not empty.

	if debug_flag > 0
		disp('- Detecting whether or not the Consistency Sets are independent enough to be found.')
	end

	for p_set_idx = 1:length(subL_p_set)
		
		if debug_flag > 1
			disp(['  + p_set_idx = ' num2str(p_set_idx) ])
		end

		c_set_p = Consist_sets{p_set_idx};
		if use_unobs_checks
			%Search through all previous c_sets
			for prev_cset_idx = (p_set_idx+1):length(subL_p_set)
				future_c_set = Consist_sets{prev_cset_idx};
				temp_diff = c_set_p \ future_c_set;
				if length(temp_diff) == 1
					if temp_diff.isEmptySet
						visible_transitions(p_set_idx) = 0;
					end
				end
			end
		end

		if c_set_p.isEmptySet
			visible_transitions(p_set_idx) = 0;
		end 

	end

	%Process visible_transitions matrix
	visible_transitions = visible_transitions(visible_transitions ~= 0);
	visible_transitions = unique(visible_transitions);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Ancestor Nodes %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if debug_flag > 0
		disp('- Creating Belief Nodes.')
	end

	ancest_nodes = [];
	for trans_idx = 1:length(visible_transitions)
		temp_L = subL_p_set(visible_transitions(trans_idx));
		temp_consist_set = Consist_sets{visible_transitions(trans_idx)};
		temp_full_set = Phi_sets{visible_transitions(trans_idx)};

		c_node = BeliefNode(temp_L,BN.t+1,temp_consist_set,temp_full_set);
		ancest_nodes = [ancest_nodes,c_node];
	end

end