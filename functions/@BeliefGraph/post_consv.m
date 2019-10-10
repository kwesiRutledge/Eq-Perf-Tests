function ancest_nodes = post_consv(varargin)
	%Description:
	%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
	%	given in lcsas.
	%
	%Usage:
	%	BG.post(BN,P_u,P_x0)
	%	BG.post(BN,P_u,P_x0,'debug',debug_flag)
	%
	%Inputs:
	%	BG - A Belief Graph object.
	%	lcsas - An array of Aff_Dyn() objects.
	%			 

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

	lcsas = BG.lcsas;

	%Get All Combinations of the node's subset
	%node_p_set = BN.idx_powerset_of_subL();
	subL_p_set = BN.subL.powerset();

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Consistency Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Phi_sets = {};
	Consist_sets = {};
	visible_transitions = [1:length(subL_p_set)];
	for p_set_idx = 1:length(subL_p_set)
		disp(['p_set_idx = ' num2str(p_set_idx)])
		
		temp_lang = subL_p_set(p_set_idx);
		temp_BN = BeliefNode(temp_lang,BN.t+1);
		if BG.find_node_idx(temp_BN) ~= -1
			%Does the belief node already exist in the tree?
			%If so, then get the consistency set from there.
			Consist_sets{p_set_idx} = BG.N( BG.find_node_idx(temp_BN) ).c_set;
		%Check first to see if any of the words in this sublanguage are too short to create something at time t+1
		elseif subL_p_set(p_set_idx).find_shortest_length() >= BN.t + 1
			if subL_p_set(p_set_idx).cardinality() == 1
				[ Consist_sets{p_set_idx} , ~ ] = lcsas.consistent_set(BN.t+1,subL_p_set(p_set_idx),P_u,P_x0,'fb_method',fb_method);
				%If any of the Phi's are empty,
				%then it is impossible for a transition to exist between the node c_level(node_ind) and the node associated with Phi
				if Consist_sets{p_set_idx}.isEmptySet
					visible_transitions(p_set_idx) = 0;
				end
			else
				%If there is more than one word in the element of the powerset, then do not call the consistent set function.
				%Instead use intersections to construct Consistency Set.

				temp_c_set_arr = [];

				%Find all p_set indices that match an element in the set.
				for c_set_idx = 1:BN.subL.cardinality()
					if subL_p_set(p_set_idx).contains( BN.subL.words{c_set_idx} )
						temp_c_set_arr = [temp_c_set_arr,Consist_sets{c_set_idx}];
					end
				end

				%Intersect this array together.
				temp_c_set = temp_c_set_arr(1);
				for c_set_idx = 2:length(temp_c_set_arr)
					temp_c_set = temp_c_set.intersect(temp_c_set_arr(c_set_idx));
				end

				%Save C_set.
				Consist_sets{p_set_idx} = temp_c_set;

				%If any of the Phi's are empty,
				%then it is impossible for a transition to exist between the node c_level(node_ind) and the node associated with Phi
				if Consist_sets{p_set_idx}.isEmptySet
					visible_transitions(p_set_idx) = 0;
				end
			end
		else
			visible_transitions(p_set_idx) = 0;
		end
			
	end 

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Identify Which Consistency Sets Can Be Independently Detected %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% %For each possible transition, see if its transition set is completely contained by another transition set
	% for ut_idx = 1:length(subL_p_set)
	% 	if debug_flag >= 1
	% 		disp(['ut_idx = ' num2str(ut_idx)])
	% 	end
	% 	if visible_transitions(ut_idx) ~= 0
	% 		for ch_idx = [ut_idx+1:length(subL_p_set)]
	% 			if debug_flag >= 1
	% 				disp(['ch_idx = ' num2str(ch_idx)])
	% 			end
	% 			temp_diff = Consist_sets{ut_idx} \ Consist_sets{ch_idx};
	% 			%Check to see if temp_diff is a single Polyhedron/Polytope (if it is multiple then the set difference is not empty)
	% 			if length(temp_diff) == 1
	% 				if (temp_diff.isEmptySet) && (~(Consist_sets{ut_idx}.isEmptySet)) && (~(Consist_sets{ch_idx}.isEmptySet))
	% 					if debug_flag >= 2
	% 						disp(['temp_diff.isEmptySet = ' num2str(temp_diff.isEmptySet) ' for:'])
	% 						disp(['- Phi1(' num2str([subL_p_set(ut_idx).words{:}]) ')' ])
	% 						disp(['- Phi2(' num2str([subL_p_set(ch_idx).words{:}]) ')' ])
	% 						disp(' ')
	% 					end
	% 					visible_transitions(ut_idx) = 0;
	% 				% elseif (ch_idx == ind_ut) && (~Projx_Phi1.isEmptySet)
	% 				% 	visible_transitions(ind_ut) = ch_idx;
	% 				end 
	% 			end
	% 		end
	% 	end
	% end
	%Process visible_transitions matrix
	visible_transitions = visible_transitions(visible_transitions ~= 0);
	visible_transitions = unique(visible_transitions);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Ancestor Nodes %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ancest_nodes = [];
	for trans_idx = 1:length(visible_transitions)
		temp_L = subL_p_set(visible_transitions(trans_idx));
		temp_consist_set = Consist_sets{visible_transitions(trans_idx)};

		c_node = BeliefNode(temp_L,BN.t+1,temp_consist_set);
		ancest_nodes = [ancest_nodes,c_node];
	end

end