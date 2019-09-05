classdef BeliefNode

	properties
		subL;
		t;
	end

	methods
		function BN = BeliefNode(subset_L,t0)
			
			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if ~(iscell(subset_L) || isa(subset_L,'Language') )
				error('Expected subset of language to be a cell array or a Language object.')
			end

			if iscell(subset_L)
				warning('The version of this class that uses L as a cell array is deprecated. We cannot guarantee that functions will work anymore as the class has been modified to support Language objects.')
				%Convert subset_L to a Language
				temp_arr = subset_L;
				clear subset_L
				subset_L = Language();
				subset_L.words = temp_arr;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Assigning Variables %%
			%%%%%%%%%%%%%%%%%%%%%%%%%

			BN.subL = subset_L;
			BN.t = t0;
		end

		function subsets = idx_powerset_of_subL(obj)
			%Description:
			%	Assigns to each word in subL.words an index.
			%	Then returns all possible subsets of the INDICES.
			%
			
			node_p_set = {};
			for comb_length = 1:length(obj.subL.words)
				temp_combs = nchoosek([1:length(obj.subL.words)],comb_length);
				for comb_ind = 1:size(temp_combs,1)
					node_p_set{end+1} = temp_combs(comb_ind,:);
				end
			end
			subsets = node_p_set;
		end

		function ancest_nodes = post(varargin)
			%Description:
			%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
			%	given in ad_arr.
			%
			%Usage:
			%	post(BN,ad_arr,P_u,P_x0)
			%	post(BN,ad_arr,P_u,P_x0,'debug',debug_flag)
			%
			%Inputs:
			%	ad_arr - An array of Aff_Dyn() objects.
			%			 

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if nargin < 4
				error('Not enough input arguments.')
			end

			BN 		= varargin{1};
			ad_arr 	= varargin{2};
			P_u 	= varargin{3};
			P_x0 	= varargin{4};

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

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			if ~exist('debug_flag')
				debug_flag = 0;
			end

			%Get All Combinations of the node's subset
			node_p_set = BN.idx_powerset_of_subL();

			n = size(ad_arr(1).A,1);
			m = size(ad_arr(1).B,2);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Define Consistency Sets %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			Phi_sets = {}; visible_transitions = [];
			for p_set_idx = 1:length(node_p_set)
				Phi_sets{p_set_idx} = consistency_set(ad_arr,(BN.t+1),{BN.subL.words{node_p_set{p_set_idx}}},P_u,P_x0);
				%If any of the Phi's are empty,
				%then it is impossible for a transition to exist between the node c_level(node_ind) and the node associated with Phi
				if ~Phi_sets{p_set_idx}.isEmptySet
					visible_transitions = [visible_transitions p_set_idx];
				end
			end 

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Identify Which Consistency Sets Can Be Independently Detected %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			%For each possible transition, see if its transition set is completely contained by another transition set
			for ut_idx = 1:length(node_p_set)
				Projx_Phi1 = [eye(n*((BN.t+1)+1)+m*(BN.t+1)) zeros(n*((BN.t+1)+1)+m*(BN.t+1),Phi_sets{ut_idx}.Dim-n*((BN.t+1)+1)-m*(BN.t+1))]*Phi_sets{ut_idx};
				%vt_tilde = visible_transitions(visible_transitions ~= ind_ut);
				for ch_idx = [ut_idx+1:length(node_p_set)]
					if debug_flag >= 1
						disp(['ch_idx = ' num2str(ch_idx)])
					end
					Projx_Phi2 = [eye(n*((BN.t+1)+1)+m*(BN.t+1)) zeros(n*((BN.t+1)+1)+m*(BN.t+1),Phi_sets{ch_idx}.Dim-n*((BN.t+1)+1)-m*(BN.t+1))]*Phi_sets{ch_idx};
					temp_diff = Projx_Phi1 \ Projx_Phi2;
					if (temp_diff.isEmptySet)
						if debug_flag >= 1
							disp(['temp_diff.isEmptySet = ' num2str(temp_diff.isEmptySet) ' for:'])
							disp(['- Phi1(' num2str([BN.subL.words{node_p_set{ut_idx}}]) ')' ])
							disp(['- Phi2(' num2str([BN.subL.words{node_p_set{ch_idx}}]) ')' ])
							disp(' ')
						end
						visible_transitions(ut_idx) = ch_idx;
					% elseif (ch_idx == ind_ut) && (~Projx_Phi1.isEmptySet)
					% 	visible_transitions(ind_ut) = ch_idx;
					end 
				end
			end
			visible_transitions = unique(visible_transitions);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Create Ancestor Nodes %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%
			ancest_nodes = [];
			for trans_idx = 1:length(visible_transitions)
				temp_L = Language();
				temp_L.words = {BN.subL.words{node_p_set{visible_transitions(trans_idx)}}};

				c_node = BeliefNode(temp_L,BN.t+1);
				ancest_nodes = [ancest_nodes,c_node];
			end

		end

		function eq_flag = is_eq(obj,BN)
			%Description:
			%	Returns true if this belief node (obj) is equal to the belief node BN.

			if obj.subL.is_eq(BN.subL) && (obj.t == BN.t)
				eq_flag = true;
			else
				eq_flag = false;
			end

		end

		function longest_T = find_longest_horizon(obj)
			%Description:
			%	Searches through all elements of the subL.words for this node and determines
			%	how much longer of a future the system can have (at maximum).

			longest_T = -1;

			for L_idx = 1:length(obj.subL.words)
				longest_T = max(length(obj.subL.words{L_idx}),longest_T) - obj.t;
			end
		end
	end

end