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

			if ~iscell(subset_L)
				error('Expected subset of language to be a cell array.')
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Assigning Variables %%
			%%%%%%%%%%%%%%%%%%%%%%%%%

			BN.subL = subset_L;
			BN.t = t0;
		end

		function subsets = idx_powerset_of_subL(obj)
			%Description:
			%	Assigns to each word in subL an index.
			%	Then returns all possible subsets of the INDICES.
			%
			
			node_p_set = {};
			for comb_length = 1:length(obj.subL)
				temp_combs = nchoosek([1:length(obj.subL)],comb_length);
				for comb_ind = 1:size(temp_combs,1)
					node_p_set{end+1} = temp_combs(comb_ind,:);
				end
			end
			subsets = node_p_set;
		end

		function ancest_nodes = post(obj,ad_arr,P_u,P_x0)
			%Description:
			%	Identifies what nodes could possibly arise after reaching the current Belief Node according to the dynamics
			%	given in ad_arr.
			%
			%Inputs:
			%	ad_arr - An array of Aff_Dyn() objects.
			%			 

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			debug_flag = 1;

			%Get All Combinations of the node's subset
			node_p_set = obj.idx_powerset_of_subL();

			n = size(ad_arr(1).A,1);
			m = size(ad_arr(1).B,2);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Define Consistency Sets %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			Phi_sets = {}; visible_transitions = [];
			for p_set_idx = 1:length(node_p_set)
				Phi_sets{p_set_idx} = consistency_set(ad_arr,(obj.t+1),{obj.subL{node_p_set{p_set_idx}}},P_u,P_x0);
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
				Projx_Phi1 = [eye(n*((obj.t+1)+1)+m*(obj.t+1)) zeros(n*((obj.t+1)+1)+m*(obj.t+1),Phi_sets{ut_idx}.Dim-n*((obj.t+1)+1)-m*(obj.t+1))]*Phi_sets{ut_idx};
				%vt_tilde = visible_transitions(visible_transitions ~= ind_ut);
				for ch_idx = [ut_idx+1:length(node_p_set)]
					if debug_flag >= 1
						disp(['ch_idx = ' num2str(ch_idx)])
					end
					Projx_Phi2 = [eye(n*((obj.t+1)+1)+m*(obj.t+1)) zeros(n*((obj.t+1)+1)+m*(obj.t+1),Phi_sets{ch_idx}.Dim-n*((obj.t+1)+1)-m*(obj.t+1))]*Phi_sets{ch_idx};
					temp_diff = Projx_Phi1 \ Projx_Phi2;
					if (temp_diff.isEmptySet)
						if debug_flag >= 1
							disp(['temp_diff.isEmptySet = ' num2str(temp_diff.isEmptySet) ' for:'])
							disp(['- Phi1(' num2str([obj.subL{node_p_set{ut_idx}}]) ')' ])
							disp(['- Phi2(' num2str([obj.subL{node_p_set{ch_idx}}]) ')' ])
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
				temp_L = {obj.subL{node_p_set{visible_transitions(trans_idx)}}};

				c_node = BeliefNode(temp_L,obj.t+1);
				ancest_nodes = [ancest_nodes,c_node];
			end

		end

		function eq_flag = is_eq(obj,BN)
			%Description:
			%	Returns true if this belief node (obj) is equal to the belief node BN.

			if language_eq(obj.subL,BN.subL) && (obj.t == BN.t)
				eq_flag = true;
			else
				eq_flag = false;
			end

		end
	end

end