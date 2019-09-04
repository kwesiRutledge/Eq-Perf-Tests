classdef BeliefGraph
	%Description:
	%
	%Member Variables:
	%	E 	- An n_E x 2 matrix of strictly positive integers, where each row represents an edge in the graph.
	%		  For a row E_i, the first entry (first column) is the INDEX of the edge's source node (i.e. N(E_i(1)) is
	%		  the source node of the edge) and the second entry (second column) is the INDEX of the edge's
	%		  destination node (i.e. N(E_i(2)) is the destination node of the edge).

	properties
		E;
		N;
		LCSS;
	end

	methods
		function [BG] = BeliefGraph(lcss,L,P_u,P_x0)
			% disp('Created empty Belief Graph. Please call the construct() function next.')
			% BT.E = [];
			% BT.N = [];

			%%
			warning('The BeliefGraph() class was designed with the assumption that the affine dynamics contain the same state, input, and disturbance size throughout the entire language.')

			%Constants

			BG.LCSS = lcss;

			%Create first node
			% node0.subset = L;
			% node0.t = 0;
			node0 = BeliefNode(L,0);
			T_max = node0.find_longest_horizon();

			c_level = [node0];
			BG.N = [node0];

			% nodes0 = [node0];
			% nodes1 = [node1];
			% edges0 = [];
			% edges1 = []; %Numerical version of the edges matrix

			for tau = 1:T_max
				%Each belief will be indexed by a time. (i.e. I hold X belieft at time t)
				for node_idx = 1:length(c_level) %Iterate through all nodes that are stored in the c_level array
					%Current node
					c_node = c_level(node_idx);

					%Calculate the ancestors of this Belief Node
					temp_post = c_node.post(BG.LCSS,P_u,P_x0);
					
					%Add the ancestors to the BeliefGraph's set of nodes if they don't already exist in the set.
					for node_idx = 1:length(temp_post)
						if BG.find_node_idx(temp_post(node_idx)) == -1
							BG.N = [BG.N,temp_post(node_idx)];
						end
					end

					% Create edges corresponding to the ancestor-descendant relationship
					for edge_idx = 1:length(temp_post)
						%Create edge using this new node and add it to the edges list
						temp_edge = [BG.find_node_idx(c_node), ...
									 BG.find_node_idx(temp_post(edge_idx))];
						BG.E = [BG.E;temp_edge];
					end

				end

				%Create next level of the tree
				c_level = BG.get_all_nodes_at_time(tau);

				disp(['There are ' num2str(length(c_level)) ' nodes at time tau = ' num2str(tau) '.' ])
			end

		end

		function [BN_idx] = find_node_idx(obj,BN)
			%Description:
			%	Identifies the index of the belief node in the set N and returns that index.
			%	If the node is not found in the set N then, this returns 0.

			BN_idx = -1;

			for node_idx = 1:length(obj.N)
				if BN.is_eq(obj.N(node_idx))
					BN_idx = node_idx;
				end 
			end
		end

		function [subN] = get_all_nodes_at_time(obj,tau)
			%Description:
			%	Retrieves all belief nodes that have time value provided.

			subN = [];
			for node_idx = 1:length(obj.N)
				temp_BN = obj.N(node_idx);
				if temp_BN.t == tau
					subN = [subN,temp_BN];
				end
			end

		end

		function [] = plot(obj)
			%Description:
			%	Uses MATLAB's Built In Graph Theory tools to visualize the belief graph.

			%% Constants

			%% Algorithm
			g0 = digraph(obj.E(:,1),obj.E(:,2));
			plot(g0)

		end

		function [leaf_node_idxs] = get_leaf_node_idxs(obj)
			%Description:
			%	Returns the indexes of the 'leaf' nodes in the graph.

			leaf_node_idxs = [1:length(obj.N)];

			for cand_leaf_idx = 1:length(obj.N)
				if any(obj.E(:,1) == cand_leaf_idx)
					leaf_node_idxs = leaf_node_idxs(leaf_node_idxs ~= cand_leaf_idx)
				end
			end
		end

		function [predecessor_idxs] = pre(varargin)
			%Description:
			%
			%Outputs:
			%	predecessor_idxs 	- Array 
			%						  Returns empty array if the input is the root node reached.

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if nargin < 2
				error('Improper number of input arguments. Need at least 2.')
			end

			obj = varargin{1};
			if isa(varargin{2},'BeliefNode')
				BN = varargin{2};
				BN_idx = obj.get_node_idx(BN);
			elseif isscalar(varargin{2})
				BN_idx = varargin{2};
				if (BN_idx < 1) || (BN_idx > length(obj.N))
					error('BeliefNode index appears to be out of range.')
				end
			else
				error('Unrecognized data type. Please provide a Belief Node or node index to the function.')
			end	

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Get all edges that have the given node as the destination %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			dest_edges = (obj.E(:,2) == BN_idx);
			predecessor_idxs = obj.E(dest_edges,1);

		end

		function [reached_root] = all_words_start_with_root(obj,L)
			%Description:
			%	This function returns true if all words in the language begin with L.

			reached_root = true;
			for word_idx = 1:length(L)
				reached_root = reached_root && (L{word_idx}(1) == 1);
			end
		end

		function [new_lang] = prepend_any_valid_node(obj,word)
			%Description:
			%	For the word 'word', look at the initial index which corresponds to the initial node in the word.
			%	For all ancestors of that node, create a new word that starts with the ancestor and is followed by
			%	the original word.
			%
			%Outputs:
			%	new_lang 	- A cell array which contains 0 elements if there are no possible completions. Otherwise, there
			%				  should be a cell array of numeric arrays in new_lang.

			new_lang = {};

			temp_pre = obj.pre(word(1));
			if ~isempty(temp_pre)
				%Create a 
				for pre_idx = 1:length(temp_pre)
					new_lang{pre_idx} = [temp_pre(pre_idx),word];
				end
			end

		end


		function [new_lang] = get_one_step_extension(obj,in_lang)
			%Description:
			%	Takes every word in the input language and 'extends' it backwards by appending to it
			%	a feasible prefix of a single symbol if possible. 

			%% Variable Definitions

			%% Algorithm

			for word_idx = 1:length(in_lang)
				%
			end

		end


		function [belief_lang,leaf_node_idxs] = get_belief_language(obj)
			%Description:
			%	Looks at the leaf nodes of the graph (nodes that do not have any descendants) and does
			%	backtracking to create a language of beliefs.

			%% Initialize variables

			belief_lang = [];

			%% Find all leaves
			leaf_node_idxs = obj.get_leaf_node_idxs(obj)

			%% Perform backtracking operation for every node 
			belief_lang = {};
			for leaf_idx = leaf_node_idxs
				sub_lang = {[leaf_idx]};
				while(~obj.all_words_start_with_root())
					%Until all of the words in this sub_language have reached the root node $n$ 
					for word_idx = 1:length(sub_lang)
						%For every word in the current sub language, expand the word by getting the pre.
						temp_word = sub_lang{word_idx};
						
					end
				end

			end


		end

	end


end