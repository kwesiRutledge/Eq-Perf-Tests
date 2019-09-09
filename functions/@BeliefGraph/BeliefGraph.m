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
		lcsas;
		L;
	end

	methods
		function [BG] = BeliefGraph(in_sys,L,P_u,P_x0)
			%Description:
			%
			%Inputs:
			%	in_sys 	- An array of Aff_Dyn() objects. May eventually become its own class/data type soon.

			% disp('Created empty Belief Graph. Please call the construct() function next.')
			% BT.E = [];
			% BT.N = [];

			%%
			warning('The BeliefGraph() class was designed with the assumption that the affine dynamics contain the same state, input, and disturbance size throughout the entire language.')

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if ~isa(L,'Language')
				error('Expected L to be a Language object.')
			end

			if ~isa(in_sys,'LCSAS')
				warning('Expected input system to be LCSAS. Checking to see if this can be salvaged...')
				if isa(in_sys,'Aff_Dyn')
					temp_sys_arr = in_sys;
					in_sys = LCSAS(temp_sys_arr);
				else
					error('System must be given either as an LCSAS object or as an array of Aff_Dyn objects.')
				end
			end

			%%Constants

			BG.lcsas = in_sys;
			BG.L = L;

			%Create first node
			% node0.subset = L;
			% node0.t = 0;
			node0 = BeliefNode(L,0);
			T_max = node0.find_longest_horizon();

			c_level = [node0];
			BG.N = [node0];

			for tau = 1:T_max
				%Each belief will be indexed by a time. (i.e. I hold X belieft at time t)
				for node_idx = 1:length(c_level) %Iterate through all nodes that are stored in the c_level array
					%Current node
					c_node = c_level(node_idx);

					%Calculate the ancestors of this Belief Node
					temp_post = c_node.post(BG.lcsas,P_u,P_x0);
					
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
					leaf_node_idxs = leaf_node_idxs(leaf_node_idxs ~= cand_leaf_idx);
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

		function [new_lang] = prepend_any_valid_node(varargin)
			%Description:
			%	Takes every word in the input language and 'extends' it backwards by prepending to it
			%	a feasible prefix of a single symbol. If no prefix is possible, the word is left unchanged. 
			%
			%Usage:
			%	new_L = prepend_any_valid_node(obj,word_in)
			%	new_L = prepend_any_valid_node(obj,L_in)
			%Outputs:
			%	new_lang 	- A cell array which contains 0 elements if there are no possible completions. Otherwise, there
			%				  should be a cell array of numeric arrays in new_lang.

			%% Input Processing %%
			obj = varargin{1};

			new_lang = Language();

			if isa(varargin{2},'Language')
				%We want to perform this operation on an entire language.
				in_lang = varargin{2};

				for word_idx = 1:length(in_lang.words)
					%For each word, perform the one step extension.
					temp_word = in_lang.words{word_idx};
					all_new_sublanguages(word_idx) = obj.prepend_any_valid_node(temp_word);
				end
				new_lang = new_lang.union(all_new_sublanguages);

			elseif isnumeric(varargin{2})
				%Hopefully this means we just have a word.
				word_in = varargin{2};

				temp_pre = obj.pre(word_in(1));
				if ~isempty(temp_pre)
					%Create a new word for each node in the pre's findings
					for pre_idx = 1:length(temp_pre)
						new_lang.words{pre_idx} = [temp_pre(pre_idx),word_in];
					end
				else
					new_lang = Language(word_in)
				end
			end

		end

		function [belief_lang,leaf_node_idxs] = get_belief_language(obj)
			%Description:
			%	Looks at the leaf nodes of the graph (nodes that do not have any descendants) and does
			%	backtracking to create a language of beliefs.

			%% Initialize variables

			belief_lang = [];

			%% Find all leaves
			leaf_node_idxs = obj.get_leaf_node_idxs();

			%% Perform backtracking operation for every node 
			belief_lang = Language();
			for leaf_idx = leaf_node_idxs
				sub_lang = Language([leaf_idx]);
				while(~sub_lang.all_words_start_with_root())
					%Until all of the words in this sub_language have reached the root node $n$ 
					temp_lang = obj.prepend_any_valid_node(sub_lang);
					sub_lang = temp_lang;
				end
				belief_lang = belief_lang.union(sub_lang);
			end
		end

	end


end