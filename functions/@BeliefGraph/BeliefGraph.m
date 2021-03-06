 classdef BeliefGraph < handle
	%Description:
	%
	%Member Variables (Protected):
	%	UsedProjection 			- Boolean.
	%					  Flag that is used to tell the BeliefGraph constructor to not use projection in the construction.
	%	FeedbackMethod 			- String.
	%					  Options: 'output' , 'state'
	%					  Feedback method used when constructing consistency sets.
	%	UsedAcceleratedAlgorithms 	- Boolean.
	%					  Uses accelerated version of the post algorithm in the construction if true.
	%	UsedUnobservabilityChecks 	- Boolean.
	%					  Does not perform obserability (containment) checks if true.
	%
	%Member Variables:
	%	E		- An n_E x 2 matrix of strictly positive integers, where each row represents an edge in the graph.
	%			  For a row E_i, the first entry (first column) is the INDEX of the edge's source node (i.e. N(E_i(1)) is
	%			  the source node of the edge) and the second entry (second column) is the INDEX of the edge's
	%			  destination node (i.e. N(E_i(2)) is the destination node of the edge).
	%	N		- 
	%	lcsas 		- Language Constrained Switched Affine System. (LCSAS object)
	%	ModeLanguage 	-
	%	BeliefLanguage 	-
	%
	%Member Functions
	%	- find_node_idx
	%	- get_all_nodes_at_time
	%	- plot
	%	- get_leaf_node_idxs
	%	- pre
	%	- prepend_any_valid_node
	%	- get_belief_language

	properties(SetAccess = protected)
		UsedProjection;
		FeedbackMethod;
		UsedAcceleratedAlgorithms;
		UsedUnobservabilityChecks;
		ConsistencySetVersion;
	end

	properties
		E;
		N;
		lcsas;
		ModeLanguage;
		BeliefLanguage;
	end

	methods
		function [BG] = BeliefGraph(varargin)
			%Description:
			%
			%Inputs:
			%	in_sys 	 	- An array of Aff_Dyn() objects. May eventually become its own class/data type soon.
			%	in_lcsas 	- An LCSAS object representing the desired system.
			%	P_u 		- A Polyhedron() object that defines the input constraints.
			%				  The input at each time must lie within the polyhedron.
			%	P_x0 		- A Polyhedron() object that defines the set of states from which the initial state is contained.
			%
			%	fb_method 	- A string indicating if the controller has access to the state or a measurement of the state.
			%				  Options include: 'output', 'state'
			%	post_flag 	- A string indicating if the Belief Graph construction algorithm uses the original version of post,
			%				  or any other.
			%				  Options include: 'minimal','no projection'
			%
			%
			%Usage:
			%	BG = BeliefGraph(in_sys,L,P_u,P_x0)
			%	BG = BeliefGraph(in_sys,L,P_u,P_x0,'verbosity',verbosity)
			%	BG = BeliefGraph(in_lcsas,P_u,P_x0,'verbosity',verbosity)
			%	BG = BeliefGraph(in_lcsas,P_u,P_x0,'verbosity',verbosity,'fb_method',fb_method)
			%	BG = BeliefGraph(in_lcsas,P_u,P_x0,'accel_flag',accel_flag)
			%	BG = BeliefGraph(in_lcsas,P_u,P_x0,'use_proj_flag',tf)
			%	BG = BeliefGraph(in_lcsas,P_u,P_x0,'use_unobservability_checks',tf)
			%	BG = BeliefGraph(in_lcsas,P_u,P_x0,'return_empty',tf)


			% disp('Created empty Belief Graph. Please call the construct() function next.')
			% BT.E = [];
			% BT.N = [];

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			in_sys = varargin{1};
			
			if ~isa(in_sys,'LCSAS')
				if nargin < 4
					error('Not enough inputs to construct BeliefGraph.')
				end

				L = varargin{2};
				P_u = varargin{3};
				P_x0 = varargin{4};

				varargin_idx = 5;
			else
				if nargin < 3
					error('Not enough inputs to construct BeliefGraph.')
				end

				L = in_sys.L;
				P_u = varargin{2};
				P_x0 = varargin{3};

				varargin_idx = 4;
			end


			while(varargin_idx <= nargin)
				switch varargin{varargin_idx}
					case 'verbosity'
						verbosity = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					case 'fb_method'
						fb_method = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					case 'post_flag'
						post_flag = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					case 'accel_flag'
						accel_flag = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					case 'use_proj_flag'
						use_proj_flag = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					case 'use_unobservability_checks'
						use_unobservability_checks = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					case 'return_empty'
						return_empty_flag = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					case 'ConsistencySetVersion'
						ConsistencySetVersion = varargin{varargin_idx+1};
						varargin_idx = varargin_idx + 2;
					otherwise
						error(['Unrecognized input to the function: ' varargin{varargin_idx}])
				end
			end

			if ~isa(L,'Language')
				error('Expected L to be a Language object.')
			end

			if ~isa(in_sys,'LCSAS')
				warning('Expected input system to be LCSAS. Checking to see if this can be salvaged...')
				if isa(in_sys,'Aff_Dyn')
					temp_sys_arr = in_sys;
					in_sys = LCSAS(temp_sys_arr,L);
				else
					error('System must be given either as an LCSAS object or as an array of Aff_Dyn objects.')
				end
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			if ~exist('verbosity')
				verbosity = 0;
			end

			if ~exist('fb_method')
				fb_method = 'output';
			end

			if ~exist('accel_flag')
				accel_flag = false;
			end

			if ~exist('use_proj_flag')
				use_proj_flag = true;
			end

			if ~exist('return_empty_flag')
				return_empty_flag = false;
			end

			if ~exist('use_unobservability_checks')
				use_unobservability_checks = true;
			end

			if ~exist('ConsistencySetVersion')
				ConsistencySetVersion = 1;
			end

			BG.UsedProjection = use_proj_flag;
			BG.FeedbackMethod = fb_method;
			BG.UsedAcceleratedAlgorithms = accel_flag;
			BG.UsedUnobservabilityChecks = use_unobservability_checks;
			BG.ConsistencySetVersion = ConsistencySetVersion;

			%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Building Belief Graph %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%

			if verbosity > 1
				warning('The BeliefGraph() class was designed with the assumption that the affine dynamics contain the same state, input, and disturbance size throughout the entire language.')
			end

			if return_empty_flag
				BG.lcsas = in_sys;
				BG.ModeLanguage = L;
				BG.E = [];
				BG.N = [];
				BG.BeliefLanguage = [];
				return
			end

			if verbosity > 1
				disp('Worked through inputs.')
			end

			in_sys.X0 = P_x0;
			in_sys.U = P_u;

			BG.lcsas = in_sys;
			BG.ModeLanguage = L;

			%Create first node
			% node0.subset = L;
			% node0.t = 0;
			%node0 = BeliefNode(L,0,in_sys.Dyn(1).C*P_x0+in_sys.Dyn(1).C_v*in_sys.Dyn(1).P_v);
			init_nodes = BG.get_initial_beliefnodes( 'X0' , P_x0 );
			temp_node = init_nodes(1);
			T_max = temp_node.find_longest_horizon();

			c_level = init_nodes;
			BG.N = init_nodes;

			for tau = 1:T_max-1
				%Each belief will be indexed by a time. (i.e. I hold X belieft at time t)

				if verbosity > 0
					disp(['- tau = ' num2str(tau) ])
				end

				for node_idx = 1:length(c_level) %Iterate through all nodes that are stored in the c_level array
					%Current node
					c_node = c_level(node_idx);

					%Calculate the ancestors of this Belief Node
					temp_post = BG.post(c_node,'debug',verbosity);
										
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

				%Prepare to create next level of the tree
				c_level = BG.get_all_nodes_at_time(tau); % Use the nodes at the current time (tau) as c_level, in the next iteration.
				if verbosity >= 1
					disp(['There are ' num2str(length(c_level)) ' nodes at time tau = ' num2str(tau) '.' ])
				end
			end

			%Create Belief Language
			BG.BeliefLanguage = BG.get_belief_language();

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

		function [ subN , bn_indices_at_t ] = get_all_nodes_at_time(obj,tau)
			%Description:
			%	Retrieves all belief nodes that have time value provided.
			%
			%Usage:
			%	[subN] = BG.get_all_nodes_at_time(tau)
			%
			%Outputs:
			%	subN - 	An array of BeliefNode objects.

			%% Input Processing %%

			%% Algorithm %%

			subN = []; bn_indices_at_t = [];
			for node_idx = 1:length(obj.N)
				temp_BN = obj.N(node_idx);
				if temp_BN.t == tau
					subN = [subN,temp_BN];
					bn_indices_at_t = [bn_indices_at_t,node_idx];
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

			belief_lang = Language();

			%% Find all leaves
			leaf_node_idxs = obj.get_leaf_node_idxs();

			%% Perform backtracking operation for every node 
			for leaf_idx = leaf_node_idxs
				sub_lang = Language([leaf_idx]);
				while ~obj.all_words_start_with_an_initial_node( sub_lang )
					%Until all of the words in this sub_language have reached the root node $n$ 
					temp_lang = obj.prepend_any_valid_node(sub_lang);
					sub_lang = temp_lang;
				end
				belief_lang = belief_lang.union(sub_lang);
			end
		end

		function tf = all_words_start_with_an_initial_node(obj,Language_in)
			%Description:
			%	Consider the words in the language (Language_in).
			%	Every word represents a sequence of nodes in this belief graph (obj).
			%	Return true if and only if all words begin with an index which is in the set of initial nodes.

			%% Constants %%

			[~,initial_nodes] = obj.get_all_nodes_at_time(0);

			%% Algorithm %%
			tf = true;

			for word_idx = 1:Language_in.cardinality()
				temp_word = Language_in.words{word_idx};
				first_idx = temp_word(1);
				if ~any(first_idx == initial_nodes)
					tf = false;
				end
			end

		end

		function tf = is_covered_by(obj,BG_in)
			%BeliefGraph.is_covered_by()
			%Description:
			%	Observes if all the nodes and edges in the current belief graph are in the input Belief Graph.
			%
			%Usage:
			%	tf = bg0.is_covered_by( bg1 )

			%% Constants %%

			my_nodes = obj.N;
			my_edges = obj.E;

			%% Algorithm %%

			tf = true;

			% Consider every node in current object.	
			for mn_idx = 1:length(my_nodes)
				temp_bn = my_nodes(mn_idx);

				bn_contained = false;
				for bgin_node_idx = 1:length(BG_in.N)
					target_bn = BG_in.N( bgin_node_idx );
					%Identify if target BeliefNode is the same as temp_bn.
					bn_match_test = temp_bn == target_bn;

					bn_contained = bn_contained | bn_match_test;
					
					%If the two nodes match, then exit this loop.
					if bn_match_test
						break;
					end
				end

				tf = tf & bn_contained;
			end

			% Consider every edge object.
			for me_idx = 1:size(my_edges,1)
				temp_edge = my_edges(me_idx,:);
				
				bn_source = my_nodes( temp_edge(1) );
				bn_dest = my_nodes( temp_edge(2) );

				edge_contained = false;
				for bgin_edge_idx = 1:size(BG_in.E,1)
					target_edge = BG_in.E( bgin_edge_idx , : );

					source_and_dest_match_test = ( bn_source == BG_in.N( target_edge(1) ) ) & ( bn_dest == BG_in.N( target_edge(2) ) );
					edge_contained = edge_contained | source_and_dest_match_test;

					%If the two edges match then exit this loop.
					if source_and_dest_match_test
						break;
					end
				end

				tf = tf & edge_contained;

			end

		end

		function tf = is_abstracted_by(obj,BG_in)
			%BeliefGraph.is_covered_by()
			%Description:
			%	Observes if, at each time, there exists a node in BG_in which ABSTRACTS the current BeliefGraph (obj).
			%
			%Usage:
			%	tf = bg0.is_abstracted_by( bg1 )

			%% Constants %%

			my_nodes = obj.N;
			my_edges = obj.E;

			%% Algorithm %%

			tf = true;

			% Consider every node in current object.	
			for mn_idx = 1:length(my_nodes)
				temp_bn = my_nodes(mn_idx);
				temp_t  = temp_bn.t;

				bn_contained = false;
				[ bgin_level_t , ~ ] = BG_in.get_all_nodes_at_time( temp_t );
				for bgin_node_idx = 1:length(bgin_level_t)
					target_bn = bgin_level_t( bgin_node_idx );
					%Identify if target BeliefNode's language contains the language of temp_bn.
					bn_match_test = temp_bn.subL.subseteq( target_bn.subL );

					bn_contained = bn_contained | bn_match_test;
					
					%If the two nodes match, then exit this loop.
					if bn_match_test
						break;
					end
				end

				tf = tf & bn_contained;
			end
			

		end



	end


end