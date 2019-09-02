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

	end
	methods
		function [BG] = BeliefGraph(ad_arr,L,P_u,P_x0)
			% disp('Created empty Belief Graph. Please call the construct() function next.')
			% BT.E = [];
			% BT.N = [];

			%%
			warning('The BeliefGraph() class was designed with the assumption that the affine dynamics contain the same state, input, and disturbance size throughout the entire language.')

			%Constants

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
					temp_post = c_node.post(ad_arr,P_u,P_x0);
					
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

		function [BG] = construct(obj,ad_arr,L,P_u,P_x0)

			%%
			warning('The BeliefGraph() class was designed with the assumption that the affine dynamics contain the same state, input, and disturbance size throughout the entire language.')

			%Constants

			BG = BeliefGraph();

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
					temp_post = c_node.post(ad_arr,P_u,P_x0);
					
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

			%Clean Up Tree

			% for vert_ind = 1:length(nodes0)
			% 	temp_node = nodes0(vert_ind);
			% 	for edge_ind = 1:size(edges0,1)
			% 		%Place a vertex number instead of a whole vertex into the node location in the edge matrix for a more compact representation
			% 		if language_eq(edges0(edge_ind,1).subset,temp_node.subset) && (edges0(edge_ind,1).t == temp_node.t)
			% 			edges1(edge_ind,1) = vert_ind;
			% 		elseif language_eq(edges0(edge_ind,2).subset,temp_node.subset) && (edges0(edge_ind,2).t == temp_node.t)
			% 			edges1(edge_ind,2) = vert_ind;
			% 		end
			% 	end
			% end	

		end

	end


end