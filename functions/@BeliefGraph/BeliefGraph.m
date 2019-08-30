classdef BeliefGraph
	%Description:
	%
	properties
		E;
		N;

	end
	methods
		function [BT] = BeliefGraph()
			disp('Created empty Belief Graph. Please call the construct() function next.')
		end

		function [BN_idx] = find_node_idx(obj,BN)
			%Description:
			%	Identifies the index of the belief node in the set N and returns that index.

			for node_idx = 1:length(obj.N)
				if BN.is_eq(obj.N(node_idx))
					BN_idx = node_idx;
				end 
			end

		end

		function [] = construct(obj,ad_arr,L,P_u,P_x0)

			%%
			warning('The BeliefGraph() class was designed with the assumption that the affine dynamics contain the same state, input, and disturbance size throughout the entire language.')

			%Constants

			%Create first node
			% node0.subset = L;
			% node0.t = 0;
			node0 = BeliefNode(L,0);

			node1.subset = [1:length(L2)];
			node1.t = 0;

			c_level = [node0];

			Phi0 = P_x0*P_u;

			nodes0 = [node0];
			nodes1 = [node1];
			edges0 = [];
			edges1 = []; %Numerical version of the edges matrix

			for tau = 1:T2
				%Each belief will be indexed by a time. (i.e. I hold X belieft at time t)
				for node_idx = 1:length(c_level) %Iterate through all nodes that are stored in the c_level array
					%Current node
					c_node = c_level(node_idx);

					%Calculate the ancestors of this Belief Node
					temp_post = c_node.post(ad_arr,P_u,P_x0);
					
					% Create edges and the next level of the tree
					for edge_ind = 1:length(visible_transitions)
						temp_node.subset = {c_level(node_ind).subset{node_p_set{visible_transitions(edge_ind)}}}; 
						temp_node.t = tau;
						
						temp_node1.subset = [];
						for subset_ind = 1:length(temp_node.subset)
							for L2_ind = 1:length(L2)
								if all(L2{L2_ind} == temp_node.subset{subset_ind})
									temp_node1.subset = [temp_node1.subset, L2_ind];
								end
							end
						end
						temp_node1.t = tau;
						%Add temporary node to the nodes list
						node_in_set_already = false;
						for clevel_ind = 1:length(nodes0)
							% disp(['language_eq(nodes0(clevel_ind).subset,temp_node.subset) = ' num2str(language_eq(nodes0(clevel_ind).subset,temp_node.subset)) ])
							if (language_eq(nodes0(clevel_ind).subset,temp_node.subset)) && (nodes0(clevel_ind).t == temp_node.t)
								node_in_set_already = true;
							end
						end
						if ~node_in_set_already
							nodes0 = [nodes0,temp_node];
							nodes1 = [nodes1,temp_node1];
						end
						%Create edge using this new node and add it to the edges list
						temp_edge = [c_level(node_ind),temp_node];
						edges0 = [edges0;temp_edge];
					end

				end

				%Create next level of the tree
				c_level = [];
				for node_ind = 1:length(nodes0)
					if nodes0(node_ind).t == tau
						c_level = [c_level,nodes0(node_ind)];
					end
				end

				disp(['There are ' num2str(length(c_level)) ' nodes at time tau = ' num2str(tau) '.' ])
			end

			%Clean Up Tree

			for vert_ind = 1:length(nodes0)
				temp_node = nodes0(vert_ind);
				for edge_ind = 1:size(edges0,1)
					%Place a vertex number instead of a whole vertex into the node location in the edge matrix for a more compact representation
					if language_eq(edges0(edge_ind,1).subset,temp_node.subset) && (edges0(edge_ind,1).t == temp_node.t)
						edges1(edge_ind,1) = vert_ind;
					elseif language_eq(edges0(edge_ind,2).subset,temp_node.subset) && (edges0(edge_ind,2).t == temp_node.t)
						edges1(edge_ind,2) = vert_ind;
					end
				end
			end

		BT.E = edges1;
		BT.N = nodes0;	

		end

	end


end