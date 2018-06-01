classdef OBSV_AUT
	%Description:
	%	This class is an observable automaton based on the distinguishable states of an automaton created as an FSM0 object,

	%Member Variables:
	%	OX - 	A set containing subsets of the original FSM state space. I may refer to this as a set of indistinguishable state sets or a set of i-sets.
	%			Data Type: |OX| x 1, Cell Array. Each cell is a list of integers (of arbitrary length)
	%
	%	OX_0 - 	The set of potential initial states of the Observable Automaton. The initial i-set.
	%			Data Type: |OX_0| x 1, Array of integers (indices of i-sets in OX).
	%
	%	Y - 	The set of possible outputs for the given i-set.
	%			Data Type: |Y| x 1, Array of integers.
	%
	%	H -		Relation between given states and output values.
	%			Data Type: |OX| x 2 Matrix. First entry of each row is the index of the i-set in OX, second entry is output from Y.
	%
	%	Delta - Transition Relation Between States
	%			Data Type: |Delta| x 2 matrix. First entry is originating i-set's index in OX (from OX), second entry is target state (from X).
	%

	properties
		OX;
		OX_0;
		Y;
		H;
		Delta;
	end

	methods
		%Constructor
		function oa = OBSV_AUT( fsm_i )
			%Description:
			%
			%Process:
			%	1. Use fsm synthesis algorithm to create all extension levels.
			%	2. Use E-levels to create nodes in OX.
			%	3. Use X0 to create initial node OX_0 from OX.
			%	4. Use E-Levels to create edges in OX.
			%
			%
			%Notes:
			%	Input is an automaton in the format of an fsm.
			%

			E_list = {};
			s_node0.xs = fsm_i.X0;
			s_node0.parent = Inf;
			s_node0.next = {};
			Ek = s_node0;

			E_list{1} = Ek;
			exp_level = 1;
			Ekp1 = [];

			while ~isempty(Ek)

				for curr_sn_ind = 1:length(Ek) %'
					%0. Clear Old Variables 

					%1. Evaluate the the succ() operator for the current super node.
					temp_succ = fsm_i.succ( Ek(curr_sn_ind).xs );
					%2. Collect all nodes in the successor set of certain labels.
					temp_n1 = [];
					temp_n0 = [];
					for node_num = temp_succ'
						if fsm_i.H_of(node_num) == 1
							temp_n1 = [temp_n1; node_num];
						elseif fsm_i.H_of(node_num) == 0
							temp_n0 = [temp_n0; node_num];
						else
							error('Unrecognized output label.');
						end	
					end
					%Save the list as a set.
					temp_n0 = unique(temp_n0);
					temp_n1 = unique(temp_n1);

					%3. Each of those clusters becomes a super node with parent being the current super node
					sn0.xs 		= temp_n0;
					sn0.parent 	= Ek(curr_sn_ind).xs;
					sn0.next 	= {};

					sn1.xs 		= temp_n1;
					sn1.parent 	= Ek(curr_sn_ind).xs;
					sn1.next 	= {};

					%3.1 Add each of these nodes to the 'next' field
					Ek(curr_sn_ind).next = { temp_n0 , temp_n1 };

					%4. Place these new nodes into the next expansion set IF THEY DO NOT EXIST IN ANY OTHER EXPANSION Ek
					sn0_match = false;
					sn1_match = false;
					for k = 1 : length(E_list)
						%Check every super node in Ek
						for ind0 = 1 : length(E_list{k})
							%Try to match the xs in the list with those in our sn0 or sn1.
							if length(unique(sn0.xs)) == length(E_list{k}(ind0).xs)
								if all( unique(sn0.xs) == unique(E_list{k}(ind0).xs) )
									sn0_match = true;
								end
							elseif length(unique(sn1.xs)) == length(E_list{k}(ind0).xs)
								if all( unique(sn1.xs) == unique(E_list{k}(ind0).xs) )
									sn1_match = true;
								end
							end
						end
					end

					if (~sn1_match) & ~isempty(sn1.xs)
						%If sn1 did not match, then add it to the next level to expand.
						Ekp1 = [Ekp1 sn1];
					end

					if (~sn0_match) & ~isempty(sn0.xs)
						%IF sn0 did not match, then add it to the next level to expand.
						Ekp1 = [Ekp1 sn0];
					end	
				end	 

				%Save things to E_list
				E_list{exp_level} = Ek;
				exp_level = exp_level + 1;
				E_list{exp_level} = Ekp1;

				%Update Ek and Ekp1;
				Ek = Ekp1;
				Ekp1 = [];

				if exp_level == 10
					break;
				end

			end
			% E_list
			% E_list{3}
			% E_list{3}(1)

			%2. Use E-levels to create nodes in OX.
			node_ind = 1;
			for k = 1:length(E_list)-1
				%Investigate every expansion level to get the unique nodes contained in it.
				for entry_ind = 1 : length(E_list{k})
					oa.OX{node_ind} = E_list{k}(entry_ind).xs;
					node_ind = node_ind + 1;
				end
			end

			%3. Use X0 to create initial node INDEX OX_0 from OX.
			for i_set_ind = 1 : length(oa.OX)
				%If there is a state that matches the initial state from the FSM, make that the initial i-set of the
				%observable automaton.
				if oa.OX{i_set_ind} == fsm_i.X0
					oa.OX_0 = i_set_ind;
				end
			end

			%4. Use E-Levels to create edges in OX.
			oa.Delta = [];
			% oa.Delta2 = [];
			for k = 1:length(E_list)-1
				%For every i-set in the expansion level, connect it to its parent (if the parent exists).
				for entry_ind = 1 : length(E_list{k})
					% %4.1 Check to see if parent exists
					% if E_list{k}(entry_ind).parent ~= Inf
					% 	%Find parent, now that we know it exists.
					% 	parent_ind = oa.get_OX_ind_of( E_list{k}(entry_ind).parent );
					% 	curr_ind = oa.get_OX_ind_of( E_list{k}(entry_ind).xs );
						
					% 	oa.Delta = [oa.Delta; parent_ind, curr_ind];
					% end
					%4.2 Check the value of the 'next' field
					for next_ind = 1:length(E_list{k}(entry_ind).next)
						%Check to see if the entry exists
						if ~isempty(E_list{k}(entry_ind).next{next_ind})
							%If the entry exists, then select that edge.
							temp_edge = [ oa.get_OX_ind_of( E_list{k}(entry_ind).xs ) , oa.get_OX_ind_of( E_list{k}(entry_ind).next{next_ind} ) ];
							if isempty(oa.Delta)
								oa.Delta = [temp_edge];
							elseif ~any( ismember(oa.Delta,temp_edge,'rows') )
								oa.Delta = [oa.Delta ; temp_edge];
							end
						end
					end
				end
			end

			%5. Use the label function from the FSM0 to make a labeling function for these nodes
			oa.H = [];
			for OX_ind = 1 : length(oa.OX)
				%Check each i-set for the label
				oa.H = [ oa.H ; OX_ind , fsm_i.H_of( oa.OX{OX_ind}(1) ) ];
			end

			%6. Set of outputs are identical to that of the FSM
			oa.Y = fsm_i.Y;

		end

		function [ox_ind] = get_OX_ind_of(obj,OX_val)
			%Description:
			%	Searches through the OX cell array to find the proper index that matches
			%	the value OX_val.
			%
			%Assumptions:
			%	Assumes OX_val is in OX
			ox_ind = Inf;
			for ind = 1 : length(obj.OX)
				%First check to see if lengths match.
				if(length(obj.OX{ind}) == length(OX_val))
					%Then check to see if the values are identical.
					if all(OX_val == obj.OX{ind})
						ox_ind = ind;
					end
				end
			end

			if isinf(ox_ind)
				error(['OX_val ' num2str(OX_val') ' does not correspond to anything in OX.'])
			end
		end

		function [pre_set] = pre(obj,ox_in)
			%Description:
			%	Calculates the set of states pre_set that could have preceded ANY d-set in ox_in.
			%
			%Output:
			%	pre_set - A vector of indices. Each index corresponds to the index of a d-set in OX.
			%
			
			pre_set = [];
			for transition_num = 1 : size(obj.Delta,1)
				%Search through every transition in the transition relation Delta
				if any(obj.Delta(transition_num,2) == ox_in )
					pre_set = [ pre_set ; obj.Delta(transition_num,1) ];
				end
			end
			pre_set = unique(pre_set);
		end

		function [succ_set] = succ(obj,ox_in)
			%Description:
			%	Calculates the set of states pre_set that could have preceded ANY d-set in ox_in.
			%
			%Output:
			%	pre_set - A vector of indices. Each index corresponds to the index of a d-set in OX.
			%
			
			succ_set = [];
			for transition_num = 1 : size(obj.Delta,1)
				%Search through every transition in the transition relation Delta
				if any(obj.Delta(transition_num,1) == ox_in )
					succ_set = [ succ_set ; obj.Delta(transition_num,2) ];
				end
			end
			succ_set = unique(succ_set);
		end

		function [q] = H_of(obj,ox_in)
			%Description:
			%	Finds the appropriate labels for every d-set INDEX that is provided in the ox_in vector.
			%
			%Inputs:
			%	ox_in:	A column vector of indices (thus positive numbers) referencing existing d-sets in OX

			q = [];
			for ind = ox_in'
				q = [q; obj.H( obj.H(:,1) == ind , 2 ) ];
			end
		end
		%%%%%%%%%%%%%%%%%%%%%%
		%% HELPER FUNCTIONS %%
		%%%%%%%%%%%%%%%%%%%%%%

		function [L,S] = find_all_2R_paths(obj)
			%Description:
			%	Find all feasible paths that start and end with i_sets containing OX_O.
			%	In other words, finds all paths in the graph such that:
			%		- the d-set at the beginning contains the index OX_0
			%		- the d-set at the end contains the index OX_0
			%		- there are no other d-sets in the path containing index OX_0

			%Create Heap of All Paths 
			L = {}; S = {};
			V = [];

			%Search starting from the initial node.
			node0.state_seq = [ obj.OX_0 ];
			node0.ind = obj.OX_0;

			E_0 = [node0];
			E_list{1} = E_0;
			E_k = E_0;
			V = [obj.OX_0];

			while(~isempty(E_k) )
				%Initialization.
				E_kp1 = [];

				%For every i-set in the current expansion level
				for node_i_ind = 1:length(E_k)
					node_i = E_k(node_i_ind);
					%1. Find successors
					succ_i = obj.succ( node_i.ind );
					%2. For each successor,
					for s_ind = succ_i'
						%2.2 If the successor is a recovery node, then add the current state sequence
						%	 to list S.
						if any( obj.OX_0 == obj.OX{s_ind} )
							S{length(S)+1} = [ node_i.state_seq ; s_ind ];

							node_ip1.ind = s_ind;
							node_ip1.state_seq = [ s_ind ];

							%add to the next expansion level if the node doesn't already exist in a previous level.
							if ~any(V == s_ind)
								E_kp1 = [E_kp1 ; node_ip1 ];
								V = [V; s_ind];
							end
						%2.3 If the successor is a standard node, then we have found a loop.
						elseif any( s_ind == node_i.state_seq )
								disp('error. Off recovery loop detected.')
						%2.1 If the successor is not in the current word label then add the new successor node to next expansion level.
						elseif ~any( s_ind == node_i.state_seq )
							%Create node
							node_ip1.ind = s_ind;
							node_ip1.state_seq = [ node_i.state_seq ; s_ind ];

							%Add node to expansion level.
							E_kp1 = [E_kp1 ; node_ip1];

							V = [V; s_ind];
							
						end
					end
				end
				%Update expansion level.
				E_list{length(E_list)+1} = E_kp1;
				E_k = E_kp1

			end


			for word_ind = 1:length(S)
				L{word_ind} = obj.H_of(S{word_ind});
			end
		end

		function [G] = gen_digraph(obj)
			%Description:
			%	Plots the automaton using matlab's inherent graph tools on the current figure.
			%

			%1. Define Adjacency Matrix
			am = zeros(length(obj.OX));
			for edge_num = 1:size(obj.Delta,1)
				%Index by the edge number and change a few entries to 1.
				am( obj.Delta(edge_num,1) , obj.Delta(edge_num,2) ) = 1;
			end
			%2. Define Node Names
			names = {};
			obj.OX
			for name_ind = 1 : length(obj.OX)
				obj.OX{name_ind}
				names{name_ind} = num2str(obj.OX{name_ind}'); %'
			end

			G = digraph(am,names);	
		end

	end

end