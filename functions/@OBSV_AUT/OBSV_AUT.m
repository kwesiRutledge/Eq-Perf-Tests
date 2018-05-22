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
								sn1.xs
								E_list{k}(ind0).xs
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

				%Update Ek and Ekp1;
				E_list{end} = Ek;
				Ek = Ekp1;
				Ekp1 = [];

				exp_level = exp_level + 1;
				E_list{exp_level} = Ek;

				if exp_level == 10
					break;
				end

			end

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
							if isempty( findrows(oa.Delta,temp_edge) )
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

		%%%%%%%%%%%%%%%%%%%%%%
		%% HELPER FUNCTIONS %%
		%%%%%%%%%%%%%%%%%%%%%%

		

	end

end