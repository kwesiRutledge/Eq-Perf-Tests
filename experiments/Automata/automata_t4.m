function [data_out] = automata_t4(varargin)
	%Description:
	%	The objective of this test is to prototype the algorithm that Prof. Ozay has mentioned, using the FSM form which I
	%	believe can better suit this example.

	disp('The objective of this test is to implement the ''Observability Automaton'' algorithm that Prof. Ozay has defined.')
	disp('Although it uses similar language to W. Wang et. al.''s 2007 article, I don''t think that this fits their framework.')

	%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Automatons %%
	%%%%%%%%%%%%%%%%%%%%%%%

	X = [0:4]';
	X0 = X(1);
	Y = [0;1];
	H = [0,1;
		 1,0;
		 2,1;
		 3,0;
		 4,1 ];
	Delta = [ 0,1; 1,2; 2,3; 3,4 ; 4,0 ];

	fsm1 = FSM0(X,X0,Y,H,Delta);

	X2   = [0:5]';
	X0_2 = [1];
	Y2   = [0;1];
	H2 = [ 	0,0 ;
			1,1 ;
			2,0 ;
			3,1 ;
			4,0 ;
			5,0 ];
	Delta2 = [	0,1;
			 	1,0;
			 	1,2;
			 	2,3;
			 	3,4;
			 	4,5;
			 	5,0 ];

	fsm3 = FSM0(X2,X0_2,Y2,H2,Delta2);

	disp('=================')
	disp('Automata created.')
	disp(' ')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Recursive Algorithm for FSM1 %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('=================================================')
	disp('Begin Attempted Automaton Construction Algorithm.')
	disp(' ')

	fsm_i = fsm3;

	E_list = {};
	s_node0.xs = fsm_i.X0;
	s_node0.parent = Inf;
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
			%3. Each of those clusters becomes a super node with parent being the current super node
			sn0.xs 		= temp_n0
			sn0.parent 	= Ek(curr_sn_ind).xs;

			sn1.xs 		= temp_n1
			sn1.parent 	= Ek(curr_sn_ind).xs;

			%4. Place these new nodes into the next expansion set IF THEY DO NOT EXIST IN ANY OTHER EXPANSION Ek
			sn0_match = false;
			sn1_match = false;
			for k = 1 : length(E_list)
				%Check every super node in Ek
				for ind0 = 1 : length(E_list{k})
					%Try to match the xs in the list with those in our sn0 or sn1.
					if length(sn0.xs) == length(E_list{k}(ind0).xs)
						if all( unique(sn0.xs) == unique(E_list{k}(ind0).xs) )
							sn0_match = true;
						end
					elseif length(sn1.xs) == length(E_list{k}(ind0).xs)
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
		Ek = Ekp1;
		Ekp1 = [];

		exp_level = exp_level + 1;
		E_list{exp_level} = Ek;

		if exp_level == 10
			break;
		end

	end

	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	data_out.fsm1 = fsm1;
	data_out.fsm3 = fsm3;

	data_out.E_list = E_list;
	data_out.final_exp_level = exp_level;
end 