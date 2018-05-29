function [data_out] = automata_t7(varargin)
	%Description:
	%	The objective of this test is to test some necessary helper functions for the OBSV_AUT class.
	%	These helper functions will create the paths that we need to use to develop prefix-based feedback
	%	that satisfies the specification when an Automaton is the constraint.

	disp('Testing new helper functions for the OBSV_AUT class.')

	%%%%%%%%%%%%%%%%%%%%%%%
	%% Define Automatons %%
	%%%%%%%%%%%%%%%%%%%%%%%

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
			 	5,1 ];
    
	fsm3 = FSM0(X2,X0_2,Y2,H2,Delta2);

	disp('=================')
	disp('Automata created.')
	disp(' ')

	%%%%%%%%%%%%%%%%%%%%
	%% Plot Automaton %%
	%%%%%%%%%%%%%%%%%%%%

	disp('===================')
	disp('Plotting Automaton.')
	disp(' ')

	G3 = fsm3.gen_digraph();

	% figure;
	% plot(G3)
	% axis off

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Observable Automaton %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	oa3 = OBSV_AUT(fsm3);

	disp('=============================')
	disp('Created Observable Automaton?')
	disp(' ')

	G3_oa = oa3.gen_digraph();

	x = [0;0;0;0;-1];
	y = [4;3;2;1;1];

	figure;
	plot(G3_oa,'XData',x,'YData',y)
	axis off

	disp('===================================')
	disp('Testing Pre() and Succ() functions.')
	disp(' ')

	disp(['oa3.pre(1) = ' num2str(oa3.pre(1)') ])
	disp(['oa3.succ(3) = ' num2str(oa3.succ(3)') ' / ' num2str(oa3.OX{ oa3.succ(3) }') ]) %'
	disp(['oa3.succ(4) = ' num2str(oa3.succ(4)') ]) %'

	disp('=====================================================================')
	disp('Searching for all paths that go from a d-set containing recovery node')
	disp('(in this case the OX_0 node) to any other d-set containing recovery')
	disp('node.')
	disp(' ')

	%Suggested Steps:
	%	1. Collect All d-set that contain recovery node.
	%	2. For each d-set,
	%	2a. Create initial expansion level, E0, using just d-set
	%	2b. Calculate the successors of every sequence node in current Expansion Level
	%	2c. Create sequence nodes for each successor
	%	2d. Add sequence node to the next expansion level, E_kp1, if:
	%			- they are novel
	%			- they are not a recovery to recovery path.

	%1. Collect all d-sets that contain recovery node.
	rec_dsets = [];
	for OX_ind = 1 : length(oa3.OX)
		if any(oa3.OX{OX_ind} == oa3.OX_0)
			rec_dsets = [rec_dsets; OX_ind];
		end
	end

	%2. For each d-set
	for dset0 = rec_dsets'
		%Define E_0 = {dset0}
		seq_node0.dset_ind = dset0;
		seq_node0.history = [];
		%2a. Create initial expansion level, E0, using just d-set
		Ek = seq_node0;
		E_list{dset_ind} = {};
		E_list{dset_ind}{1} = Ek;

		exp_level = 1;

		while ~isempty(Ek)
			%2a. For every sequence node in the current expansion level
			for seq_node_ind = 1:length(Ek)

				parent_seq_node = Ek(seq_node_ind);

				%2b. Calculate the successors of every sequence node in current Expansion Level
				succ_k = oa3.succ(parent_seq_node.dset_ind);

				%2c. Create sequence nodes for each successor
				%	 Create 2 new sequence nodes (nodes that keep track of the sequence that led to the current node).
				%	 (There can only be a maximum of 2 new nodes per parent.)
				seq_node0 = [];
				seq_node1 = [];

				for succ_ind = succ_k
					if oa3.H(succ_ind,2) == 0
						seq_node0.dset_ind = succ_ind;
					else
						seq_node1.dset_ind = succ_ind;
					end
				end

				%Add the histories.
				if ~isempty(seq_node0)
					seq_node0.history = [ parent_seq_node.history ; parent_seq_node.dset_ind ];
				end
				if ~isempty(seq_node1)
					seq_node1.history = [ parent_seq_node.history ; parent_seq_node.dset_ind ];
				end

				%2d. Add sequence node to the next expansion level, E_kp1, if:
				%- Newest dset_ind is novel, but not a recovery d-set's index


				E_list{dset_ind}{exp_level+1} = E_kp1;

			end
				 	 
		end
	end


	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	data_out.fsm3 = fsm3;
	data_out.oa3 = oa3;

	data_out.pre_1 = oa3.pre(1);
	data_out.succ_3 = oa3.succ(3);
	data_out.succ_4 = oa3.succ(4);

	data_out.rec_dsets = rec_dsets;
	
end 