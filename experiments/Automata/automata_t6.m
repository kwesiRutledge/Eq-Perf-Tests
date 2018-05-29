function [data_out] = automata_t6(varargin)
	%Description:
	%	The objective of this test is to prototype the algorithm that Prof. Ozay has mentioned, using the FSM form which I
	%	believe can better suit this example.

	disp('The objective of this test is to implement plotting functions for these automata using MATLAB''s digraph() function.')

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

	figure;
	plot(G3)
	axis off

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Observable Automaton %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	oa3 = OBSV_AUT(fsm3);
	G3_oa = oa3.gen_digraph();

	x = [0;0;0;0;-1];
	y = [4;3;2;1;1];


	figure;
	plot(G3_oa,'XData',x,'YData',y)
	axis off

	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	data_out.fsm3 = fsm3;
	data_out.oa3 = oa3;
	
end 