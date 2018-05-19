function [data_out] = automata_t5(varargin)
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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Observable Automaton %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	oa1 = OBSV_AUT(fsm1);
	oa3 = OBSV_AUT(fsm3);


	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	data_out.fsm1 = fsm1;
	data_out.fsm3 = fsm3;

	data_out.oa1 = oa1;
	data_out.oa3 = oa3;
end 