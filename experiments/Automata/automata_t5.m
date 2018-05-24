
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
			 	5,1];

	fsm3 = FSM0(X2,X0_2,Y2,H2,Delta2);

	X4 = [1:15]';
	X4_0 = 5;
	Y4 = [0;1];
	H4 = [ 	1,1;
			2,0;
			3,1;
			4,0;
			5,1;
			6,0;
			7,0;
			8,1;
			9,1;
			10,0;
			11,0;
			12,0;
			13,0;
			14,1;
			15,1 ];
	Delta4 = [ 	1,2;
				2,3;
				3,4;
				4,5;
				5,1;
				5,6;
				6,7;
				7,8;
				8,9;
				9,10;
				10,5;
				5,11;
				11,12;
				12,13;
				13,5;
				13,14;
				14,15;
				15,13 ];

	fsm4 = FSM0(X4,X4_0,Y4,H4,Delta4);

    fsm5 = FSM0(X4,X4_0,Y4,H4,[Delta4; 5,5]);
    
    X6 = [1:7]';
    X6_0 = 2;
    Y6 = [0;1];
    H6 = [  1,0;
            2,1;
            3,0;
            4,0;
            5,0;
            6,1;
            7,0];
    Delta6 = [  1,2;
                2,1;
                2,2;
                2,3;
                2,4;
                4,5;
                5,6;
                6,7;
                7,2;
                3,6 ];
            
    fsm6 = FSM0(X6,X6_0,Y6,H6,Delta6);
    
    fsm7 = FSM0(X6,X6_0,Y6,H6,[Delta6;6,3;3,2]);
    
	disp('=================')
	disp('Automata created.')
	disp(' ')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Observable Automaton %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	oa1 = OBSV_AUT(fsm1);
	oa3 = OBSV_AUT(fsm3);
	oa4 = OBSV_AUT(fsm4);
    oa5 = OBSV_AUT(fsm5);
    oa6 = OBSV_AUT(fsm6);
    oa7 = OBSV_AUT(fsm7);
	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	data_out.fsm1 = fsm1;
	data_out.fsm3 = fsm3;
	data_out.fsm4 = fsm4;
    data_out.fsm5 = fsm5;
    data_out.fsm6 = fsm6;
    data_out.fsm7 = fsm7;

	data_out.oa1 = oa1;
	data_out.oa3 = oa3;
	data_out.oa4 = oa4;
    data_out.oa5 = oa5;
    data_out.oa6 = oa6;
    data_out.oa7 = oa7;
end 