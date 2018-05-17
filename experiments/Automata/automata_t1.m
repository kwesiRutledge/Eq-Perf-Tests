function [data_out] = automata_t1(varargin)
	%Description:
	%	The objective of this test is to prototype the functionality developed in the Python Notebook/described in De Santis et. al.'s 
	%	2017 Automatica paper.


	disp('Testing the various helper functions that we need to calculate indistinguishability of Finite State Machines (FSM''s)')
	disp('according to De Santis et. al. in their 2017 Automatica Paper.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Automaton 1, the automaton that was used in the IPython Notebook. %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
	
	disp('Automaton/FSM 1 created.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Testing Output Function %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp(' ')
	disp('Testing Output Function.')
	disp(' ')

	disp([ 'H(X[0]) = ' num2str(fsm1.H_of(X(1))) ])

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Testing Special Matrix Constructors %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp(' ')
	disp('Testing Matrix Constructors')
	disp(' ')

	disp('Pi =')
	disp(num2str(fsm1.create_Pi))

	disp(' ')
	disp('Theta = ')
	disp(num2str(fsm1.create_Theta))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Testing The Pre and Succ, Set-Based Operators %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('=====================================')
	disp('Testing Cart_Prod(), Pre() and Succ()')
	disp(' ')

	disp(['Succ([X0]) = ' num2str(fsm1.succ(fsm1.X0)) ])
	disp(['Succ([X4]) = ' num2str(fsm1.succ(fsm1.X(5)) )])

	disp(['Pre([X0]) = ' num2str(fsm1.pre(fsm1.X0))])
	disp(['Pre([X4]) = ' num2str(fsm1.pre(fsm1.X(5)))])

	disp(['{X0,X2} x {X1,X3} = '])
	disp(num2str(fsm1.cart_prod([X(1),X(3)]',[X(2),X(4)]')))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Testing Backwards Indistinguishability Recursion %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('======================================')
	disp('Testing Backwards Indistinguishability')
	disp(' ')
	disp('A1.Pi \ A1.Theta = ')
	disp(num2str(setdiff(fsm1.Pi,fsm1.Theta,'rows')))
	disp(' ')
	disp(['B_{1}( A1.Pi \ A1.Theta ) = '])
	disp(num2str(setdiff(fsm1.Pi,fsm1.Theta,'rows')))
	disp(' ')
	disp('B_2( A1.Pi \ A1.Theta ) = ')
	disp(num2str( fsm1.B_recur(setdiff(fsm1.Pi,fsm1.Theta,'rows')) ))
	disp(' ')
	disp('B_3( A1.Pi \ A1.Theta ) = ')
	disp(num2str( fsm1.B_recur(fsm1.B_recur(setdiff(fsm1.Pi,fsm1.Theta,'rows'))) ))
	disp(' ')
	disp('B_4( A1.Pi \ A1.Theta ) = ')
	disp(num2str( fsm1.B_recur(fsm1.B_recur(fsm1.B_recur(setdiff(fsm1.Pi,fsm1.Theta,'rows')))) ))

	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	data_out.fsm1 = fsm1;
	data_out.output_tests = [];
	for i = fsm1.X' %'
		data_out.output_tests = [ data_out.output_tests ; fsm1.H_of(i) ];
	end
	data_out.Pi = fsm1.create_Pi;
	data_out.Theta = fsm1.create_Theta;
	data_out.operator_tests.succ_0 = fsm1.succ(fsm1.X0);
	data_out.operator_tests.succ_4 = fsm1.succ(fsm1.X(5));
	data_out.operator_tests.pre_0 = fsm1.pre(fsm1.X0);
	data_out.operator_tests.pre_4 = fsm1.pre(fsm1.X0);
	data_out.operator_tests.cp = fsm1.cart_prod([X(1),X(3)]',[X(2),X(4)]');

	data_out.backwards_tests.B1 = setdiff(fsm1.Pi,fsm1.Theta,'rows');
	data_out.backwards_tests.B2 = fsm1.B_recur(setdiff(fsm1.Pi,fsm1.Theta,'rows'));
	data_out.backwards_tests.B3 = fsm1.B_recur(fsm1.B_recur(setdiff(fsm1.Pi,fsm1.Theta,'rows')));
	data_out.backwards_tests.B4 = fsm1.B_recur(fsm1.B_recur(fsm1.B_recur(setdiff(fsm1.Pi,fsm1.Theta,'rows'))));
end 