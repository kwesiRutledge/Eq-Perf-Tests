function [data_out] = automata_t3(varargin)
	%Description:
	%	The objective of this test is to prototype the functionality developed in the Python Notebook/described in De Santis et. al.'s 
	%	2017 Automatica paper.
	%	Uses a different FSM than test # 1.

	disp('Attempting to Gain insights on what an observability automaton will do for the 2 examples that were so difficult before.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create an Example Automaton %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Creating Example 1 which should also be in Test 1
	X = [0:4]';
	E = [0,1]';
	f = [ 0,1,1;
		  1,0,2;
		  2,1,3;
		  3,0,4;
		  4,1,0];
	Gamma = [ 0,1 ;
			  1,0 ;
			  2,1 ;
			  3,0 ;
			  4,1 ];
	x0 = 0;

	fsa0 = FSA(X,E,f,Gamma,x0);
	disp('===============================')
	disp('Finite State Automaton created.')
	disp(' ')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define a custom index function %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('======================================')
	disp('Also defining a custom index function.')
	disp('Following the paper, this will be in terms of observable and unobservable words, but this is flexible.')
	
	E_o = [1];
	TR_G = fsa0.TR;
	cust_I = [];
	for tr_ind = 1 : size(TR_G,1)
		cust_I = [cust_I; TR_G(tr_ind,:), TR_G(tr_ind,2) ];
	end
	disp('cust_I = ')
	
	% fsa0.def_script_I(cust_I);
	fsa0.ind_func = cust_I;
	disp(num2str(fsa0.ind_func))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Information Mapping/Projection %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('===============================')
	disp('Testing the Projection Function')
	disp('Will require recursion.')
	disp(' ')
	disp('Test eval_I function.')

	xt = 1;
	et = 0;

	disp(['I(' num2str(xt) ',' num2str(et) ') = ' num2str(fsa0.eval_I(xt,et)) ])

	disp(' ')
	disp(['theta([1,0,1,0,1]) = ' num2str(fsa0.info_map([1,0,1,0,1]))  ]) 
	disp(['theta([1,0,1,0,1,1]) = ' num2str(fsa0.info_map([1,0,1,0,1,1])) ])
	disp(['theta([1,0,1,0]) = ' num2str(fsa0.info_map([1,0,1,0]))  ])  


	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	disp(' ')
	disp('Saving Results...')

	data_out.fsa = fsa0;
	data_out.E_o = E_o;
	data_out.cust_I = cust_I;
end 