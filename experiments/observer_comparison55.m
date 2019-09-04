function [results] = observer_comparison55( varargin )
	%Description:
	%	Tests the new class for Languages.

	L1 = Language([1,2,3],[1:4],[1,2]);
	L2 = Language([1,3,4],[1:4]);

	results.parameters.L1 = L1;
	results.parameters.L2 = L2;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Experiment 1: Testing Containment %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Experiment 1: Testing Containment')
	disp(' ')

	temp_L = L1.contains([1,2]);

	disp(['L1.contains([1,2]) ? ' num2str(L1.contains([1,2])) ])
	disp(['L2.contains([1,2]) ? ' num2str(L2.contains([1,2])) ])

	results.exp1.contains = temp_L;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Experiment 2: Testing Union %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Experiment 2: Testing Containment')
	disp(' ')

	temp_L = L1.union([L2]);

	results.exp2.L_union = temp_L;




end
