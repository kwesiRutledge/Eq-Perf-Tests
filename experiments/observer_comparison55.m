function [results] = observer_comparison55( varargin )
	%Description:
	%	Tests the new class for Languages.

	L1 = Language([1,2,3],[1:4],[1,2]);
	L2 = Language([1,3,4],[1:4]);
	L3 = Language([1:2],[1:3],[1:4]);
	L4 = Language([1:2],[1:3],[1:4],[1,3,4]);
	L5 = Language([1:2]);
	L6 = Language([1:3]);
	L7 = Language([1:4]);

	results.parameters.L1 = L1;
	results.parameters.L2 = L2;
	results.parameters.L3 = L3;
	results.parameters.L4 = L4;

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

	disp('Experiment 2: Testing Union and Equality')
	disp(' ')

	temp_L = L1.union([L2]);

	disp(['Is the union of L1 and L2 what we expect? ' num2str(temp_L.is_eq(L4)) ])
	temp_L = L5.union([L6,L7]);
	disp(['Is the union of L5,L6,and L7 what we expect? ' num2str(temp_L.is_eq(L1))])
	results.exp2.L_union = temp_L;






end
