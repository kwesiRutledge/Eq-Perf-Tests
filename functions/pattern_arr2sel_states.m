function [ missing_bG ] = pattern_arr2sel_states(single_pattern,n,T)
	%Description:
	%	Converts a single missing data pattern (represented as an arrat)
	%	into a constraint for the design of	Finite Horizon Affine Estimators
	%	that achieve Equalized Recovery.
	%
	%Inputs:
	%	single_pattern - 	A binary arrat of length T that describes whether or not
	%						an observation or data point is missing at a given time. 
	%	n - Dimension of state space.


	%% Constants

	missing_bG = [];

	%% Create matrix
	for p_ind = 0 : T
		temp_row = [zeros(n,p_ind*n) single_pattern(p_ind+1)*eye(n) zeros(n,(T-p_ind)*n)];
		missing_bG = [missing_bG; temp_row];
	end

end