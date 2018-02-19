function [ select_influenced_states_mat ] = pattern_arr2sel_states(pattern_array,n,T)
	%Description:
	%	Converts a matrix of possible data patterns into a constraint for the design of
	%	Finite Horizon Affine Estimators that achieve Equalized Recovery.
	%
	%Inputs:
	%
	%

	%% Constants

	select_influenced_states_mat = [];

	%% Create matrix
	for pattern_num = 1 : size(pattern_array,1)
		temp_row = [];
		for t = 0 : T
			temp_row = [temp_row pattern_array(pattern_num,t+1)*eye(n)];
		end
		select_influenced_states_mat = [select_influenced_states_mat; temp_row];
	end

end