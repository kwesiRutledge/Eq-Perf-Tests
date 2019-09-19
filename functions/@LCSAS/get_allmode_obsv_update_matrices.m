function [ A_Lt , B_Lt , k_Lt , L_Lt , C_Lt ] = get_allmode_obsv_update_matrices(obj,t,obsv_gains)
	%Description:
	%	This function constructs the matrices that are used in the "All Mode" Observer to calculate the
	%	the next state of all modes given the current input.
	%
	%Inputs:
	%	t 			- The current time of the LCSAS.
	%	obsv_gains 	- An n_x x n_y x n_modes array that stores the observer gains for each mode observer.
	%
	%Usage:
	%	[ A_Lt , B_Lt , k_Lt , L_Lt , C_Lt ] = lcsas.get_allmode_obsv_matrices(t,obsv_gains)

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%Check the value of t
	if (t < 0) || (t > obj.L.find_longest_length()-1 )
		error(['t must be between 0 and ' num2str(obj.L.find_longest_length()) ' (the longest length of a word in L).' ] )
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(obj.Dyn(1).A,1);
	m = size(obj.Dyn(1).B,2);
	p = size(obj.Dyn(1).C,1);
	wd = size(obj.Dyn(1).B_w,2); %Assume that the size of the disturbance doesn't change for any of the included dynamics
	vd = size(obj.Dyn(1).C_v,2);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	A_Lt = []; B_Lt = []; C_Lt = [];
	L_Lt = []; k_Lt = [];

	for lang_idx = 1:length(obj.L.words)
		%Get Current Word and Mode at t
		temp_word = obj.L.words{lang_idx};
		temp_mode = temp_word(t+1);

		%Edit the desired matrices
		A_Lt(end+[1:n],end+[1:n]) = obj.Dyn(temp_mode).A;
		B_Lt(end+[1:n],[1:m]) = obj.Dyn(temp_mode).B;
		k_Lt(end+[1:n],1) = obj.Dyn(temp_mode).f;

		L_Lt(end+[1:n],end+[1:p]) = obsv_gains(:,:,temp_mode);
		C_Lt(end+[1:p],end+[1:n]) = obj.Dyn(temp_mode).C;

	end

end