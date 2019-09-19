function [ S_obsv , H_obsv , Cbar_obsv , J_obsv , k_obsv ] = get_allmode_obsv_mpc_matrices(obj,obsv_gains)
	%Description:
	%	Creates the MPC Matrices for the all mode observer.
	%
	%Inputs:
	%	obsv_gains 	- An n_x x n_y x n_modes array that stores the observer gains for each mode observer.
	%
	%Usage:
	%	[ S_obsv , H_obsv , Cbar_obsv , J_obsv , k_obsv ] = lcsas.get_allmode_obsv_mpc_matrices(obsv_gains)

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(obj.Dyn(1).A,1);
	m = size(obj.Dyn(1).B,2);
	p = size(obj.Dyn(1).C,1);
	wd = size(obj.Dyn(1).B_w,2); %Assume that the size of the disturbance doesn't change for any of the included dynamics
	vd = size(obj.Dyn(1).C_v,2);

	Pw_Lts = {}; Pv_Lt = {};
	for t = 0:obj.L.find_longest_length()-1
		temp_Pw_Lt = 1; temp_Pv_Lt = 1;
		for sigma_idx = 1:length(obj.L.words)
			temp_word = obj.L.words{sigma_idx};
			temp_Pw_Lt = temp_Pw_Lt * obj.Dyn(temp_word(t+1)).P_w;
			temp_Pv_Lt = temp_Pv_Lt * obj.Dyn(temp_word(t+1)).P_v;
		end
		%Save in a cell matrix
		Pw_Lts{t+1} = temp_Pw_Lt;
		Pv_Lts{t+1} = temp_Pv_Lt;
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Compute the All Mode Observer Matrices for Each Mode and Store them in a Dyn Object
	temp_dyn_list = [];
	for t = 0:obj.L.find_longest_length()-1
		[ A_Lt , B_Lt , k_Lt , L_Lt , C_Lt ] = obj.get_allmode_obsv_update_matrices(t,obsv_gains)
		temp_dyn = Aff_Dyn(A_Lt,B_Lt,k_Lt, C_Lt,Pw_Lts{t+1},Pv_Lts{t+1});

		%Append new dynamics object to the 
		temp_dyn_list = [temp_dyn_list,temp_dyn]
	end

	%Now create the mpc matrices using the built in get_mpc_matrices() function.
	allmode_sys = LCSAS(temp_dyn_list,Language([1:obj.L.find_longest_length()]));
	[H_obsv,S_obsv,Cbar_obsv,J_obsv,k_obsv] = obj.get_mpc_matrices('word',obj.L.words{1});

end