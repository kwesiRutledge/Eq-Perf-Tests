function [l_diag_constr] = create_lower_diagonal_constraint_on_gains( lcsas , Gain_set , fb_type )
	%get_causal_constr_on_extd_gains.m
	%Description:
	%	This function should use the dimensions of the lcsas to create causality constraints
	%	on the transformed feedback gains Q.
	%
	%Usage:
	%	[constraints] = lcsas.create_lower_diagonal_constraint_on_gains( K_set )
	%
	%Inputs:
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();
	T = length(lcsas.L.words{1});

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	l_diag_constr = [];

	switch fb_type
	case 'Disturbance'

		for gain_idx = 1:length(Gain_set)
			%*Strictly* Lower Diagonal Constraint
            l_diag_constr = l_diag_constr + [ Gain_set{gain_idx}([1:n_u],:) == 0 ]; 
			for bl_row_num = 1 : T-1
				l_diag_constr = l_diag_constr + [ Gain_set{gain_idx}(	[bl_row_num*n_u+1:(bl_row_num+1)*n_u], ...
																		[bl_row_num*n_w+1:end] ) == 0 ];
			end
		end

	otherwise
		error(['Unexpected feedback type given to create_lower_diagonal_constraint_on_gain(): ' fb_type ])
	end
	
end