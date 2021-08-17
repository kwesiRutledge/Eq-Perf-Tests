function [l_diag_constr] = get_causal_constr_on_extd_gains( obj , lcsas , Q_set )
	%get_causal_constr_on_extd_gains.m
	%Description:
	%	This function should use the dimensions of the lcsas to create causality constraints
	%	on the transformed feedback gains Q.
	%
	%Usage:
	%	[constraints] = cg.get_causal_constr_on_extd_gains( lcsas , Q_set )
	%
	%Inputs:
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if ~isa(lcsas,'LCSAS')
		error('The second input must be a LCSAS object.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);
	p = size(lcsas.Dyn(1).C,1);
	wd = size(lcsas.Dyn(1).B_w,2); %Assume that the size of the disturbance doesn't change for any of the included dynamics
	vd = size(lcsas.Dyn(1).C_v,2);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	l_diag_constr = [];

	for gain_idx = 1:length(Q_set)
		T_i = size(Q_set{gain_idx},2)/p;
		%Lower Diagonal Constraint
		for bl_row_num = 1 : T_i-1
			l_diag_constr = l_diag_constr + [ Q_set{gain_idx}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
																[bl_row_num*p+1:end] ) == 0 ];
		end
	end
	
end