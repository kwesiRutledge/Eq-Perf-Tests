function [dual_vars, constr] = create_sadraddini_AH_inclusion_constr( obj , x_bar , X , H_x , h_x , y_bar , Y , H_y , h_y )
	%create_sadraddini_cond.m
	%Description:
	%	This function creates the AH-polytope includion in AH-Polytope constraint from
	%	Sadraddini et al. 2019.
	%
	%	The constraint is a sufficient condition for: X \subseteq Y
	%	where 	X = x_bar + X * P_x
	%			Y = y_bar + Y * P_y
	%			P_x = { x \in R^{n_x} | H_x x  <= h_x }
	%			P_y = { y \in R^{n_y} | H_y y  <= h_y }
	%	
	%Usage:
	%	[dual_vars, constr] = cg.create_sadraddini_AH_inclusion_constr( x_bar , X , H_x , h_x , y_bar , Y , H_y , h_y )

	%% Constants

	[q_y,n_y] = size(H_y);
	[q_x,n_x] = size(H_x);

	%% Create Dual Variables

	Gamma0 	= sdpvar(n_y,n_x,'full');
	beta0 	= sdpvar(n_y,1,'full');
	Lambda0 = sdpvar(q_y,q_x,'full');

	dual_vars = {Gamma0,beta0,Lambda0};

	%% Create Constraint

	constr = [];

	constr = constr + [Lambda0 >= 0];

	constr = constr + [X == Y*Gamma0] + [y_bar - x_bar == Y * beta0];
	constr = constr + [Lambda0*H_x == H_y * Gamma0] + [Lambda0*h_x <= h_y + H_y*beta0];

end