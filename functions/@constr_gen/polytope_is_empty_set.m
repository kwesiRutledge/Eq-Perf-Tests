function [ constraint_out , y ] = polytope_is_empty_set( cg , H , h )
	%Description:
	%
	%
	%Usage:
	%	[ constr , y ] = cg.polytope_is_empty_set( H , h )
	%
	%Todo:
	%	Eventually make eps0 a potential input to the function.

	%% Constants
	n_H = size(H,1);

	eps0 = 10^(-4);

	%% Create Constraints in YALMIP
	y = sdpvar(n_H,1);

	constraint_out = [ H'*y == 0] + [ h'*y <= -eps0 ] + [y >= 0];

end