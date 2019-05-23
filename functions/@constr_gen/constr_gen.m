classdef constr_gen
	%Description:
	%	This class abstracts away how we make constraints for the equalized recovery problems.
	properties

	end

	methods
		%Constructor
		function [generator] = constr_gen()
			%Description:
			%	The constructor for this class does nothing.

			disp('constr_gen() does not initialize any member variables.')
			% generator = [];

		end

		function [ dual_vars , constraint ] =  get_zonot_inclusion_constr(obj,Z_in,Z_circ,A,b,theta_2)
			%get_zonot_inclusion_constr.m
			%Description:
			%	Creates a zonotope inclusion constraint based on the desired property:
			%		A*Z_in + b \subseteq Z_circ
			%	This uses the interpretation of inclusion that our group developed in the Jupyter Notebook XXX.ipynb.
			%
			%Inputs:
			%	Z_in

			%Create Constraint for Upper Inclusion (theta_2 respects upper bound)
			Nu1 = sdpvar(Z_in.dim, Z_circ.num_g, 'full')
			Nu2 = sdpvar(Z_circ.dim, Z_circ.num_g, 'full');
			Lambda1 = sdpvar( 2*Z_in.num_g , Z_circ.num_g , 'full' );

			constr_contain1 = [ Nu1'*Z_in.c - Nu2'*(b-Z_circ.c) + Lambda1'*ones(2*Z_in.num_g,1) <= theta_2*ones(Z_circ.num_g,1) ];
			constr_contain1 = constr_contain1 + [ [zeros(Z_circ.num_g,Z_in.dim+Z_in.num_g) eye(Z_circ.num_g)] == Nu1'*[eye(Z_in.dim), -Z_in.G, zeros(Z_in.dim,Z_circ.num_g)] + ...
																													Nu2'*[A zeros(Z_circ.dim,Z_in.num_g) -Z_circ.G] + ...
																													Lambda1'*[zeros(Z_in.num_g*2,Z_in.dim) [eye(Z_in.num_g); -eye(Z_in.num_g) ] zeros(Z_in.num_g*2,Z_circ.num_g) ] ];

			Nu3 = sdpvar(Z_in.dim, Z_circ.num_g, 'full');
			Nu4 = sdpvar(Z_circ.dim, Z_circ.num_g, 'full');
			Lambda2 = sdpvar( 2*Z_in.num_g , Z_circ.num_g , 'full' );

			constr_contain2 = [Nu3'*Z_in.c - Nu4'*(b-Z_circ.c) + Lambda2'*ones(2*Z_in.num_g,1) <= theta_2*ones(Z_circ.num_g,1)];
			constr_contain2 = constr_contain2 + [ - [zeros(Z_circ.num_g,Z_in.dim+Z_in.num_g) eye(Z_circ.num_g)] == Nu3'*[eye(Z_in.dim) -Z_in.G zeros(Z_in.dim,Z_circ.num_g)] + ...
																													Nu4'*[A zeros(Z_circ.dim,Z_in.num_g) -Z_circ.G] + ...
																													Lambda2'*[zeros(Z_in.num_g*2,Z_in.dim) [eye(Z_in.num_g); -eye(Z_in.num_g) ] zeros(Z_in.num_g*2,Z_circ.num_g) ] ];

			pos_constr = [Lambda1 >= 0] + [Lambda2 >= 0];

			%Results
			dual_vars = {Nu1,Nu2,Nu3,Nu4,Lambda1,Lambda2};
			constraint = pos_constr + constr_contain1 + constr_contain2;

		end

	end
end