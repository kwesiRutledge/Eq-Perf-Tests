classdef Zonotope
	%Description:
	%
	properties
		G;
		c;
		dim;
		num_g;
	end

	methods
		function [Z] = Zonotope(G_in,c_in)
			%Description:
			%	Defines the set of generators and the center associated with the desired zonotope.
			
			%%%%%%%%%%%%%%%%%%
			%%Manage Inputs %%
			%%%%%%%%%%%%%%%%%%

			Z.G = G_in;
			Z.c = c_in;

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Management %%
			%%%%%%%%%%%%%%%%%%%%%%

			if size(G_in,1) ~= size(c_in)
				error('The size of G and c do not match.')
			end

			%Define variables on the periphery.
			[Z.dim,Z.num_g] = size(G_in);
		end

		function [poly_out] = to_poly(obj,method_num)
			%to_poly_v2.m
			%	Description:
			%		This function is meant to transform a Zonotope into a polyhedron.

			%% Constants
			% num_g = size(obj.G,2);

			if size(obj.G,1) ~= size(obj.G,2)
				warning('We are unsure how this algorithm works when G is not square.')
			end

			if ~exist('method_num')
				method_num = 2;
			end

			%% Output

			switch method_num
			case 1
				V = []; 						%Set of vertices
				for vert_ind = 0:2^obj.num_g-1
					gen_wghts = dec2bin(vert_ind,obj.num_g)';
					gen_wghts = double(gen_wghts == '1') + (-1)*double(gen_wghts == '0');

					V(:,vert_ind+1) = obj.G*gen_wghts;
				end

				poly_out = Polyhedron('V',V');

			case 2
				%Compute Polyhedron using the affine map function of MPT3
				poly_out = obj.G * Polyhedron('lb',-ones(1,obj.dim),'ub',ones(1,obj.dim)) + obj.c;
			otherwise
				error('Unrecognized method number.')
			end
		end

		function [] = plot(obj)
			%Description:
			%	

			%% Constants
			temp_poly = Z.to_poly();

			%% Plot Polyhedron
			plot(temp_poly);
		end
	end
end