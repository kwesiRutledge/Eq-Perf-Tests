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

		function [poly_out] = to_poly(obj)
			%to_poly_v2.m
			%	Description:
			%		This function is meant to transform a Zonotope into a polyhedron.

			%% Constants
			% num_g = size(obj.G,2);

			warning('We are unsure how this algorithm works when G is not square.')

			%% Output

			V = []; 						%Set of vertices
			for vert_ind = 0:2^obj.num_g-1
				gen_wghts = dec2bin(vert_ind,obj.num_g)';
				gen_wghts = double(gen_wghts == '1') + (-1)*double(gen_wghts == '0');

				V(:,vert_ind+1) = obj.G*gen_wghts;
			end

			poly_out = Polyhedron('V',V')

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