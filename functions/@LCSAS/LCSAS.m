classdef LCSAS
	%Definition:
	%	Object for a Language Constrained Switched Affine System.
	%
	%Requires:
	%	Aff_Dyn() class.

	properties
		Dyn;
		n_modes;
		L;
	end

	methods

		function out_sys = LCSAS(ad_arr,L)
			%Description:
			%	Constructor for the LCSAS object.
			%
			%Inputs:
			%	ad_arr 	- Array of Aff_Dyn() objects.

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if ~isa(ad_arr,'Aff_Dyn')
				error('Expected the input to LCSAS to be an array of Aff_Dyn objects.')
			end

			if ~isa(L,'Language')
				error('Expected the input to be a Language object.')
			end

			%Algorithm
			out_sys.Dyn = ad_arr;
			out_sys.n_modes = length(ad_arr);

			out_sys.L = L;

		end
	end

end