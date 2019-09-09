classdef LCSAS
	%Definition:
	%	Object for a Language Constrained Switched Affine System.
	%
	%Requires:
	%	Aff_Dyn() class.

	properties
		Dyn;
		n_modes;
	end

	methods

		function out_sys = LCSAS(ad_arr)
			%Description:
			%	Constructor for the LCSAS object.
			%
			%Inputs:
			%	ad_arr 	- Array of Aff_Dyn() objects.

			if ~isa(ad_arr,'Aff_Dyn')
				error('Expected the input to LCSAS to be an array of Aff_Dyn objects.')
			end

			%Algorithm
			out_sys.Dyn = ad_arr;
			out_sys.n_modes = length(ad_arr);

		end
	end

end