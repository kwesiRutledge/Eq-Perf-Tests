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
		domain;
	end

	methods

		function out_sys = LCSAS(varargin)
			%Description:
			%	Constructor for the LCSAS object.
			%
			%Inputs:
			%	ad_arr 	- Array of Aff_Dyn() objects.
			%	Domain  - A 
			%
			%Usage:
			%	out_sys = LCSAS(ad_arr,L)
			%	out_sys = LCSAS(ad_arr,L,Domain)

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			switch nargin
			case 2
				ad_arr = varargin{1};
				L = varargin{2};
			case 3
				ad_arr = varargin{1};
				L = varargin{2};
				domain = varargin{3};
			otherwise
				error(['Unexpected number of input arguments: ' num2str(nargin) ])
			end

			if ~isa(ad_arr,'Aff_Dyn')
				error('Expected the input to LCSAS to be an array of Aff_Dyn objects.')
			end

			if ~isa(L,'Language')
				error('Expected the input to be a Language object.')
			end

			if exist('domain')
				if ~isa(domain,'Polyhedron')
					error('Domain must be a Polyhedron')
				end
			end


			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			out_sys.Dyn = ad_arr;
			out_sys.n_modes = length(ad_arr);

			out_sys.L = L;
			if exist('domain')
				out_sys.domain = domain;
			end

		end
	end

end