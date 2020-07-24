classdef LCSAS
	%Definition:
	%	Object for a Language Constrained Switched Affine System.
	%
	%Construction:
	%	lcsas = LCSAS(ad_arr,L);
	%	lcsas = LCSAS(ad_arr,L,Domain);	
	%
	%Member functions:
	%	- consistency_set
	%	- get_mpc_matrices
	%
	%Member Variables:
	%	- Dyn
	%	- n_modes
	%	- L
	%	- domain (*)
	%	- X0 (*)
	%	* = This member variable may or may not be defined in the variable instance.
	%		You should check to see if it exists before using it.
	%
	%Requires:
	%	Aff_Dyn() class.

	properties
		Dyn;
		n_modes;
		L;
		domain;
		X0; %X0
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
			%	out_sys = LCSAS(ad_arr,L,'Domain',Domain)
			%	out_sys = LCSAS(ad_arr,L,'X0',X0)

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			if nargin >= 2
				ad_arr = varargin{1};
				L = varargin{2};
			end

			argin_idx = 3;

			if (nargin == 3) & (isa(varargin{3},'Polyhedron'))
				%Detected a Domain input
				domain = varargin{3};
				argin_idx = argin_idx + 1;
			end

			while argin_idx <= nargin
				switch varargin{argin_idx}
				case 'Domain'
					domain = varargin{argin_idx+1};
					argin_idx = argin_idx + 2;
				case 'X0'
					X0_in = varargin{argin_idx+1};
					argin_idx = argin_idx + 2;
				otherwise
					error(['Unexpected input to LCSAS: See Argument ' num2str(argin_idx) '.' ])
				end
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

			out_sys.X0 = [];
			if exist('X0_in')
				out_sys.X0 = X0_in;
			end

		end
	end

end