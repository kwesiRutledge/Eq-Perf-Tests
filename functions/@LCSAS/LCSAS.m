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
	%	- check
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
		U;	%Input Set
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
				case 'U'
					U_in = varargin{argin_idx+1};
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

			out_sys.U = [];
			if exist('U_in')
				out_sys.U = U_in;
			end

		end

		function n_x = Dim_x(lcsas)
			%Description:
			%	Computes the dimension of the state from the LCSAS.

			n_x = size(lcsas.Dyn(1).A,1);

		end

		function n_u = Dim_u(lcsas)
			%Description:
			%	Computes the dimension of the input from the LCSAS.
			%
			%Assumptions:
			%	This function assumes that the input is the same dimension for all modes in the system.
			
			n_u = size(lcsas.Dyn(1).B,2);

		end

		function n_y = Dim_y(lcsas)
			%Description:
			%	Computes the dimension of the output from the LCSAS.
			
			n_y = size(lcsas.Dyn(1).C,1);

		end

		function n_w = Dim_w(lcsas)
			%Description:
			%	Computes the dimension of the process noise/disturbance from the LCSAS.
			
			n_w = size(lcsas.Dyn(1).B_w,2);

		end

		function n_v = Dim_v(lcsas)
			%Description:
			%	Computes the dimension of the output noise/disturbance from the LCSAS.
			
			n_v = size(lcsas.Dyn(1).C_v,2);

		end

		function varargout = Dimensions(lcsas)
			%Description:
			%	Collects all of the interesting dimensions of the system.
			%Usage:
			%	[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions()	

			switch nargout
			case 1
				%Create struct which will contain all values
				varargout{1}.x = lcsas.Dim_x();
				varargout{1}.u = lcsas.Dim_u();
				varargout{1}.y = lcsas.Dim_y();
				varargout{1}.w = lcsas.Dim_w();
				varargout{1}.v = lcsas.Dim_v();

			case 5
				%Save each dimension to a new variable.
				varargout{1} = lcsas.Dim_x();
				varargout{2} = lcsas.Dim_u();
				varargout{3} = lcsas.Dim_y();
				varargout{4} = lcsas.Dim_w();
				varargout{5} = lcsas.Dim_v();
				
			otherwise
				error(['Called ''Dimensions'' with an improper number of output arguments: ' num2str(nargout) ' (expected 1 or 5).' ])
			end

		end

		function check(varargin)
			%check.m
			%Description:
			%	Checks certain aspects of the LCSAS object to make sure that they are defined properly.

			lcsas = varargin{1};

			for argin_idx = 2:nargin
				switch varargin{argin_idx}
					case 'U'
						if isempty(lcsas.U)
							error(['The field U is empty. Please define it as a Polyhedron.'])
						end

						if ~isa(lcsas.U,'Polyhedron')
							error(['The field U must be a Polyhedron.'])
						end

						argin_idx = argin_idx + 1;

					case 'X0'

						if isempty(lcsas.X0)
							error(['The field X0 is empty. Please define it as a Polyhedron.'])
						end

						if ~isa(lcsas.X0,'Polyhedron')
							error(['The field X0 must be a Polyhedron.'])
						end

						argin_idx = argin_idx + 1;

					case 'L'

						if ~isa(lcsas.L,'Language')
							error(['The field L is not a Language object in input LCSAS.'])
						end
						
					otherwise
						error(['Unexpected input to check: ''' varargin{argin_idx} '''!'])
				end
			end

			if nargin == 1
				disp('Nothing was given to check.')
			end

		end

		function n_modes = NumberOfModes(lcsas)
			%Description:
			%	Identifies the number of modes in the lcsas, by investigating each word in the language.

			%% Constants
			L = lcsas.L;

			%% Algorithm
			n_modes = -1;
			for word_index = 1:L.cardinality()
				%Find the maximum in each word.
				temp_word = L.words{word_index};

				max_mode_in_tw = max(temp_word);

				if max_mode_in_tw > n_modes 
					n_modes = max_mode_in_tw;
				end
			end

		end

		function W_T = ToWSequence(lcsas,temp_word)
			%Description:
			%

			%% Constants
			TimeHorizon = length(temp_word); 

			%% Algorithm
			W_T = 1;
			for t = 0:TimeHorizon-1
				q_t = temp_word(t+1);
				W_T = W_T * lcsas.Dyn(q_t).P_w;
			end

		end

		function PlotAllReachableSets(varargin)
			%Description:
			%	This function tries to plot all reachable sets for the LCSAS.

			%% Input Processing
			lcsas = varargin{1};
			lcsas.check('U','X0')

			%% Constants
			X0 = lcsas.X0;
			U  = lcsas.U;

			L = lcsas.L;
			T = length(L.words{1});

			[ n_x , n_u , n_y , n_w , n_v ] = lcsas.Dimensions();

			U_T = 1;
			for t = 0:T-1
				U_T = U_T * U;
			end

			%% Algorithm
			R_0_T = {};
			for word_index = 1:L.cardinality()
				% get mpc matrices for target_word
				target_word = L.words{word_index};
				TimeHorizon = length(target_word);
				[S_w_i,S_u_i,~,J_i,f_bar_i,B_w_block,~] = lcsas.get_mpc_matrices('word',target_word);

				% Create Reachable Set
				W_T = lcsas.ToWSequence(target_word);

				R_0_T{word_index} = S_w_i*B_w_block * W_T + S_u_i * U_T + J_i * X0 + S_w_i*f_bar_i;
			end

			% Plot
			for word_index = 1:L.cardinality()
				% get mpc matrices for target_word
				target_word = L.words{word_index};
				TimeHorizon = length(target_word);

				figure;
				for t = 0:TimeHorizon
					subplot(1,TimeHorizon+1,t+1)
					if t == 0
						plot(X0)
					else
						plot(R_0_T{word_index}.projection(n_x*t+[1:n_x]))
					end

				end

			end



		end

	end

end