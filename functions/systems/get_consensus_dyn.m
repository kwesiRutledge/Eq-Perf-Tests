function [consensus_lcsas,des_pos] = get_consensus_dyn(varargin)
	%Description:
	%	Creates the dynamics for a set of point agents that:
	%	- Each are SINGLE Integrators
	%	- Have a set of consensus positions that are predefined.
	%
	%Example Usage:
	%	out_lcsas = get_consensus_dyn(n_x,n_y,dt)
	%	out_lcsas = get_consensus_dyn(n_x,n_y,dt,'disturb_params',P_w,P_v)
	%	out_lcsas = get_consensus_dyn(n_x,n_y,dt,'disturb_params',eta_w,eta_v)
	%
	%Notes:
	%	- Agent i is represented by the (2i-1) and the 2ith coordinates.
	%		+ The (2i-1)th coordinate is the error in the x position of agent i from its desired position.
	%		+ The (2i)th coordinate is the error in the y position of agent i from its desired position.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 3
		error('Need at least 3 arguments!')
	end

	n_x = varargin{1};
	n_y = varargin{2};
	dt = varargin{3};

	arg_idx = 4;
	while arg_idx < nargin
		switch varargin{arg_idx}
			case 'disturb_params'
				%User Provides Disturbance Set Info
				if isa(varargin{arg_idx+1},'Polyhedron')
					P_w = varargin{arg_idx+1};
					P_v = varargin{arg_idx+2}
				elseif isscalar(varargin{arg_idx+1})
					eta_w = varargin{arg_idx+1};
					eta_v = varargin{arg_idx+2};

					P_w = eta_w * Polyhedron('lb',-ones(1,2),'ub',ones(1,2));
					P_v = eta_v * Polyhedron('lb',-ones(1,n_x*n_y*2),'ub',ones(1,n_x*n_y*2));
				else
					error('Unrecognized datatype given for disturbance set.')
				end

				arg_idx = arg_idx+3;
			case 'L'
				L = varargin{arg_idx+1};
				arg_idx = arg_idx + 2;
				if ~isa(L,'Language')
					error('The input L is not a Lanugage object.')
				end
			otherwise
				error(['Unexpected input: ' varargin{arg_idx} ])
		end
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%
	num_agents = n_x*n_y;
	n = num_agents*2;

	dim = 2;

	unit_box_x = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));
	unit_box_w = Polyhedron('lb',-ones(1,dim),'ub',ones(1,dim));
	unit_box_v = unit_box_x;

	if ~exist('P_w')
		eta_w = 0.5; eta_v = 0.2;
		P_w = eta_w * unit_box_w;
		P_w2 = P_w + eta_w*ones(P_w.Dim,1);
		P_v = eta_v * unit_box_v;
	end

	r=2; %Spacing Between every row.

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Desired Positions %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	des_width = (n_x-1)*r;
	des_height = (n_y-1)*r;

	des_pos = [repmat(linspace(-des_width/2,des_width/2,n_x),1,n_y);
			   kron(linspace(-des_height/2,des_height/2,n_y),ones(1,n_x) )];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Dynamics Array %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	A = eye(n);

	B = dt*eye(n);
	B_w = repmat(eye(2),num_agents,1);

	C = eye(n);
	C_v = eye(n);

	consens_dyn1 = Aff_Dyn(A,B,zeros(n,1),C,P_w,P_v,B_w,C_v);
	consens_dyn2 = consens_dyn1;
	consens_dyn2.P_w = P_w2;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Creating the Mode Language L %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ~exist('L')

		lang_length = 6;
		strongw_start = 4;
		word_arr = {};
		for word_idx = 1:(lang_length-strongw_start+1)
			word_arr{word_idx} = [ones(1,strongw_start-1+(word_idx-1)),2*ones(1,lang_length-strongw_start+1-(word_idx-1))];
		end
		L = Language(word_arr);

	end

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	consensus_lcsas = LCSAS([consens_dyn1,consens_dyn2],L);

end
