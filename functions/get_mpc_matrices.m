function [varargout] = get_mpc_matrices(varargin)
	% get_mpc_matrices
	%	Description:
	%		Creates MPC matrices that are needed for the affine switched systems in our work.
	%		This is currently designed to handle:
	%		- A missing observation system being part of the switched system
	%		- A system with unique disturbance matrices as being part of the switched system
	%
	%	Usage:
	%		[H,S,C_bar,J,f_bar] = get_mpc_matrices(sys_arr,T)
	%		[H,S,C_bar,J,E_bar,G_bar] = get_mpc_matrices(sys_arr,T)
	%		[H,S,C_bar,J,f_bar] = get_mpc_matrices(sys_arr,L)
	%
	%	Inputs:
	%		sys_arr - 	An array of structs, each containing system matrices and an initial 
	%					condition for the desired system.
	%					Required fields: .A,.B,.C,.x0
	%
	%		T - 		Time horizon for the MPC matrices.
	%					When T is given (i.e. a single integer is given as the second input),
	%					it is assumed that the switching sequence is: 
	%		
	%		L -			This is the sequence of discrete state values which determines the
	%					the affine update rule that is used in the 
	%		
	%
	%	Outputs:
	%		H - 	This is the matrix which defines how the process disturbances (w(t)) affect the
	%				system's state trajectory.
	%
	%		S - 	This is the matrix which defines how the inputs ( u(t) ) affect the system's
	%				state trajectory.
	%
	%		C_bar - ... defines how the state of the system is transmitted to the measurement trajectory.
	%
	%		J - 	... defines how the initial state (x(t_0)) affects the system's state trajectory.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if (nargin ~= 2)
		error(['Inappropriate number of arguments. (Received ' num2str(nargin) ')'])
	end

	sys_arr = varargin{1};
	L 		= varargin{2};

	%Verify that all systems in the system array are Aff_Dyn objects
	for sys_num = 1:length(sys_arr)
		if ~isa(sys_arr,'Aff_Dyn')
			error(['Entry #' num2str(sys_num) ' of the input system array is not of the type Aff_Dyn.'])
		end
	end

	if isscalar(L)
		L = ones(L,1);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Defining the Matrices %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Find H, S, and J Matrices
	H = calc_w_effect_mat(sys.A,T);
	S = calc_u_effect_mat(sys.A,sys.B,T);
	J = calc_x0_mat(sys.A,sys.x0,T);
	f_bar = kron(ones(T+1,1),sys)

	%Calculate Big C Matrix
	C_at_each_n = {};
	for i = 1:T
		if any(m_locs == (i-1))
			C_at_each_n{i} = zeros(size(sys.C));
		else
			C_at_each_n{i} = sys.C;
		end 
	end

	C_bar = [ blkdiag(C_at_each_n{:}) zeros( size(sys.C) * [ T 0 ; 0 1 ] ) ];

	%%%%%%%%%%%%%%%%%%%%%%
	%% Optional Outputs %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargout >= 5
		% Check to see if E/B_w or G/C_v exist
		% If not just error out.

		if (~isfield(sys,'E')) & (~isfield(sys,'B_w')) & (~isa(sys, 'Aff_Dyn'))
			error('E or B_w does not exist!')
		end

		if isfield(sys,'B_w') || isa(sys, 'Aff_Dyn')
			E = sys.B_w;
		else
			E = sys.E;
		end

		%Create E_bar
		E_at_each_n = {};
		for i = 1:T
			E_at_each_n{i} = E;
		end
		E_bar = blkdiag(E_at_each_n{:});

	end

	if nargout >= 6

		if ~isfield(sys,'G') & ~isfield(sys,'C_v') & (~isa(sys, 'Aff_Dyn'))
			error('G or C_v does not exist!')
		end

		if isfield(sys,'C_v') || isa(sys, 'Aff_Dyn')
			C_v = sys.C_v;
		else
			C_v = sys.G;
		end

		%Create G_bar
		G_at_each_n = {};
		for i = 1:T
			if any(m_locs == (i-1))
				G_at_each_n{i} = zeros(size(C_v));
			else
				G_at_each_n{i} = C_v;
			end 
		end
		G_bar = blkdiag(G_at_each_n{:});	

	end

	%%%%%%%%%%%%%%%%%%%%
	%% Create Outputs %%
	%%%%%%%%%%%%%%%%%%%%

	varargout{1} = H;
	varargout{2} = S;
	varargout{3} = C_bar;
	varargout{4} = J;

	if nargout >= 5
		varargout{5} = E_bar;
	end

	if nargout >= 6
		varargout{6} = G_bar;
	end
end