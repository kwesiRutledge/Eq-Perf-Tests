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
	%		[H,S,C_bar,J,f_bar,B_w_bar,C_v_bar] = get_mpc_matrices(sys_arr,T)
	%		[H,S,C_bar,J,f_bar] = get_mpc_matrices(sys_arr,sigma)
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
	%		sigma -		An array of integers.
    %                   This should be a single sequence of integer values
    %                   which defines the value of the discrete state at
    %                   each time state (and thus which dynamics object in
    %                   the sys_arr should be used).
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
	%
	%	Assumptions:
	%		We assume that the piecewise affine systems that are given as input maintain constant dimension of
	%		their disturbance (w), input (u), etc. That means that if u(t) is a 2 dimensional input at time 1
	%		then it always is.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if (nargin ~= 2)
		error(['Inappropriate number of arguments. (Received ' num2str(nargin) ')'])
	end

	sys_arr = varargin{1};
	sigma 	= varargin{2};

	if iscell(sigma)
		error(['Do not give input word sigma as a cell array.'])
	end

	%Verify that all systems in the system array are Aff_Dyn objects
	for sys_num = 1:length(sys_arr)
		if ~isa(sys_arr,'Aff_Dyn')
			error(['Entry #' num2str(sys_num) ' of the input system array is not of the type Aff_Dyn.'])
		end
	end

	if isscalar(sigma)
		sigma = ones(sigma,1);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Defining the Matrices %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Find H, S, and J Matrices
	H = calc_w_effect_mat(sys_arr,sigma);
	S = calc_u_effect_mat(sys_arr,sigma);
	J = calc_J_mat(sys_arr,sigma);
	
	%Calculate the big f matrix
	f_bar = [];
	for i = 1:length(sigma)
		f_bar = [f_bar; sys_arr(sigma(i)).f];
	end

	%Calculate Big C Matrix
	C_at_each_n = {};
	for i = 1:length(sigma)
		C_at_each_n{i} = sys_arr(sigma(i)).C; 
	end

	C_bar = [ blkdiag(C_at_each_n{:}) zeros( size(sys_arr(1).C) * [ length(sigma) 0 ; 0 1 ] ) ];

	%%%%%%%%%%%%%%%%%%%%%%
	%% Optional Outputs %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargout >= 6

		%Create E_bar
		E_at_each_n = {};
		for i = 1:length(sigma)
			E_at_each_n{i} = sys_arr(sigma(i)).B_w;
        end
        B_w_bar = blkdiag(E_at_each_n{:});

	end

	if nargout >= 7

		%Create G_bar
		G_at_each_n = {};
		for i = 1:length(sigma)
			G_at_each_n{i} = sys_arr(sigma(i)).C_v;
		end
		C_v_bar = blkdiag(G_at_each_n{:});	

	end

	%%%%%%%%%%%%%%%%%%%%
	%% Create Outputs %%
	%%%%%%%%%%%%%%%%%%%%

	varargout{1} = H;
	varargout{2} = S;
	varargout{3} = C_bar;
	varargout{4} = J;
	varargout{5} = f_bar;

	if nargout >= 6
		varargout{6} = B_w_bar;
	end

	if nargout >= 7
		varargout{7} = C_v_bar;
	end
end