function [varargout] = get_mpc_matrices(varargin)
	% get_mpc_matrices
	%	Description:
	%		Creates MPC matrices that are needed for the affine switched systems in our work.
	%		This is currently designed to handle:
	%		- A missing observation system being part of the switched system
	%		- A system with unique disturbance matrices as being part of the switched system
	%		Modified slightly to work with LCSAS objects.
	%
	%	Usage:
	%		[ H , S , C_bar , J , f_bar ] = get_mpc_matrices(lcsas,'word',sigma)
	%		[ H , S , C_bar , J , f_bar , B_w_bar , C_v_bar ] = get_mpc_matrices(lcsas,'word',sigma)
	%		[ H_cell , S_cell , C_bar_cell , J_cell , f_bar_cell , B_w_bar_cell , C_v_bar_cell ] = get_mpc_matrices(lcsas,'All Words')
	%
	%	Inputs:
	%		T - 		Time horizon for the MPC matrices.
	%					When T is given (i.e. a single integer is given as the second input),
	%					it is assumed that the switching sequence is: 
	%		
	%		sigma -		An array of integers.
    %                   This should be a single sequence of integer values
    %                   which defines the value of the discrete state at
    %                   each time state (and thus which dynamics object in
    %                   the lcsas should be used).
    %
    %		lcsas -		A Language Constrained, Switched Affine System object representing the
    %					the switched system.
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

	if ~any(nargin == [2,3])
		error(['Inappropriate number of arguments. (Received ' num2str(nargin) ')'])
	end

	lcsas = varargin{1};
	in_str  = varargin{2};

	if nargin >= 3
		sigma = varargin{3};

		if iscell(sigma)
			error(['Do not give input word sigma as a cell array.'])
		end
	end

	switch in_str
	case 'All Words'
		varargout = get_set_of_mpc_matrices(lcsas,nargout);
		return;
	case 'time_horizon'
		sigma = ones(sigma,1);
	case 'word'
		sigma; %Do nothing
	otherwise
		error('Unexpected input string. Available options: ''All Words'', ''time_horizon'' and ''word''.')
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Defining the Matrices %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Find H, S, and J Matrices
	H = calc_w_effect_mat(lcsas.Dyn,sigma);
	S = calc_u_effect_mat(lcsas.Dyn,sigma);
	J = calc_J_mat(lcsas.Dyn,sigma);
	
	%Calculate the big f matrix
	f_bar = [];
	for i = 1:length(sigma)
		f_bar = [f_bar; lcsas.Dyn(sigma(i)).f];
	end

	%Calculate Big C Matrix
	C_at_each_n = {};
	for i = 1:length(sigma)
		C_at_each_n{i} = lcsas.Dyn(sigma(i)).C; 
	end

	C_bar = [ blkdiag(C_at_each_n{:}) zeros( size(lcsas.Dyn(1).C) * [ length(sigma) 0 ; 0 1 ] ) ];

	%%%%%%%%%%%%%%%%%%%%%%
	%% Optional Outputs %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargout >= 6

		%Create E_bar
		E_at_each_n = {};
		for i = 1:length(sigma)
			E_at_each_n{i} = lcsas.Dyn(sigma(i)).B_w;
        end
        B_w_bar = blkdiag(E_at_each_n{:});

	end

	if nargout >= 7

		%Create G_bar
		G_at_each_n = {};
		for i = 1:length(sigma)
			G_at_each_n{i} = lcsas.Dyn(sigma(i)).C_v;
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

function [ output_cell_array ] = get_set_of_mpc_matrices(lcsas,num_outputs)
	%Description:
	%	Returns collections of MPC Matrices for the entire set of words in the language lcsas.L

	% Find H S and J matrices
	H_cell = {}; S_cell = {}; J_cell = {}; f_bar_cell = {}; C_bar_cell = {};
	B_w_bar_cell = {}; C_v_bar_cell = {};

	for word_index = 1:lcsas.L.cardinality()

		switch num_outputs
			case 5
				[ H_cell{word_index} , S_cell{word_index} , C_bar_cell{word_index} , J_cell{word_index} , f_bar_cell{word_index} ] = get_mpc_matrices(lcsas,'word',lcsas.L.words{word_index});
			case 6
				[ H_cell{word_index} , S_cell{word_index} , C_bar_cell{word_index} , J_cell{word_index} , f_bar_cell{word_index} , B_w_bar_cell{word_index} ] = get_mpc_matrices(lcsas,'word',lcsas.L.words{word_index});
			case 7
				[ H_cell{word_index} , S_cell{word_index} , C_bar_cell{word_index} , J_cell{word_index} , f_bar_cell{word_index} , B_w_bar_cell{word_index} , C_v_bar_cell{word_index} ] = get_mpc_matrices(lcsas,'word',lcsas.L.words{word_index});
			otherwise
				error(['Unexpected number of outputs. Requested ' num2str(num_outputs) ' but can only provide 5,6, or 7.'])
		end
	end

	%Create output
	output_cell_array{1} = H_cell;
	output_cell_array{2} = S_cell;
	output_cell_array{3} = C_bar_cell;
	output_cell_array{4} = J_cell;
	output_cell_array{5} = f_bar_cell;

	if num_outputs >= 6
		output_cell_array{6} = B_w_bar_cell;
	end

	if num_outputs >= 7
		output_cell_array{7} = C_v_bar_cell;
	end



end