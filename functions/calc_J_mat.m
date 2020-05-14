function [J] = calc_J_mat(varargin)
	%calc_J_mat
	%	Description:
	%		Calculates the matrix that describes how the initial state x(t_0)
	%		affects the state trajectory [x(t_0);x(t_0+1); ... ;x(T-1);x(T)]
	%
	%	Usage:
	%		J = calc_x0_mat(sys.A,T)
	%		J = calc_x0_mat(sys_arr,sig);

	%Constants

	%%Input Checking

	if nargin ~= 2
		error('Algorithm was made to accept only 2 arguments.')
	end

	%% Algorithm

	if isa(varargin{1},'Aff_Dyn')

		%Process inputs
		sys_arr = varargin{1};
		sig = varargin{2};

		%Constants
		n = size(sys_arr(1).A,1);

		%Creating matrix
		J = eye(n);

		for i = 1:length(sig)
			J = [J; sys_arr(sig(i)).A*J(end-n+1:end,:)];
		end

	elseif isnumeric(varargin{1})

		%Process Inputs
		A = varargin{1};
		T = varargin{2};

		%Creating matrix.
		J = eye(size(A,1));

		for i = 1 : (T-1)
			J = [J; (A^i)];
		end
	else
		error('Unexpected data type given as first input.')
	end
end