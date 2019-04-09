function [x0_mat] = calc_x0_mat(varargin)
	%calc_x0_mat
	%	Description:
	%		Calculates the matrix that describes how the initial state x(t_0)
	%		affects the state trajectory [x(t_0);x(t_0+1); ... ;x(T-1);x(T)]
	%
	%	Usage:
	%		J = calc_x0_mat(sys.A,sys.x0,T)
	%		J = calc_x0_mat(sys_arr,L);

	%Constants

	if isa(varargin{1},'Aff_Dyn')

		%Process inputs
		sys_arr = varargin{1};
		L = varargin{2};

		%Constants
		n = size(sys_arr(1).A,1);

		%Creating matrix
		x0_mat = sys_arr(L(1)).x0;

		for i = 1:length(L)
			x0_mat = [x0_mat; sys_arr(L(i)).A*x0_mat(end-n+1:end,1)];
		end

	elseif isnumeric(varargin{1})

		%Process Inputs
		A = varargin{1};
		x0 = varargin{2};
		T = varargin{3};

		%Creating matrix.
		x0_mat = x0;

		for i = 1 : T
			x0_mat = [x0_mat; (A^i)* x0];
		end
	else
		error('Unexpected data type given as first input.')
	end
end