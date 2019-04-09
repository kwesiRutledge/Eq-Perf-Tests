function [S] = calc_u_effect_mat(varargin)
	%calc_u_effect_mat
	%	Description:
	%		Calculates the matrix that describes how the trajectory of inputs
	%		[u(t_0),u(t_0+1),...,u(T-1)] affects the state trajectory
	%		[x(t_0);x(t_0+1); ... ;x(T-1);x(T)]
	%
	%	Usage:
	%		S = calc_u_effect_mat(sys.A,sys.B,T);
	%		S = calc_u_effect_mat(sys_arr,L);

	%Select which version of the system to use by checking the type of the first argument.
	if isa(varargin{1},'Aff_Dyn')

		%Input Processing
		sys_arr = varargin{1};
		L = varargin{2};

		%Constants
		n = size(sys_arr(1).A,1);
		m = size(sys_arr(1).B,2);

		T = length(L);

		S = zeros(n*(T+1),m*T);

		nonzero_part = []; %Nonzero part of the row.

		for i = 1:T
			if i == 1 
				nonzero_part = [ sys_arr(L(i)).B ];
			else
				nonzero_part = [  sys_arr(L(i-1)).A*nonzero_part, sys_arr(L(i)).B ];
			end

			%Update G Matrix
			S([i*n+1:(i+1)*n],[1:m*i]) = nonzero_part;

		end


	elseif isnumeric(varargin{1})		

		%Input Processing
		A = varargin{1};
		B = varargin{2};
		T = varargin{3};

		%Constants
		S = zeros(size(A,1)*(T+1),size(B,2)*T);

		for i = 1: T

			temp = [];
			for k = 1:i
				temp = [temp A^(i-k)*B];
			end
			
			S([i*size(A,1)+1:(i+1)*size(A,1)],[1:size(B,2)*i]) = temp;
		
		end

	else

		error('The first input to ''calc_u_effect_mat()'' is of unrecognized type.')

	end

end