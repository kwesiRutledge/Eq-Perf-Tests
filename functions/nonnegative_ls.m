function [x,fval,exitflag,output] = nonnegative_ls( varargin )
%nonnegative_ls.m
%Description:
%	Solves the nonnegative least squares problem, i.e.
%		arg min_x || Ax - y ||_2 subject to x >= 0
%	where the inequality on x holds elementwise.
%
%Usage:
%	[x,fval,exitflag,output] = nonnegative_ls( A , y )
%
%Output:
%	x 			- The argument that minimizes the objective.
%	fval 		- The value of the objective when the minimizer is used.
%	exitflag	- This flag gives information about the reason why quadprog() stopped
%				  while attempting to solve your problem.
%				  1 = Function converged to the solution x
%				  -2 = Problem is infeasible
%				  (See quadprog() documentation for more.)
%	output 		- Information about the optimization process

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	A = varargin{1};
	y = varargin{2};

	if nargin > 2
		argin_idx = 3;
		while argin_idx < nargin
			switch varargin{argin_idx}
				case 'verbosity'
					verbosity = varargin{argin_idx+1};
					argin_idx = argin_idx + 2;
				otherwise
					error('Unrecognized input.')
			end

		end
	end

	if ~exist('verbosity')
		verbosity = 0;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n_x = size(A,2);

	if verbosity == 0
		ops0 = optimoptions('quadprog','Display','none');
	end

	%% Reformulating Variables 

	Q = A'*A;
	c = -A'*y;

	[x, fval , exitflag, output] = quadprog( Q , c , -eye(n_x) , zeros(n_x,1) , ...
											 [] , [] , [] , [] , [] , ops0 );

	fval = norm(A*x-y,2);

end


