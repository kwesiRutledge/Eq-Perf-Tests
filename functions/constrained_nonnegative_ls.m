function [x,fval,exitflag,output] = constrained_nonnegative_ls( varargin )
%nonnegative_ls_withconstraints.m
%Description:
%	Solves the nonnegative least squares problem, i.e.
%		arg min_x || Ax - y ||_2
%		subject to 	x >= 0 
%					(some constraints on x)
%	where the inequality on x holds elementwise.
%
%Usage:
%	[x,fval,exitflag,output] = nonnegative_ls_withconstraints( A , y , H_x , h_x , He_x , he_x )
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
	H_x = varargin{3};
	h_x = varargin{4};
	He_x = varargin{5};
	he_x = varargin{6};

	if nargin > 6
		argin_idx = 7;
		while argin_idx < nargin
			switch varargin{argin_idx}
				case 'verbosity'
					verbosity = varargin{argin_idx+1};
					argin_idx = argin_idx + 2;
				otherwise
					error('Unrecognized input.')
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
		ops_infty = optimoptions('quadprog','Display','none');
	else
		ops_infty = optimoptions('quadprog');
	end

	%% Reformulating Variables 

	Q = A'*A;
	c = -A'*y;

	[x, ~ , exitflag, output] = quadprog( Q , c , [-eye(n_x) ; H_x ] , [zeros(n_x,1); h_x] , ...
											 [] , [] , [] , [] , [] , ops_infty );

	if verbosity ~= 0
		disp(['exitflag = ' num2str(exitflag) ])
	end

	fval = norm(A*x-y,2);

end