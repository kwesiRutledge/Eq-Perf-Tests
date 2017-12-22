function [ results ] = observer_comparison5( varargin )
	%	observer_comparison5.m
	%		The goal of this test is to observe the effect of using the
	%		dynamic output feedback controller's initial condition as an
	%		optimization variable in yong's method.

	%%%%%%%%%%%%%%%%%%%
	%% Manage inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin == 0
		verbosity	= 0;
	elseif nargin == 1
		verbosity	= varargin{1};
	elseif nargin == 2
		verbosity	= varargin{1};
		dim_s		= varargin{2};
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	t_horizon = 7;
	
	perf_level = 1;

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	n = size(acc.A,1);

	if ~exist('dim_s')
		dim_s = n;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Optimize, with s0 as an optimization variable. %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('++++++++++++++++++++++++++++++++++++++++++++++++++')
	disp('Running MODIFIED Dynamic Output Feedback Observer.')
	disp(' ')

	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));	%For error system, the input we are interested in
									% isn't modified by a B matrix
	%acc_e.x0 = zeros(n,1);		%The initial condition needs to be set due to my rigid,
								% old way of doing things. It will not be used when pl is given.

	acc_e.x0 = sdpvar(n,1,'full');

	%Construct Optimization Matrices
	%-------------------------------

	s0 = sdpvar(dim_s,1,'full');
	dyn_obs_sys = dyn_obs_ify(acc_e,dim_s,s0);

	[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(dyn_obs_sys,t_horizon);

	if verbosity >= 1
		disp('Created Skaf Constants')
	end

	if verbosity >= 1
		disp('YALMIP Robust Optimization: ')
	end

	Q = sdpvar(size(H,2),size(Cm,1),'full');
	r = sdpvar(size(H,2),1,'full');

	w_y = sdpvar(n,t_horizon,'full');
	v_y = sdpvar(size(acc_e.C,1),t_horizon,'full');

	%Make w and v into vectors
	w = []; v = [];
	for t = 1 : t_horizon
		w = [ w ; w_y(:,t) ; zeros(dim_s,1) ];
		v = [ v ; zeros(dim_s,1) ; v_y(:,t) ];
	end

	if verbosity >= 1
		disp('- Variables Created')
	end

	% Create Compound Expressions with YALMIP Variables
	%--------------------------------------------------
	Pxw = (eye((n+dim_s)*(t_horizon+1))+H*Q*Cm)*G;
	Pxv = H*Q;
	x_tilde = (eye((n+dim_s)*(t_horizon+1)) + H*Q*Cm)*x0m + H*r;

	R = [zeros(n,(n+dim_s)*t_horizon) eye(n) zeros(n,dim_s)];

	objective = norm( R*(x_tilde + Pxw * w + Pxv * v) , Inf );

	if verbosity >= 1
		disp('- Objective Created')
	end

	% Create YALMIP Optimization's Constraints
	%-----------------------------------------

	l_diag_constr = []; robust_constr = []; epi_constr = []; disturb_constrs = [];

	%Lower Diagonal Constraint

	for bl_row_num = 1 : t_horizon-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(dyn_obs_sys.B,2)+1:bl_row_num*size(dyn_obs_sys.B,2)], ...
												[bl_row_num*size(dyn_obs_sys.C,1)+1:end] ) == 0 ];
	end

	if verbosity >= 2
		disp('Lower Diagonal Constraint:')
		l_diag_constr
	end

	%Robustifying against w0 and v0
	for t = 1 : t_horizon
		robust_constr = robust_constr + [ -acc_e.m <= v_y(:,t) <= acc_e.m , uncertain(v_y(:,t)) ];
		robust_constr = robust_constr + [ -acc_e.d <= w_y(:,t) <= acc_e.d , uncertain(w_y(:,t)) ];
	end

	robust_constr = robust_constr + [ -perf_level <= dyn_obs_sys.x0 <= perf_level , uncertain(dyn_obs_sys.x0) ];

	if verbosity >= 2
		disp('Robust Constraints:')
		robust_constr
	end

	%Epigraph Constraint
	alpha0 = sdpvar(1,1,'full');
	epi_constr = [ objective <= alpha0 ];

	if verbosity >=2 
		disp('Epigraph Constraint:')
		epi_constr
	end

	if verbosity >= 1
		disp('- Constraints Created')
	end

	% Solve Optimization
	%-------------------

	ops = sdpsettings('verbose',verbosity);
	results.design_s0.opt_sol = optimize(l_diag_constr+robust_constr+epi_constr,alpha0,ops);

	if verbosity >= 1
		if results.design_s0.opt_sol.problem == 0
			disp(['YALMIP Optimization Solved'])
		else
			error(['YALMIP Optimization NOT Solved.'])
		end
	end

	disp(' ')
	disp('Finished MODIFIED Dynamic Output Feedback Observer.')
	disp('++++++++++++++++++++++++++++++++++++++++++++++++++')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform Similar Optimization with Standard Yong Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	results.curr_fcn = generate_yong_controller( acc_e , t_horizon , verbosity , 'PL' , perf_level );

	%%%%%%%%%%%%%%%%%%%%%
	%% Compare Results %%
	%%%%%%%%%%%%%%%%%%%%%

	disp(' ')
	disp('=======')
	disp('Results')
	disp(' ')
	disp([ 's0 as a variable: ' num2str(value(alpha0)) ])
	disp([ 's0 = ' num2str(value(s0')) ])
	disp([ 'Function (circa 2017.10.21): ' num2str(results.curr_fcn.opt_obj) ])

	disp(' ')
	disp('=======')

	%%%%%%%%%%%%%%%%%
	%% Saving Data %%
	%%%%%%%%%%%%%%%%%

	results.sys = acc_e;

	results.design_s0.G = G;
	results.design_s0.H = H;
	results.design_s0.Cm = Cm;
	results.design_s0.x0m = x0m;

	% Saving Results of Optimization
	%-------------------------------

	results.design_s0.Q = value(Q);
	results.design_s0.Q( isnan(value(Q)) ) = 0;
	results.design_s0.r = value(r);
	results.design_s0.r( isnan(value(r)) ) = 0;
	results.design_s0.opt_obj = value(alpha0);
	value(alpha0);

	if any(isnan(value(Q))) | any(isnan(value(r))) | any(isnan(value(alpha0)))
		%Warn the user
		warning('There are NaNs in Q,r,or opt_obj. This can be because of bad inputs or because some part of the optimization variables are unused. Setting NaNs to zero.')
	end

	results.design_s0.F = value( (pinv(value(eye(size(results.design_s0.Q,1)) + results.design_s0.Q*Cm*H)) ) * results.design_s0.Q);
	results.design_s0.u0 = value((eye(size(results.design_s0.F,1)) + results.design_s0.F*Cm*H) * results.design_s0.r);

	%Value of s0
	results.design_s0.s0 = value(s0);

end