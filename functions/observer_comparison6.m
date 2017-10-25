function [ results ] = observer_comparison6( varargin )
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

	T_max = 10;
	
	perf_level = 1;

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	n = size(acc.A,1);

	if ~exist('dim_s')
		dim_s = n;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Set Up Loops for Different Time Horizons %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));	%For error system, the input we are interested in
									% isn't modified by a B matrix
	%acc_e.x0 = zeros(n,1);		%The initial condition needs to be set due to my rigid,
								% old way of doing things. It will not be used when pl is given.

	acc_e.x0 = sdpvar(n,1,'full');

	for t_horizon = 1 : T_max

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Perform Similar Optimization with Standard Yong Controller %%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		results.curr_fcn(t_horizon) = generate_skaf_controller( acc_e , t_horizon , verbosity , 'PL' , perf_level );

		results.std_obsv_guarantees(t_horizon) = results.curr_fcn(t_horizon).opt_obj;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Applying ROO Requires A Slightly Modified Optimization %%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		%Create Reduced Order Observer's Error Dynamics
		acc_roo.A = acc.A(3,3);
		acc_roo.B = eye(1);
		acc_roo.C = acc.A([1:2],3);
		acc_roo.G = acc.E([1:2],1);
		% acc_roo.m = acc.m;
		acc_roo.d = acc.d; 
		acc_roo.x0 = sdpvar(1,1,'full');

		n_roo = size(acc_roo.A,1);

		%Create Skaf Matrices
		[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(acc_roo,t_horizon);

		Q = sdpvar(size(H,2),size(Cm,1),'full');
		r = sdpvar(size(H,2),1,'full');

		w = sdpvar(n_roo*t_horizon,1,'full');

		% Create Compound Expressions with YALMIP Variables
		%--------------------------------------------------

		%Calculate Big C Matrix
		G_at_each_n = {};
		for i = 1:t_horizon
			G_at_each_n{i} = acc_roo.G; 
		end

		G_big = [ blkdiag(G_at_each_n{:}) ];

		Pxw = (eye(n_roo*(t_horizon+1))+H*Q*Cm)*G + H*Q*G_big;
		%Pxv = H*Q;
		x_tilde = (eye(n_roo*(t_horizon+1)) + H*Q*Cm)*x0m + H*r;

		R = [ zeros(n_roo,n_roo*t_horizon) eye(n_roo) ];

		objective = norm( R*(x_tilde + Pxw * w) , Inf );

		% Create YALMIP Optimization's Constraints
		%-----------------------------------------

		l_diag_constr = []; robust_constr = []; epi_constr = [];

		%Lower Diagonal Constraint

		for bl_row_num = 1 : t_horizon-1
			l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_roo.B,2)+1:bl_row_num*size(acc_roo.B,2)], ...
													[bl_row_num*size(acc_roo.C,1)+1:end] ) == 0 ];
		end

		%Robustifying against w0 and v0
		robust_constr = robust_constr + [ -acc_roo.d <= w <= acc_roo.d , uncertain(w) ];
		robust_constr = robust_constr + [ -perf_level <= acc_roo.x0 <= perf_level , uncertain(acc_roo.x0) ];

		%Epigraph Constraint
		alpha0 = sdpvar(1,1,'full');
		epi_constr = [ objective <= alpha0 ];

		ops = sdpsettings('verbose',verbosity);
		results.roo(t_horizon) = optimize(l_diag_constr+robust_constr+epi_constr,alpha0,ops);
		results.roo_guarantees(t_horizon) = value(alpha0);
		if verbosity >= 1
			if results.roo.problem == 0
				disp(['YALMIP Optimization Solved'])
			else
				error(['YALMIP Optimization NOT Solved.'])
			end
		end

	end

	figure;
	hold on;
	plot(results.roo_guarantees)
	plot(results.std_obsv_guarantees)

	% %%%%%%%%%%%%%%%%%%%%%
	% %% Compare Results %%
	% %%%%%%%%%%%%%%%%%%%%%

	% disp(' ')
	% disp('=======')
	% disp('Results')
	% disp(' ')
	% disp([ 's0 as a variable: ' num2str(value(alpha0)) ])
	% disp([ 's0 = ' num2str(value(s0')) ])
	% disp([ 'Function (circa 2017.10.21): ' num2str(results.curr_fcn.opt_obj) ])

	% disp(' ')
	% disp('=======')

	% %%%%%%%%%%%%%%%%%
	% %% Saving Data %%
	% %%%%%%%%%%%%%%%%%

	% results.sys = acc_e;

	% results.design_s0.G = G;
	% results.design_s0.H = H;
	% results.design_s0.Cm = Cm;
	% results.design_s0.x0m = x0m;

	% % Saving Results of Optimization
	% %-------------------------------

	% results.design_s0.Q = value(Q);
	% results.design_s0.Q( isnan(value(Q)) ) = 0;
	% results.design_s0.r = value(r);
	% results.design_s0.r( isnan(value(r)) ) = 0;
	% results.design_s0.opt_obj = value(alpha0);
	% value(alpha0);

	% if any(isnan(value(Q))) | any(isnan(value(r))) | any(isnan(value(alpha0)))
	% 	%Warn the user
	% 	warning('There are NaNs in Q,r,or opt_obj. This can be because of bad inputs or because some part of the optimization variables are unused. Setting NaNs to zero.')
	% end

	% results.design_s0.F = value( (pinv(value(eye(size(results.design_s0.Q,1)) + results.design_s0.Q*Cm*H)) ) * results.design_s0.Q);
	% results.design_s0.u0 = value((eye(size(results.design_s0.F,1)) + results.design_s0.F*Cm*H) * results.design_s0.r);

	% %Value of s0
	% results.design_s0.s0 = value(s0);

end