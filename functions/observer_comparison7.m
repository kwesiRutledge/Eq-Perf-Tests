function [ results ] = observer_comparison7( varargin )
	%	observer_comparison7.m
	%		The goal of this test is to discover

	%%%%%%%%%%%%%%%%%%%
	%% Manage inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin == 0
		verbosity	= 0;
	elseif nargin == 1
		verbosity	= varargin{1};
	elseif nargin == 2
		verbosity	= varargin{1};
		T_max		= varargin{2};
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%t_horizon = 7;
	if ~exist('T_max')
		T_max = 5;
	end

	perf_level = 1;

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	n = size(acc.A,1);

	acc_e = acc;
	acc_e.B = eye(n);

	for t_horizon = 1 : T_max

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Worst Case Estimation when No Input given to Observer System %%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		%Create Variables
		xi_0 	= sdpvar(n,1,'full');
		delta 	= sdpvar(n,t_horizon,'full');
		mu		= sdpvar(n,t_horizon,'full');
		alpha1 	= sdpvar(1,1,'full');

		if verbosity >= 1
			disp('- Created Variables')
		end

		% Create Constraints
		%-------------------

		%Objective

		partial_obj = (acc_e.A^t_horizon)*xi_0;
		for k = 0:t_horizon-1
			partial_obj = partial_obj + (acc_e.A^k)*delta(:,k+1);
		end
		objective = norm( partial_obj , Inf );

		objective_constr = [ objective <= alpha1 ];

		%Disturbance and ICs

		disturb_constrs = [];
		for k = 1 : t_horizon
			disturb_constrs = disturb_constrs + [ -acc_e.d <= delta(:,k) <= acc_e.d , uncertain(delta(:,k)) ];
		end

		ic_constrs = [ -perf_level <= xi_0 <= perf_level, uncertain(xi_0) ];

		if verbosity >= 1
			disp('- Created Constraints')
		end

		sol = optimize( objective_constr + disturb_constrs + ic_constrs , alpha1 , sdpsettings('verbose',verbosity) );

		%%%%%%%%%%%%%%%%%%
		%% Save Results %%
		%%%%%%%%%%%%%%%%%%

		results.sol(t_horizon) = sol;
		results.opt_val(t_horizon) = value(alpha1);

	end

	results.sys = acc_e;

	%%%%%%%%%%%%%%%%%%
	%% Plot Results %%
	%%%%%%%%%%%%%%%%%%
	figure;
	plot(results.opt_val)

	xlabel('Number of Packets Lost (N)')
	ylabel('Error Magnitude $||e(N)||_{\infty}$','Interpreter','latex')

	title('Worst Case Error when Observer gives no Input when Packets Lost')

end
