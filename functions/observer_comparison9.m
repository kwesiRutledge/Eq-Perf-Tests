function [ results ] = observer_comparison9( varargin )
	%	observer_comparison9.m
	%		The goal of this test is to provide an example for which the triangle inequality
	%		cannot find a controller which achieves Equalized Performance while the Skaf method
	%		can find it.

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

	perf_level = 1;

	%Define our example system
	a = [2 3 1]

	example_sys.A =	[	0		1		0	;
						0		0		1	;
						-a(1)	-a(2)	-a(3)	];
	n = size(example_sys.A,1);

	example_sys.B = eye(n);
	example_sys.C = eye(n);
	example_sys.E = [0 0 1]';%example_sys.B; 		%Noise occurs on the input signal

	example_sys.m = 0;
	example_sys.d = perf_level;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Attempt to Find Controller that Guarantees Equalized Performance [TRIANGLE INEQUALITY] %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	L_ti = sdpvar(size(example_sys.C,1))

	obj_ti = norm( example_sys.A - L_ti*example_sys.C , Inf) * perf_level + norm(L_ti,Inf)*example_sys.m + norm(example_sys.E,Inf)*example_sys.d;

	ti_sol = optimize([] , obj_ti , sdpsettings('verbose',verbosity) );

	if ti_sol.problem == 0
		disp(['Triangle Inequality Optimization Solved'])
	else
		error(['Triangle Inequality Optimization NOT Solved.'])
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Attempt to Find Controller that Guarantees Equalized Performance [SKAF AND BOYD] %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	opt_obsv = generate_skaf_controller( example_sys , 1 , verbosity , 'PL' , perf_level );

	%%%%%%%%%%%%%%%%%%
	%% Plot Results %%
	%%%%%%%%%%%%%%%%%%
	

	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	results.ti_sol = ti_sol;
	results.opt_obsv_sol = opt_obsv;

end
