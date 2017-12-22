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
	a = [0.5 3.5 0.25];

	example_sys.A =	[	0		1		0	;
						0		0		1	;
						-a(1)	-a(2)	-a(3)	];
	n = size(example_sys.A,1);

	example_sys.B = eye(n);
	example_sys.C = [ 0 1 0 ; 0 0 1 ];
	example_sys.E = [1 0 0]';%example_sys.B; 		%Noise occurs on the input signal

	example_sys.G = eye(size(example_sys.C,1));

	example_sys.m = 0; %0.5*perf_level;
	example_sys.d = 0.75*perf_level;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Attempt to Find Controller that Guarantees Equalized Performance [TRIANGLE INEQUALITY] %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	L_ti = sdpvar(size(example_sys.A,1),size(example_sys.C,1),'full');

	obj_ti = norm( example_sys.A + L_ti*example_sys.C , Inf) * perf_level + norm(L_ti * [ 0 0 ; 0 1],Inf)*example_sys.m + norm(example_sys.E,Inf)*example_sys.d;

	ti_sol.optim_flags = optimize([] , obj_ti , sdpsettings('verbose',verbosity) );

	if ti_sol.optim_flags.problem == 0
		disp(['Triangle Inequality Optimization Solved'])
	else
		error(['Triangle Inequality Optimization NOT Solved.'])
	end

	%Save Optimal Objective's Value
	ti_sol.opt_obj = value(obj_ti);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Attempt to Find Controller that Guarantees Equalized Performance [SKAF AND BOYD, G MATRIX INCLUDED] %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	n = size(example_sys.A,1);
	t_horizon = 1;

	% Constraints 
	%------------

	l_diag_constr = [];
	robust_constr = [];
	epi_constr = [];

	%Make x0 an optimization Variable
	example_sys.x0 = sdpvar(size(example_sys.A,1),1,'full');

	%Create constraints on x0 (based on performance level)
	robust_constr = robust_constr + [ -perf_level <= example_sys.x0 <= perf_level , uncertain(example_sys.x0) ];

	% Create Constant Matrices Based On Model
	%----------------------------------------

	[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(example_sys,t_horizon);

	if verbosity >= 1
		disp('Created Skaf and Boyd Matrices.')
	end

	% Perform Robust Optimization Using YALMIP "uncertain"
	%-----------------------------------------------------
	if verbosity >= 1
		disp('YALMIP Robust Optimization: ')
	end

	Q = sdpvar(size(H,2),size(Cm,1),'full');
	r = sdpvar(size(H,2),1,'full');

	%Disturbance vectors
	w = sdpvar(n*t_horizon,1,'full');
	v = sdpvar(size(example_sys.C,1)*t_horizon,1,'full');

	if verbosity >= 1
		disp('- Created Optimization Variables.')
	end

	%Create Expressions Containing Q,r for Optimization
	Pxw = (eye(n*(t_horizon+1))+H*Q*Cm)*G * diag(example_sys.E) ;
	Pxv = H*Q * example_sys.G;
	x_tilde = (eye(n*(t_horizon+1)) + H*Q*Cm)*x0m + H*r;

	%Create Objective
	if any( strcmp(varargin,'R') )
		R_string_loc = find( strcmp(varargin,'R') );
		R = varargin{R_string_loc+1};
	else
		R = [zeros(n,n*t_horizon) eye(n)]; %Create the standard selection matrix
	end

	size((x_tilde + Pxw * w + Pxv * v))

	objective = norm( R*(x_tilde + Pxw * w + Pxv * v) , Inf );

	if verbosity >= 1
		disp('- Created Objective.')
	end

	%Create Constraints

	%Q is lower diagonal.
	for bl_row_num = 1 : t_horizon-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(sys.B,2)+1:bl_row_num*size(sys.B,2)], ...
												[bl_row_num*size(sys.C,1)+1:end] ) == 0 ];
	end

	if verbosity >= 2
		l_diag_constr
	end

	%Robustifying against w and v
	robust_constr = robust_constr + [ -example_sys.m <= v <= example_sys.m , uncertain(v) ];
	robust_constr = robust_constr + [ -example_sys.d <= w <= example_sys.d , uncertain(w) ];

	if verbosity >= 2
		robust_constr
	end

	%Epigraph Constraints
	alpha0 = sdpvar(1,1,'full');
	epi_constr = [ objective <= alpha0 ];

	if verbosity >=2 
		epi_constr
	end

	if verbosity >= 1
		disp('- Created Constraints.')

		%Solve Optimization
		disp('Solving YALMIP Robust Optimization...')
	end

	op_num = verbosity;
	ops = sdpsettings('verbose',op_num);
	robust_w_g.optim_flags = optimize(l_diag_constr+robust_constr+epi_constr,alpha0,ops);

	if verbosity >= 1
		if robust_w_g.optim_flags.problem == 0
			disp(['YALMIP Robust Optimization Solved'])
		else
			%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
			disp(['YALMIP Robust Optimization NOT Solved.'])
		end
	end

	%Save results
	robust_w_g.Q = value(Q);
	robust_w_g.r = value(r);
	robust_w_g.opt_obj = value(alpha0);

	robust_w_g.F = value( (pinv(value(eye(size(Q,1)) + Q*Cm*H)) ) * Q);
	robust_w_g.u0 = value((eye(size(robust_w_g.F,1)) + robust_w_g.F*Cm*H) * r);	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Attempt to Find Controller that Guarantees Equalized Performance [SKAF AND BOYD] %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	opt_obsv = generate_skaf_controller( example_sys , 1 , verbosity , 'PL' , perf_level );

	%%%%%%%%%%%%%%%%%%
	%% Show Results %%
	%%%%%%%%%%%%%%%%%%
	disp('=====================')
	disp('Experiment 9 Results:')
	disp(' ')
	disp(['Triangle Inequality Optimal Objective:               ' num2str(ti_sol.opt_obj) ])
	disp(['Optimal Observer (with G?) Optimal Objective:        ' num2str(robust_w_g.opt_obj)])
	disp(['Function (Doesn''t consider Measurement Noise Shape): ' num2str(opt_obsv.opt_obj)])

	%%%%%%%%%%%%%%%%%%
	%% Save Results %%
	%%%%%%%%%%%%%%%%%%

	results.ti_sol = ti_sol;
	results.opt_obsv_sol = opt_obsv;
	results.robust_w_g = robust_w_g;

end
