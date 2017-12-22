function [ results ] = observer_comparison12( varargin )
%	observer_comparison12.m
%		Description:
%			The objective of this experiment is to observe what happens for an arbitrary system
%			(where Equalized Performance is achievable) where we apply a very similar constraint
%			related to Equalized Recovery.
%
%			This Equalized Recovery Problem:
%				M1 = M2, T
%				Let ||xi(0)||<= M1. What is the minimum value for M1 such that,
%				- ||xi(t)||<= M2 \forall t \in [1,T-1]
%				- ||xi(T)||<= M1 
%				OR
%				- ||xi(t)||<= M1 \forall t
%
%		Inputs:


	%%%%%%%%%%%%%%%%%%%
	%% Manage Inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin == 0
		verbosity	= 0;
	elseif nargin == 1
		verbosity	= varargin{1};
	elseif nargin == 3
		verbosity	= varargin{1};
		T_missing	= varargin{2};
		T_available = varargin{3};
	elseif nargin == 4
		verbosity	= varargin{1};
		T_missing	= varargin{2};
		T_available = varargin{3};
		perf_level 	= varargin{4};
	else
		error('Unacceptable number of arguments.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%Default Values for observer_comparison
	if nargin < 3
		T_missing = 3;
		T_available = 3;
	end
	
	T = T_missing + T_available;

	%default Value for perf_level
	if ~exist('perf_level')
		perf_level = 1;
	end

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	n = size(acc.A,1);
	p = size(acc.C,1);

	%Create Error System
	acc_e = acc;
	acc_e.B = eye(n);

	d_u = size(acc_e.B,2);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Variables
	delta 		= sdpvar(n*T,1,'full');
	mu 			= sdpvar(p*T,1,'full');
	acc_e.x0 	= sdpvar(n,1,'full');

	pl_t		= sdpvar(1,1,'full');

	alpha_l 	= sdpvar(T+1,1,'full');

	% Create Trajectory Matrices
	[S,H,Cm,xi0m] = create_skaf_n_boyd_matrices(acc_e,T);

	% Create Objective (Epigraph constraint)
	Q = sdpvar(size(H,2),size(Cm,1),'full');
	r = sdpvar(size(H,2),1,'full');

	Pxd = (eye(n*(T+1))+S*Q*Cm)*S ;%E_bar*[ eye(b_dim*t_horizon) zeros(b_dim*t_horizon,b_dim) ];
	Pxm = S*Q;
	xi_tilde = (eye(n*(T+1)) + S*Q*Cm)*xi0m + S*r;

	R = [ zeros(n,n*T) eye(n) ];

	objective = norm( R*(xi_tilde + Pxd * delta + Pxm * mu) , Inf );
	interm_norms = norm([ eye(n*T) zeros(n*T,n) ]*( xi_tilde + Pxd * delta + Pxm * mu ) , Inf);
	% interm_vals = [ eye(n*T_available) zeros(n*T_available,n) ]*( xi_tilde + Pxd * delta(n*T_missing+1:end,1) + Pxm * mu(p*T_missing+1:end,1) );

	epi_constr2 = [ objective <= perf_level , interm_norms <= perf_level ];

	% Feasibility Constraints
	feasib_constrs = [];
	feasib_constrs = feasib_constrs + [-perf_level <= alpha_l <= perf_level];

	% Create Robustification Constraints

	robust_constrs = [];
	robust_constrs = robust_constrs + [ -acc_e.d <= delta <= acc_e.d , uncertain(delta) ];
	robust_constrs = robust_constrs + [ -acc_e.m <= mu <= acc_e.m , uncertain(mu)];
	robust_constrs = robust_constrs + [ -perf_level <= acc_e.x0 <= perf_level , uncertain(acc_e.x0) ];

	% Create Causality (Lower Diagonal) Constraint
	l_diag_constr = [];
	for bl_row_num = 1 : T_available-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
												[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
	end

	% Create Graph Constraints
	graph_constrs = [];
	for i = 0:T
		graph_constrs = graph_constrs + [ norm(select_m(i,T)*( xi_tilde + Pxd * delta + Pxm * mu ) , Inf) <= alpha_l(i+1) ];
	end

	% Optimize!
	ops = sdpsettings('verbose',verbosity);
	optim1 = optimize(epi_constr2+robust_constrs+l_diag_constr+feasib_constrs+graph_constrs ,pl_t+sum(alpha_l),ops);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Generate the same controller with Function %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[~,optim2] = achieve_eq_recovery_for(acc_e,T,perf_level,perf_level,verbosity);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Compare Results to Function %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('Experiment 1')
	disp('Do the control matrices synthesized by each method look alike?')

	if sum(sum(value(Q) == optim2.Q)) == prod(size(value(Q)))
		disp('Q matrices: IDENTICAL')
	else
		disp('Q matrices: DIFFERENT')
	end

	if sum(value(r) == optim2.r) == prod(size(value(r)))
		disp('r matrices: IDENTICAL')
	else
		disp('r matrices: DIFFERENT')
	end

	disp(' ')
	disp('I currently believe that the function is performing correctly.')
	disp('Now, we will do the comparison with Eqaulized Performance to find whether or not a line search')
	disp('over M1 will yield M1(affine) < M1(Luenberger) ')
	disp(' ')
	disp('Experiment 2')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Luenberger Line Search %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear n p d_u

	%Load a new system that CAN achieve Equalized Performance
	%load('data/system_examples/random1.mat')
	% load('data/system_examples/mIp1.mat')

	tricky_sys1.A = [0		1		0;
					 0 		0		1;
					 -0.5	-3.5	-0.25];
	tricky_sys1.B = [1;0;0];
	tricky_sys1.C = [1 0 0; 0 1 0];

	curr_sys = tricky_sys1;

	curr_sys.B = eye(size(curr_sys.A,1));

	%Create Constants for new system
	n = size(curr_sys.A,1);
	p = size(curr_sys.C,1);
	d_u = size(curr_sys.B,2);

	M_L = 1000;
	deltaM = M_L/2;
	finish_flag_L = false;
	while ~finish_flag_L

		%
		curr_sys.m = 0;
		curr_sys.d = 0.4;

		%Attempt to do the Luenberger Test
		clear delta mu
		delta 		= sdpvar(n,1,'full');
		mu 			= sdpvar(p,1,'full');
		curr_sys.x0 = sdpvar(n,1,'full');
		L 			= sdpvar(n,p,'full');

		%Create robust optimization Constraints
		robust_constrs = [];
		robust_constrs = robust_constrs + [-curr_sys.m <= mu <= curr_sys.m , uncertain(mu)];
		robust_constrs = robust_constrs + [-curr_sys.d <= delta <= curr_sys.d , uncertain(delta)];
		robust_constrs = robust_constrs + [ -M_L<= curr_sys.x0 <= M_L, uncertain(curr_sys.x0) ];

		%Define objective
		objective = [];
		objective = norm((curr_sys.A+L*curr_sys.C)*curr_sys.x0 + L*mu + delta,Inf);

		%Define Objective
		sol_attempt = optimize(robust_constrs+[objective <= M_L],objective , ops );

		% Check Termination Condition
		% +++++++++++++++++++++++++++
		if verbosity >= 0
			disp(['deltaM = ' num2str(deltaM) ])
		end

		if deltaM <= 0.01
			finish_flag_L = true;
		end

		% Adjust M_L
		% ++++++++++
		if sol_attempt.problem == 0		%If the problem was successfully solved, then decrease the proposed performance level.
			M_L = M_L - deltaM;
			deltaM = deltaM/2;
		else
			M_L = M_L + deltaM;
		end

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Finite Horizon Affine Estimator Line Search %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Define Constants
	T0 = T;

	M_A = 1000;
	deltaA = M_A/2;
	finish_flag_A = false;

	while ~finish_flag_A

		%
		curr_sys.m = 0;
		curr_sys.d = 0.4*M_A;

		%Attempt to do the FHAE Synthesis
		[~,sol_attempt] = achieve_eq_recovery_for(curr_sys,T0,M_A,M_A,verbosity);
		
		disp(sol_attempt.sol.info)

		% Check Termination Condition
		% +++++++++++++++++++++++++++
		if verbosity >= 0
			disp(['deltaA = ' num2str(deltaA) ])
		end

		if deltaA <= 0.01
			finish_flag_A = true;
		end

		% Adjust M_L
		% ++++++++++
		if sol_attempt.sol.problem == 0		%If the problem was successfully solved, then decrease the proposed performance level.
			M_A = M_A - deltaA;
			deltaA = deltaA/2;
		else
			M_A = M_A + deltaA;
		end

	end


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %% Run Controller with Synthesized Noise %%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% num_rollouts = 10^5;

	% %Save Results of Previous Optimizations
	% u_ol = value(u_open); 	%Open Loop Input during missing data period
	
	% Q_cl = value(Q);
	% r_cl = value(r);
	% F_cl = value( (inv(value(eye(size(Q,1)) + Q*Cm*S)) ) * Q);
	% %u0_cl = value((eye(size(F_cl,1)) + F_cl*Cm*S) * r);
	% u0_cl = value( inv( value(eye(size(Q,1)) + Q*Cm*S)) * r );

	% %Create noise
	% mu_t = unifrnd(-acc_e.m,acc_e.m,p*(T_missing+T_available),num_rollouts);
	% delta_t = unifrnd(-acc_e.d,acc_e.d,n*(T_missing+T_available),num_rollouts);
	% x0_t = unifrnd(-perf_level,perf_level,n,num_rollouts);

	% %Create Constants
	% [S1,~,~,~] = create_skaf_n_boyd_matrices(acc_e,T_missing);
	% [S2,H2,Cm2,~] = create_skaf_n_boyd_matrices(acc_e,T_available);

	% %Create Trajectories
	% diagA_str = '[';
	% for i = 0:T_missing
	% 	diagA_str = [diagA_str 'acc_e.A^' num2str(i) ];
	% 	if i ~= (T_missing)
	% 		diagA_str = [diagA_str ';'];
	% 	end
	% end
	% diagA_str = [ diagA_str ']'];
	% A_col = eval(diagA_str);
	% xi1 = A_col*x0_t + S1*(u_ol+delta_t(1:n*T_missing,1));

	% % Available data part
	% xi2(1:n,:) = xi1(end-n+1:end,:);
	% for t = T_missing:(T_missing+T_available-1)

	% 	%Create C diag, matrix
	% 	diagC_str = [];
	% 	diagC_str = '(';
	% 	for i = 1:t-T_missing+1
	% 		diagC_str = [ diagC_str 'acc_e.C' ];
	% 		if i ~= (t-T_missing+1)
	% 			diagC_str = [ diagC_str ',' ];
	% 		end
	% 	end
	% 	diagC_str = [ diagC_str ')' ];
	% 	diagC = eval(['blkdiag' diagC_str ]);

	% 	%Create Feedback Matrices
	% 	xi2_m1 		= xi2((t-T_missing)*n+1:(t-T_missing+1)*n,:);
	% 	F_cl_temp 	= F_cl((t-T_missing)*d_u+1:(t-T_missing+1)*d_u,1:(t-T_missing+1)*p);
	% 	mu_t_m1		= mu_t(T_missing*p+1:(t+1)*p,:);
	% 	u0_temp 	= u0_cl( (t-T_missing)*d_u+1:(t-T_missing+1)*d_u,1);

	% 	xi2((t-T_missing+1)*n+1:(t-T_missing+2)*n,:) = acc_e.A*xi2_m1 + F_cl_temp*( diagC * xi2(1:(t-T_missing+1)*n,:) + mu_t_m1) + repmat(u0_temp,1,num_rollouts) + delta_t( t*n+1:(t+1)*n,:);

	% 	% xi1_temp(t*)
	% end

	% % Find Norms
	% % ++++++++++

	% xi_t = [ xi1 ; xi2(n+1:end,:) ];

	% if (size(xi_t,1)/n) ~= (T_missing+T_available+1)
	% 	error('xi_t appears to be incorrectly sized.')
	% end

	% for test_num = 1:num_rollouts
	% 	for t = 0:T_missing+T_available

	% 		exp_norms(t+1,test_num) = norm( xi_t(t*n+1:(t+1)*n,test_num) , Inf );

	% 	end
	% end

	%%%%%%%%%%%%%%
	%% PLOTTING %%
	%%%%%%%%%%%%%%

	% figure;
	% hold on;
	% % bar([0:T_missing+T_available],value(alpha_l),'w');
	% bar([0:T_missing+T_available],[ perf_level 5*ones(1,T_missing+T_available-1) perf_level ],'w');
	% for test_num = 1:100
	% 	plot([0:T_missing+T_available],exp_norms(:,test_num));
	% end

	% xlabel('Time')
	% ylabel('\infty Norm of the Estimation Error')
	% title(['Estimator''s Error when data is missing until t = ' num2str(T_missing)])

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc;
	results.experim_params.T_missing = T_missing;
	results.experim_params.T_available = T_available;
	results.experim_params.perf_level = perf_level;
	% results.experim_params.num_rollouts = num_rollouts;

	results.func_comparison.full_optim = optim1;
	results.func_comparison.func_optim = optim2;
	results.func_comparison.error_bound_at_time_t = value(alpha_l);

	results.line_search.sys 			 = curr_sys;
	results.line_search.luenberger_min_M = M_L;

end