function [ results ] = observer_comparison11( varargin )
%	observer_comparison11.m
%		The objective of this experiment is to synthesize a design that achieves
%		Equalized Recovery for a given T_m and T_a on the ACC (or other) system
%		and then test the controller using noise that follows the constraints given
%		in the problem.

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
	delta 		= sdpvar(n*(T_missing+T_available),1,'full');
	u_open 		= sdpvar(d_u*T_missing,1,'full');
	acc_e.x0 	= sdpvar(n,1,'full');

	alpha1		= sdpvar(1,1,'full');

	alpha_l 	= sdpvar(T_missing+T_available+1,1,'full');

	% Missing Observations for T_missing time steps
	% +++++++++++++++++++++++++++++++++++++++++++++

	% Create Trajectory Matrices
	[S,~,~,xi0m] = create_skaf_n_boyd_matrices(acc_e,T_missing);

	% Create Objective + Constraints
	obj = norm( select_m(T_missing,T_missing)*(xi0m + S*(u_open + delta(1:n*T_missing,1)) ), Inf );
	epi_constr1 = [ obj <= alpha1 ];

	robust_constrs = [ -acc_e.d <= delta <= acc_e.d , uncertain(delta)] + [ -perf_level <= acc_e.x0 <= perf_level , uncertain(acc_e.x0) ];

	graph_constrs = [];
	for i = 0:T_missing
		graph_constrs = graph_constrs + [ norm(select_m(i,T_missing)*(xi0m + S*(u_open + delta(1:n*T_missing,1)) ),Inf) <= alpha_l(i+1) ];
	end

	%Optimize
	ops = sdpsettings('verbose',verbosity);
	optim1 = optimize(epi_constr1+robust_constrs+graph_constrs,alpha1 + sum(alpha_l(1:T_missing+1)),ops);

	% Availaible Observations for T_available time steps
	% ++++++++++++++++++++++++++++++++++++++++++++++++++

	% Variables
	mu 			= sdpvar(p*(T_missing+T_available),1,'full');

	alpha2		= sdpvar(1,1,'full');
	alpha3 		= sdpvar(1,1,'full');

	% Create Trajectory Matrices
	[S,H,Cm,xi0m] = create_skaf_n_boyd_matrices(acc_e,T_available);

	% Create Objective (Epigraph constraint)

	Q = sdpvar(size(H,2),size(Cm,1),'full');
	r = sdpvar(size(H,2),1,'full');

	Pxd = (eye(n*(T_available+1))+S*Q*Cm)*S ;%E_bar*[ eye(b_dim*t_horizon) zeros(b_dim*t_horizon,b_dim) ];
	Pxm = S*Q;
	xi_tilde = (eye(n*(T_available+1)) + S*Q*Cm)*xi0m + S*r;

	R = [ zeros(n,n*T_available) eye(n) ];

	objective = norm( R*(xi_tilde + Pxd * delta(n*T_missing+1:end,1) + Pxm * mu(p*T_missing+1:end,1)) , Inf );
	interm_norms = norm([ eye(n*T_available) zeros(n*T_available,n) ]*( xi_tilde + Pxd * delta(n*T_missing+1:end,1) + Pxm * mu(p*T_missing+1:end,1) ) , Inf);
	interm_vals = [ eye(n*T_available) zeros(n*T_available,n) ]*( xi_tilde + Pxd * delta(n*T_missing+1:end,1) + Pxm * mu(p*T_missing+1:end,1) );

	epi_constr2 = [ objective <= alpha2 , interm_norms <= alpha3 ];
	% epi_constr2 = [ objective <= alpha2, interm_norms <= 20 ];
	% epi_constr2 = [ objective <= alpha2];

	% Feasibility Constraints
	feasib_constrs = [];
	feasib_constrs = feasib_constrs + [alpha2 <= perf_level];

	% Create Robustification Constraints

	robust_constrs = [];
	robust_constrs = robust_constrs + [ -acc_e.d <= delta <= acc_e.d , uncertain(delta) ];
	robust_constrs = robust_constrs + [ -acc_e.m <= mu <= acc_e.m , uncertain(mu)];
	robust_constrs = robust_constrs + [ -value(alpha1) <= acc_e.x0 <= value(alpha1) , uncertain(acc_e.x0) ];

	% Create Causality (Lower Diagonal) Constraint
	l_diag_constr = [];
	for bl_row_num = 1 : T_available-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
												[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
	end

	% Create Graph Constraints
	graph_constrs = [];
	for i = T_missing+1:T_missing+T_available
		graph_constrs = graph_constrs + [ norm(select_m(i-T_missing,T_available)*( xi_tilde + Pxd * delta(n*T_missing+1:end,1) + Pxm * mu(p*T_missing+1:end,1) ) , Inf) <= alpha_l(i+1) ];
	end

	% Optimize!
	optim2 = optimize(epi_constr2+robust_constrs+l_diag_constr+feasib_constrs+graph_constrs ,alpha2+alpha3+sum(alpha_l(T_missing+2:T_missing+T_available+1)),ops);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Run Controller with Synthesized Noise %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	num_rollouts = 10^5;

	%Save Results of Previous Optimizations
	u_ol = value(u_open); 	%Open Loop Input during missing data period
	
	Q_cl = value(Q);
	r_cl = value(r);
	F_cl = value( (inv(value(eye(size(Q,1)) + Q*Cm*S)) ) * Q);
	%u0_cl = value((eye(size(F_cl,1)) + F_cl*Cm*S) * r);
	u0_cl = value( inv( value(eye(size(Q,1)) + Q*Cm*S)) * r );

	%Create noise
	mu_t = unifrnd(-acc_e.m,acc_e.m,p*(T_missing+T_available),num_rollouts);
	delta_t = unifrnd(-acc_e.d,acc_e.d,n*(T_missing+T_available),num_rollouts);
	x0_t = unifrnd(-perf_level,perf_level,n,num_rollouts);

	%Create Constants
	[S1,~,~,~] = create_skaf_n_boyd_matrices(acc_e,T_missing);
	[S2,H2,Cm2,~] = create_skaf_n_boyd_matrices(acc_e,T_available);

	%Create Trajectories
	diagA_str = '[';
	for i = 0:T_missing
		diagA_str = [diagA_str 'acc_e.A^' num2str(i) ];
		if i ~= (T_missing)
			diagA_str = [diagA_str ';'];
		end
	end
	diagA_str = [ diagA_str ']'];
	A_col = eval(diagA_str);
	xi1 = A_col*x0_t + S1*(u_ol+delta_t(1:n*T_missing,1));

	% Available data part
	xi2(1:n,:) = xi1(end-n+1:end,:);
	for t = T_missing:(T_missing+T_available-1)

		%Create C diag, matrix
		diagC_str = [];
		diagC_str = '(';
		for i = 1:t-T_missing+1
			diagC_str = [ diagC_str 'acc_e.C' ];
			if i ~= (t-T_missing+1)
				diagC_str = [ diagC_str ',' ];
			end
		end
		diagC_str = [ diagC_str ')' ];
		diagC = eval(['blkdiag' diagC_str ]);

		%Create Feedback Matrices
		xi2_m1 		= xi2((t-T_missing)*n+1:(t-T_missing+1)*n,:);
		F_cl_temp 	= F_cl((t-T_missing)*d_u+1:(t-T_missing+1)*d_u,1:(t-T_missing+1)*p);
		mu_t_m1		= mu_t(T_missing*p+1:(t+1)*p,:);
		u0_temp 	= u0_cl( (t-T_missing)*d_u+1:(t-T_missing+1)*d_u,1);

		xi2((t-T_missing+1)*n+1:(t-T_missing+2)*n,:) = acc_e.A*xi2_m1 + F_cl_temp*( diagC * xi2(1:(t-T_missing+1)*n,:) + mu_t_m1) + repmat(u0_temp,1,num_rollouts) + delta_t( t*n+1:(t+1)*n,:);

		% xi1_temp(t*)
	end

	% Find Norms
	% ++++++++++

	xi_t = [ xi1 ; xi2(n+1:end,:) ];

	if (size(xi_t,1)/n) ~= (T_missing+T_available+1)
		error('xi_t appears to be incorrectly sized.')
	end

	for test_num = 1:num_rollouts
		for t = 0:T_missing+T_available

			exp_norms(t+1,test_num) = norm( xi_t(t*n+1:(t+1)*n,test_num) , Inf );

		end
	end

	%%%%%%%%%%%%%%
	%% PLOTTING %%
	%%%%%%%%%%%%%%

	figure;
	hold on;
	% bar([0:T_missing+T_available],value(alpha_l),'w');
	bar([0:T_missing+T_available],[ perf_level 5*ones(1,T_missing+T_available-1) perf_level ],'w');
	for test_num = 1:100
		plot([0:T_missing+T_available],exp_norms(:,test_num));
	end

	xlabel('Time')
	ylabel('\infty Norm of the Estimation Error')
	title(['Estimator''s Error when data is missing until t = ' num2str(T_missing)])

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc;
	results.experim_params.T_missing = T_missing;
	results.experim_params.T_available = T_available;
	results.experim_params.perf_level = perf_level;
	results.experim_params.num_rollouts = num_rollouts;

	results.synth_data.xi_t = xi_t;
	results.synth_data.mu_t = mu_t;
	results.synth_data.delta_t = delta_t;
	results.synth_data.xi0_t = xi0m

	results.missing_optim = optim1;
	results.missing_optim.opt_obj = value(alpha1);
	results.missing_optim.u_open = value(u_open);

	results.available_optim = optim2;
	results.available_optim.opt_final_error = value(alpha2);
	results.available_optim.opt_interm_error = value(alpha3);
	results.available_optim.Q = value(Q);
	results.available_optim.r = value(r);
	results.available_optim.F = F_cl;
	results.available_optim.u0 = u0_cl;
	results.available_optim.error_bound_at_time_t = value(alpha_l);

end