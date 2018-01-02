function [ results ] = observer_comparison13( varargin )
%	observer_comparison13.m
%		Description:
%			The objective of this experiment is to design a finite horizon, affine
%			estimator (fhae) that is robust against the possibility of 1 observation
%			missing in the entire sequence of length T.
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
	elseif nargin == 2
		verbosity	= varargin{1};
		T 			= varargin{2};
	elseif nargin == 4
		verbosity	= varargin{1};
		T	= varargin{2};
		M1 = varargin{3};
		M2 	= varargin{4};
	else
		error('Unacceptable number of arguments.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%Default Values for observer_comparison
	if nargin < 2
		T = 6;
		M1 = 1;
		M2 = 4;
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

	% Announce New Experiments
	disp('Experiment 1: Synthesis Attempts')

	% Variables
	delta 		= sdpvar(n*T,1,'full');
	mu 			= sdpvar(p*T,1,'full');
	acc_e.x0 	= sdpvar(n,1,'full');

	pl_t		= sdpvar(1,1,'full');

	alpha_l 	= sdpvar(T+1,1,'full');

	% Create Objective (Epigraph constraint)
	Q = sdpvar(d_u*T,p*T,'full');
	r = sdpvar(d_u*T,1,'full');

	% Creating Constraints assuming that 1 piece of data is missing
	% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	epi_constrs = []; graph_constrs = [];

	Pxd = []; Pxm = []; xi_tilde = [];

	for missing_loc = 0:(T-1)-1
		% Create Trajectory Matrices
		[S,H,Cm,xi0m] = create_skaf_n_boyd_matrices(acc_e,T,'missing',missing_loc);

		%Create Special selection matrices for selecting the proper variables
		R = [ zeros(n,n*T) eye(n) ];
		mu_select = [];
		for i = 1:T
			if any(missing_loc == (i-1))
				mu_select = [ mu_select ; zeros(p,p*T) ];
			else
				mu_select = [ mu_select ; [ zeros(p,p*(i-1)) eye(p) zeros(p,p*(T-i)) ] ];
			end
		end

		Pxd = [ Pxd ; (eye(n*(T+1))+S*Q*Cm)*S ];%E_bar*[ eye(b_dim*t_horizon) zeros(b_dim*t_horizon,b_dim) ];
		Pxm = [ Pxm ; S*Q*mu_select ];
		xi_tilde = [ xi_tilde ; (eye(n*(T+1)) + S*Q*Cm)*xi0m + S*r];

		% Formalize Trajectories for objective/intermediate level
		% objective_term = [ objective_term ; R*(xi_tilde{missing_loc+1} + Pxd{missing_loc+1} * delta + Pxm{missing_loc+1} * mu_select * mu) ];
		% interm_norms{missing_loc+1} = norm([ eye(n*T) zeros(n*T,n) ]*( xi_tilde{missing_loc+1} + Pxd{missing_loc+1} * delta + Pxm{missing_loc+1} * mu_select * mu ) , Inf);
		% interm_vals = [ eye(n*T_available) zeros(n*T_available,n) ]*( xi_tilde + Pxd * delta(n*T_missing+1:end,1) + Pxm * mu(p*T_missing+1:end,1) );

		% epi_constrs = epi_constrs + [ objective{missing_loc+1} <= M1 , interm_norms{missing_loc+1} <= M2 ];

		% When considering the 
		% for t = 0 : T
		% 	graph_constrs = graph_constrs + [ norm(select_m(t,T)*( xi_tilde{missing_loc+1} + Pxd{missing_loc+1} * delta + Pxm{missing_loc+1} * mu_select * mu ) , Inf) <= alpha_l(t+1,missing_loc+1) ];
		% 	% disp(t)
		% end
	end

	% Create Trajectory Matrices
	[S0,H0,Cm0,xi0m] = create_skaf_n_boyd_matrices(acc_e,T);

	Pxd0 = (eye(n*(T+1))+S0*Q*Cm0)*S0 ;%E_bar*[ eye(b_dim*t_horizon) zeros(b_dim*t_horizon,b_dim) ];
	Pxm0 = S0*Q;
	xi_tilde0 = (eye(n*(T+1)) + S0*Q*Cm0)*xi0m + S0*r;

	%
	Pxd = [ Pxd ; Pxd0 ];
	Pxm = [ Pxm ; Pxm0 ];
	xi_tilde = [xi_tilde ; xi_tilde0];

	%Create selection matrices
	R = [];
	R2 = [];
	for s_ind = 1:T
		R = [ R ; [ zeros(n,n*(T+1)*(s_ind-1)) select_m(T,T) zeros(n,n*(T+1)*(T-s_ind)) ] ];
		R2 = [R2 ; [ zeros(n*(T-1),n*(T+1)*(s_ind-1)) zeros(n*(T-1),n) eye(n*(T-1)) zeros(n*(T-1),n) zeros(n*(T-1),n*(T+1)*(T-s_ind)) ] ];
	end

	objective = norm( R*(xi_tilde + Pxd * delta + Pxm * mu) , Inf );
	interm_norms = norm(R2*( xi_tilde + Pxd * delta + Pxm * mu ) , Inf);

	epi_constrs = epi_constrs + [objective <= M1 , interm_norms <= M2];

	% Create Graph Constraints
	% graph_constrs = [];
	comm_time_ind = {};
	for time_ind = 0:T
		temp_select_mat = [];
		for s_ind = 1:T
			temp_select_mat = [ temp_select_mat ; [ zeros(n,n*(T+1)*(s_ind-1)) select_m(time_ind,T) zeros(n,n*(T+1)*(T-s_ind)) ] ];
			comm_time_ind{time_ind+1} = temp_select_mat;
		end
	end

	for i = 0:T
		graph_constrs = graph_constrs + [ norm(comm_time_ind{i+1}*( xi_tilde + Pxd * delta + Pxm * mu ) , Inf) <= alpha_l(i+1) ];
	end

	% Feasibility Constraints
	% +++++++++++++++++++++++

	feasib_constrs = [];
	feasib_constrs = feasib_constrs + [alpha_l <= M2];

	% Create Robustification Constraints
	% ++++++++++++++++++++++++++++++++++

	robust_constrs = [];
	robust_constrs = robust_constrs + [ -acc_e.d <= delta <= acc_e.d , uncertain(delta) ];
	robust_constrs = robust_constrs + [ -acc_e.m <= mu <= acc_e.m , uncertain(mu)];
	robust_constrs = robust_constrs + [ -M1 <= acc_e.x0 <= M1 , uncertain(acc_e.x0) ];

	% Create Causality (Lower Diagonal) Constraint
	% ++++++++++++++++++++++++++++++++++++++++++++

	l_diag_constr = [];
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
												[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
	end

	% Optimize!
	ops = sdpsettings('verbose',verbosity);
	optim1 = optimize(epi_constrs+robust_constrs+l_diag_constr+feasib_constrs+graph_constrs ,sum(alpha_l),ops);

	if optim1.problem ~= 0
		warning('The design problem was not completely solved.')
	end

	% Save Results
	Q_cl = value(Q);
	r_cl = value(r);

	F_cl 	= value( (inv(value(eye(size(Q,1)) + Q*Cm*S)) ) * Q);
	u0_cl 	= value( inv(value(eye(size(Q,1)) + Q*Cm*S)) * r );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Apply Controller to Synthesized Data %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	contr1.F = F_cl;
	contr1.u0 = u0_cl;

	num_rollouts = 10^4;

	[xi_t,xi_mag_t] 	= apply_controller_to_rollouts(acc_e,contr1,T,num_rollouts,M1);
	[xi_t2,xi_mag_t2] 	= apply_controller_to_rollouts(acc_e,contr1,T,num_rollouts,M1,'missing',3);
	[xi_t3,xi_mag_t3] 	= apply_controller_to_rollouts(acc_e,contr1,T,num_rollouts,M1,'missing',4);

	%%%%%%%%%%%
	%% Plots %%
	%%%%%%%%%%%

	rollouts_per_plot = 1000;
	m_bounds = [M1 ones(1,T-1)*M2 M1];

	num_plots = 4;

	xi_mag{1} = xi_mag_t; 

	figure;
	for i = 1 : num_plots
		if i > 1
			[xi{i},xi_mag{i}] = apply_controller_to_rollouts(acc_e,contr1,T,num_rollouts,M1,'missing',i-1);
		end

		subplot(2,2,i)
		hold on;

		errorbar([0:T],m_bounds/2,m_bounds/2)
		for r_num = 1:rollouts_per_plot
			plot([0:T],xi_mag{i}(:,r_num))
		end
		xlabel('Time Step t')
		ylabel('$||\xi(t)||\infty$, Norm of Estimation Error','Interpreter','latex')
		if i > 1
			title(['Norm of the Estimation Error (Missing Data Occurs at t=' num2str(i-1) ')'])
		else
			title(['Norm of the Estimation Error (Data Always Available)'])
		end
		legend('Guarantees')
	end

	% subplot(1,3,1)
	% hold on;
	
	% errorbar([0:T],alpha_l/2,alpha_l/2)
	% for i = 1:rollouts_per_plot
	% 	plot([0:T],xi_mag_t(:,i))
	% end

	% subplot(1,3,2)
	% hold on;
	% errorbar([0:T],alpha_l/2,alpha_l/2)
	% for i = 1:rollouts_per_plot
	% 	plot([0:T],xi_mag_t2(:,i))
	% end	

	% subplot(1,3,3)
	% hold on;
	% errorbar([0:T],alpha_l/2,alpha_l/2)
	% for i = 1:rollouts_per_plot
	% 	plot([0:T],xi_mag_t3(:,i))
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
	results.experim_params.T = T;
	results.experim_params.M1 = M1;

	results.robustifying.full_optim = optim1;
	results.robustifying.error_bound_at_time_t = value(alpha_l);
	results.robustifying.Q = Q_cl;
	results.robustifying.r = r_cl;
	results.robustifying.F = F_cl;
	results.robustifying.u0 = u0_cl;

	results.rollouts.xi1 = xi_t;
	results.rollouts.xi1_mag = xi_mag_t;
	results.rollouts.xi2 = xi_t2;
	results.rollouts.xi2_mag = xi_mag_t2;

end