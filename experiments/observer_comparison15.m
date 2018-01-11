function [ results ] = observer_comparison15( varargin )
%	observer_comparison15.m
%		Description:
%			The objective of this experiment is to design a finite horizon, affine
%			estimator (fhae) that is robust against the possibility of 1 observation
%			missing in the entire sequence of length T. THIS USES THE ROBUST OPTIMIZATION
%			FRAMEWORK PROPOSED BY SZE ZHENG YONG. The last sentence explains how this
%			differs from experiment set #13.
%
%			Qualities of this Problem:
%			- Robust Optimization proposed by Sze Zheng Yong
%			- Minimizing mu_1 and mu_2.
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
%			verbosity - 
%			T -			Time Horizon

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
	else
		error('Unacceptable number of arguments.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%Default Values for observer_comparison
	if nargin < 2
		T = 6;
	end

	%default Value for perf_level
	% if ~exist('perf_level')
	% 	perf_level = 1;
	% end

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	n = size(acc.A,1);
	p = size(acc.C,1);

	%Create Error System
	acc_e = acc;
	acc_e.B = eye(n);

	wd = size(acc_e.E,2);
	vd = size(acc.C,1);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Announce New Experiments
	disp('=================================================')
	disp('Experiment 1: Observing the Effect of Changing M1')
	disp(' ')

	% Constants
	% +++++++++

	% Optimization Variables
	% ++++++++++++++++++++++
	w 			= sdpvar(wd*T,1,'full');
	v 			= sdpvar(vd*T,1,'full');
	acc_e.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');
	M1_list		= [ 0.6 0.8 1 10 ];

	alpha_l 	= sdpvar(T+1,1,'full');

	% Feedback Variables
	Q = sdpvar(n*T,p*T,'full');
	r = sdpvar(n*T,1,'full');

	% Dual Variables
	Pi_1 = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
	Pi_2 = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');

	for M1_ind = 1 : length(M1_list)

		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,E_big] = create_skaf_n_boyd_matrices(acc_e,T);

		positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
		end

		noise_constrs = [ Pi_1 * [ acc_e.d * ones(2*wd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1_list(M1_ind) * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
		noise_constrs = noise_constrs + [ Pi_2 * [ acc_e.d * ones(2*wd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1_list(M1_ind) * ones(2*n,1) ] <= M1_list(M1_ind) * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T
			pre_xi = [ pre_xi ; acc_e.A^i];
		end

		G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*E_big S0*Q (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,vd*T+n) ;
									zeros(2*vd*T,wd*T) [ eye(vd*T) ; -eye(vd*T) ] zeros(2*vd*T,n) ;
									zeros(2*n,(vd+wd)*T) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

		%Lower Diagonal Constraint
		l_diag_constr = [];
		for bl_row_num = 1 : T-1
			l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
													[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
		end

		% OPTIMIZATION
		% ++++++++++++
		ops = sdpsettings('verbose',verbosity);
		optim1{M1_ind} = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr, ...
									alpha_2, ...
									ops);

		if optim1{M1_ind}.problem ~= 0
			error(['The design problem was not completely solved.' optim1.info ])
		end

		% Save Feedback Matrices
		% ++++++++++++++++++++++
		Q1{M1_ind} = value(Q);
		r1{M1_ind} = value(r);
		F1{M1_ind} = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
		u0_1{M1_ind} 		= value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

		opt_obj1{M1_ind} 	= value(alpha_2);

		disp(['Synthesized Controller for M1 = ' num2str(M1_list(M1_ind)) '. Found minimum M2 = ' num2str(value(alpha_2)) ])

	end

	% Save Results
	results.m1_changed.m1 = M1_list;


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot Results for First Experiment Test %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	rollouts_per_plot = 1000;
	xi_t1 = {}; xi_mag_t1 = {};

	figure;

	for controller_num = 1 : length(M1_list)

		contr1.F = F1{controller_num};
		contr1.u0 = u0_1{controller_num};

		num_rollouts = 10^4;

		[xi_t1{controller_num},xi_mag_t1{controller_num}] 	= apply_controller_to_rollouts(acc_e,contr1,T,num_rollouts,M1_list(controller_num));

		subplot(2,2,controller_num);
		hold on;
		bar([0:T],[ M1_list(controller_num) opt_obj1{controller_num}*ones(1,T-1) M1_list(controller_num) ],'w')
		for i = 1:rollouts_per_plot
			plot([0:T],xi_mag_t1{controller_num}(:,i))
		end

		xlabel('Time')
		ylabel('\infty Norm of the Estimation Error')
		title(['Estimator Performance when $M_1$=' num2str(M1_list(controller_num))],'Interpreter','latex')

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller for 1 missing Data Case %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	clear delta
	clear mu
	clear alpha_l

	% clear M1_list
	% clear M1_ind

	%Announce new experiment
	disp('====================================================')
	disp('Experiment 2: Synthesis For "At Most 1" Data Missing')
	disp('              Applied ''continuously'' on a trajectory')
	disp(' ')

	% Constants
	% +++++++++

	M1 = 1;
	perf_level = M1;

	% Creating Constraints
	% ++++++++++++++++++++

	[S0,H0,Cm0,xi0m,E_big] = create_skaf_n_boyd_matrices(acc_e,T);

	positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	noise_constrs = [ Pi_1 * [ acc_e.d * ones(2*wd*T,1) ; acc_e.m * ones(2*p*T,1) ; perf_level * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
	noise_constrs = noise_constrs + [ Pi_2 * [ acc_e.d * ones(2*wd*T,1) ; acc_e.m * ones(2*p*T,1) ; perf_level * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; acc_e.A^i];
	end

	G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*E_big S0*Q (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

	bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,vd*T+n) ;
								zeros(2*vd*T,wd*T) [ eye(vd*T) ; -eye(vd*T) ] zeros(2*vd*T,n) ;
								zeros(2*n,(vd+wd)*T) [ eye(n) ; -eye(n) ] ];

	dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

	%Lower Diagonal Constraint
	l_diag_constr = [];
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
												[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
	end

	% Can we simply try to add further constraints?
	% +++++++++++++++++++++++++++++++++++++++++++++

	for missing_loc = 1:(T-1)-1
		% Create Trajectory Matrices
		[~,~,Cm,~] = create_skaf_n_boyd_matrices(acc_e,T,'missing',missing_loc);

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

		% Pxd = [ Pxd ; (eye(n*(T+1))+S*Q*Cm)*S ];%E_bar*[ eye(b_dim*t_horizon) zeros(b_dim*t_horizon,b_dim) ];
		% Pxm = [ Pxm ; S*Q*mu_select ];
		% xi_tilde = [ xi_tilde ; (eye(n*(T+1)) + S*Q*Cm)*xi0m + S*r];

		G = [ (eye(n*(T+1))+S0*Q*Cm)*S0*E_big S0*Q*mu_select (eye(n*(T+1))+S0*Q*Cm)*pre_xi ];

		%Add to the constraint set
		dual_equal_constrs = dual_equal_constrs + [Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G];
		dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim2 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr, ...
			alpha_2, ...
			ops);

	if optim2.problem ~= 0
		error(['The design problem was not completely solved.' optim2.info ])
	end

	% Save Feedback Matrices
	% ++++++++++++++++++++++
	Q2 = value(Q);
	r2 = value(r);
	F2 = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	u0_2 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

	opt_obj2 = value(alpha_2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot the Performance of this Designed Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Constants
	T_s = 0.5;

	%Controller Definition
	contr2.F = F2;
	contr2.u0 = u0_2;

	num_plots = 4;
	num_rollouts = 10;

	A_col = [];
	for  temp_ind = 0:T
		A_col = [ A_col ; acc_e.A^temp_ind ];
	end 

	E_bar = [];
	for temp_ind = 1:T
		E_bar = blkdiag(E_bar,acc_e.E);
	end

	%Plot Constants
	[xi_t0,xi_mag_t0] = apply_controller_to_rollouts(acc_e,contr2,T,num_rollouts,M1);
	xi = xi_t0;
	xi_mag = xi_mag_t0;

	for missing_ob_num = 1:3
		[S,H,Cm,~] = create_skaf_n_boyd_matrices(acc_e,T,'missing',missing_ob_num);

		% Create some trajectory matrices
		Pxd = (eye(n*(T+1))+H*contr2.F*inv(eye(p*T)-Cm*H*contr2.F)*Cm)*S*E_big;
		Pxm = H*contr2.F*inv(eye(p*T)-Cm*H*contr2.F);
		xi_factor = H*contr2.F*inv(eye(p*T)-Cm*H*contr2.F )*Cm;

		delta 	= unifrnd(-acc_e.d,acc_e.d,wd*T,num_rollouts);
		mu 		= unifrnd(-acc_e.m,acc_e.m,vd*T,num_rollouts);
		xi_0 	= xi(end-n+1:end,:);

		mu([missing_ob_num*p+1:(missing_ob_num+1)*p],:) = 0;

		xi_temp =  A_col*xi_0 + H*repmat(contr2.u0,1,num_rollouts) + ...
	        xi_factor*(A_col*xi_0 + H*repmat(contr2.u0,1,num_rollouts)) + ...
	        Pxd * delta + Pxm * mu;

	    xi = [xi;xi_temp([n+1:end],:)];
	    for rollout_ind = 1 : num_rollouts 
			for t = 1:T
				xi_mag(1+missing_ob_num*T+t,rollout_ind) = norm(select_m(t,T)*xi_temp(:,rollout_ind),Inf);
			end
		end
	end

	figure;
	hold on;
	bar([0:T*(3+1)]*T_s,[ M1 opt_obj2*ones(1,T-1) M1 opt_obj2*ones(1,T-1) M1 opt_obj2*ones(1,T-1) M1 opt_obj2*ones(1,T-1) M1 ],'w')
	for rollout_ind = 1 : num_rollouts
		plot([0:T*(3+1)]*T_s,xi_mag(:,rollout_ind)' )
	end

	legend('Guarantees')
	axis([-0.5*T_s (T*(3+1)+0.5)*T_s 0 3.5])

	xlabel('Time [sec]')
	ylabel('$||x(t)-\hat{x}(t)||_{\infty}$','Interpreter','latex')
	% title(['Estimator Performance when $M_1$=' num2str(M1_list(controller_num))],'Interpreter','latex')

	xt = get(gca, 'XTick');
	set(gca, 'FontSize', 16)

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc;
	results.experim_params.T = T;
	results.experim_params.M1s =M1_list;

	results.m1_changed.optimization = optim1;
	results.m1_changed.opt_obj = opt_obj1;
	results.m1_changed.Q = Q1;
	results.m1_changed.r = r1;
	results.m1_changed.F = F1;
	results.m1_changed.u0 = u0_1;
	results.m1_changed.xi = xi_t1;
	results.m1_changed.xi_mag = xi_mag_t1;

	results.missing1.optimization = optim2;
	results.missing1.opt_obj = opt_obj2;
	results.missing1.Q = Q2;
	results.missing1.r = r2;
	results.missing1.F = F2;
	results.missing1.u0 = u0_2;

end