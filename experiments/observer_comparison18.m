function [ results ] = observer_comparison18( varargin )
%	observer_comparison18.m
%		Description:
%			The objective of this experiment is to design a finite horizon, affine
%			estimator (fhae) that is robust against the possibility of 1 observation
%			missing in the entire sequence of length T.
%
%			What makes this experiment different from 13, 14, and 16 is that:
%				- It uses the recent optimization presented by Prof. Yong.
%				- We make explicit use of an "E" matrix.
%				- The system chosen here is NOT ACC.
%
%			This Equalized Recovery Problem:
%				Let M1 T
%				Let ||xi(0)||<= M1. What is the minimum value for M1 such that,
%				- ||xi(t)||<= M2 \forall t \in [1,T-1]
%				- ||xi(T)||<= M1 
%				OR
%				- ||xi(t)||<= M1 \forall t
%
%		Inputs:
%			verbosity - 
%			T -			Time Horizon
%			M1 - 		


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
	elseif nargin == 3
		verbosity	= varargin{1};
		T	= varargin{2};
		M1 = varargin{3};
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
	end

	%Define our example system
	%a = [0.5 3.5 0.25];
	a = [0.5 0.4 0.25];

	example_sys.A =	[	0		1		0	;
						0		0		1	;
						-a(1)	-a(2)	-a(3)	];
	n = size(example_sys.A,1);

	example_sys.B = eye(n);
	example_sys.C = [ 0 1 0 ; 0 0 1 ];
	example_sys.E = [1 0 0]';%example_sys.B; 		%Noise occurs on the input signal

	example_sys.G = eye(size(example_sys.C,1));

	example_sys.m = 0; %0.5*perf_level;
	example_sys.d = 0.75*M1;

	test_sys = example_sys;

	n = size(test_sys.A,1);
	p = size(test_sys.C,1);

	m = size(test_sys.B,2);
	dd = size(test_sys.E,2); %Delta Dimension

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Announce New Experiments
	disp('===========================================')
	disp('Experiment 1: Synthesis For No Data Missing')
	disp(' ')
	disp('Assuming that M1 and T are given (they are constants), find the minimal M2.')
	disp(' ')

	% Constants
	% +++++++++

	E_big = [];
	for i = 1:T
		E_big = blkdiag(E_big,test_sys.E);
	end

	% Optimization Variables
	% ++++++++++++++++++++++
	delta 		= sdpvar(dd*T,1,'full');
	mu 			= sdpvar(p*T,1,'full');
	test_sys.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');

	alpha_l 	= sdpvar(T+1,1,'full');

	% Feedback Variables
	Q = sdpvar(m*T,p*T,'full');
	r = sdpvar(m*T,1,'full');

	% Dual Variables
	Pi_1 = sdpvar(2*n*T,2*(dd+p)*T+2*n,'full');
	Pi_2 = sdpvar(2*n,2*(dd+p)*T+2*n,'full');

	% Creating Constraints
	% ++++++++++++++++++++

	[S0,H0,Cm0,xi0m] = create_skaf_n_boyd_matrices(test_sys,T);

	positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	noise_constrs = [ Pi_1 * [ test_sys.d * ones(2*dd*T,1) ; test_sys.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
	noise_constrs = noise_constrs + [ Pi_2 * [ test_sys.d * ones(2*dd*T,1) ; test_sys.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; test_sys.A^i];
	end

	G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*E_big S0*Q (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

	bounded_disturb_matrix = [ [ eye(dd*T) ; -eye(dd*T) ] zeros(2*dd*T,p*T+n) ;
								zeros(2*p*T,dd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
								zeros(2*n,(p+dd)*T) [ eye(n) ; -eye(n) ] ];

	dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

	%Lower Diagonal Constraint
	l_diag_constr = [];
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(test_sys.B,2)+1:bl_row_num*size(test_sys.B,2)], ...
												[bl_row_num*size(test_sys.C,1)+1:end] ) == 0 ];
	end

	% OPTIMIZATION
	% ++++++++++++
	ops = sdpsettings('verbose',verbosity);
	optim1 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr, ...
			alpha_2, ...
			ops)

	if optim1.problem ~= 0
		error(['The design problem was not completely solved.' optim1.info ])
	end

	% Save Feedback Matrices
	% ++++++++++++++++++++++
	Q1 = value(Q);
	r1 = value(r);
	F1 = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	u0_1 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

	opt_obj1 = value(alpha_2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot Results for First Experiment Test %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	contr1.F = F1;
	contr1.u0 = u0_1;

	num_rollouts = 10^5;

	[xi_t,xi_mag_t] 	= apply_controller_to_rollouts(test_sys,contr1,T,num_rollouts,M1);

	rollouts_per_plot = 1000;

	figure;
	hold on;
	bar([0:T],[ M1 value(alpha_2)*ones(1,T-1) M1 ],'w')
	for i = 1:rollouts_per_plot
		plot([0:T],xi_mag_t(:,i))
	end

	xlabel('Time')
	ylabel('\infty Norm of the Estimation Error')
	legend('Guarantees')
	title('Estimator''s Error when ALL Data is available')
	xt = get(gca, 'XTick');
	set(gca, 'FontSize', 16)

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %% Synthesize Controller for 1 missing Data Case %%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	clear delta
	clear mu
	clear alpha_l

	%Announce new experiment
	disp('====================================================')
	disp('Experiment 2: Synthesis For "At Most 1" Data Missing')

	% Can we simply try to add further constraints?
	% +++++++++++++++++++++++++++++++++++++++++++++

	for missing_loc = 0:(T-1)
		% Create Trajectory Matrices
		[~,~,Cm,~] = create_skaf_n_boyd_matrices(test_sys,T,'missing',missing_loc);

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
			ops)

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

	%Controller Definition
	contr2.F = F2;
	contr2.u0 = u0_2;

	num_plots = 4;

	%Plot Constants
	[xi_t0,xi_mag_t0] = apply_controller_to_rollouts(test_sys,contr2,T,num_rollouts,M1);
	xi{1} = xi_t0;
	xi_mag{1} = xi_mag_t0;

	figure;
	for i = 1 : num_plots
		if i > 1
			[xi{i},xi_mag{i}] = apply_controller_to_rollouts(test_sys,contr2,T,num_rollouts,M1,'missing',i-1);
		end

		subplot(2,2,i)
		hold on;

		bar([0:T],[ M1 value(alpha_2)*ones(1,T-1) M1 ],'w')
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
		xt = get(gca, 'XTick');
		set(gca, 'FontSize', 16)
	end


	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = test_sys;
	results.experim_params.T = T;
	results.experim_params.M1 = M1;

	results.none_missing.optimization = optim1;
	results.none_missing.opt_obj = opt_obj1;
	results.none_missing.Q = Q1;
	results.none_missing.r = r1;
	results.none_missing.F = F1;
	results.none_missing.u0 = u0_1;

	results.missing1.optimization = optim2;
	results.missing1.opt_obj = opt_obj2;
	results.missing1.Q = Q2;
	results.missing1.r = r2;
	results.missing1.F = F2;
	results.missing1.u0 = u0_2;

end