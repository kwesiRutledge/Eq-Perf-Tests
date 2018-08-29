function [ results ] = observer_comparison25( varargin )
%	observer_comparison25.m
%		Description:
%			The objective of this experiment is to design a finite horizon, affine
%			estimator (fhae) that is robust against the possibility of 1 observation
%			missing in the entire sequence of length T.
%			We will compare: Worst Case Language method v. My proposal
%
%			What makes this experiment different from 13 and 14 is that:
%				- It uses the recent optimization presented by Prof. Yong.
%				- We make explicit use of the "E" matrix in the ACC example.
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
	T_s = 0.5;

	%Default Values for observer_comparison
	if nargin < 2
		T = 6;
		M1 = 1;
	end

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	n = size(acc.A,1);
	p = size(acc.C,1);

	%Create Error System
	acc_e = acc;
	acc_e.B = eye(n);

	m = size(acc_e.B,2);
	dd = size(acc_e.E,2); %Delta Dimension

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %% Synthesize Controller for 1 missing Data Case %%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	clear delta
	clear mu
	clear alpha_l

	%Announce new experiment
	disp('=============================================================================')
	disp('Experiment 1: Synthesis For "At Most 1" Data Missing with worst case Language')

	L = [];
	L = [	1,0,1,1,1,1;
			1,1,0,1,1,1	];

	L_star = ones(1,6);
	for word_in_L = 1:size(L,1)
		L_star = bitand(L_star,L(word_in_L,:));
	end

	% Optimization Variables
	% ++++++++++++++++++++++
	delta 		= sdpvar(dd*T,1,'full');
	mu 			= sdpvar(p*T,1,'full');
	acc_e.x0 	= sdpvar(n,1,'full');

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

	[S0,H0,Cm0,xi0m,E_big] = create_skaf_n_boyd_matrices(acc_e,T,'missing',find(L_star == 0)-1);

	positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	noise_constrs = [ Pi_1 * [ acc_e.d * ones(2*dd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
	noise_constrs = noise_constrs + [ Pi_2 * [ acc_e.d * ones(2*dd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; acc_e.A^i];
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
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
												[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
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

	% 

	% Save Feedback Matrices
	% ++++++++++++++++++++++
	worst_case_L.Q = value(Q);
	worst_case_L.r2 = value(r);
	worst_case_L.F2 = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	worst_case_L.u0_2 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

	worst_case_L.opt_obj = value(alpha_2);

	%Announce new experiment
	disp('=================================================================================')
	disp('Experiment 2: Synthesis For "At Most 1" Data Missing with switching observer_gain')

	clear delta
	clear Q
	clear alpha_2
	clear G

	clear noise_constrs
	clear positive_constr
	clear dual_equal_constrs
	clear l_diag_constr

	% Optimization Variables
	% ++++++++++++++++++++++
	delta 		= sdpvar(dd*T,1,'full');
	mu 			= sdpvar(p*T,1,'full');
	acc_e.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');

	alpha_l 	= sdpvar(T+1,1,'full');

	% Feedback Variables
	Q1 = sdpvar(m*T,p*T,'full');
	Q2 = sdpvar(m*T,p*T,'full');
	r1 = sdpvar(m*T,1,'full');
	r2 = sdpvar(m*T,1,'full');

	% Dual Variables
	Pi_1_1 = sdpvar(2*n*T,2*(dd+p)*T+2*n,'full');
	Pi_2_1 = sdpvar(2*n,2*(dd+p)*T+2*n,'full');

	Pi_1_2 = sdpvar(2*n*T,2*(dd+p)*T+2*n,'full');
	Pi_2_2 = sdpvar(2*n,2*(dd+p)*T+2*n,'full');

	% Creating Constraints
	% ++++++++++++++++++++

	[S0,H0,Cm0_1,xi0m,E_big] = create_skaf_n_boyd_matrices(acc_e,T,'missing',find(L(1,:) == 0)-1);
	[~,~,Cm0_2,~,~] = create_skaf_n_boyd_matrices(acc_e,T,'missing',find(L(2,:) == 0)-1);

	positive_constr = [ Pi_1_1 >= 0, Pi_2_1 >= 0, Pi_1_2 >= 0 , Pi_2_1 >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	noise_constrs = [ Pi_1_1 * [ acc_e.d * ones(2*dd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r1 ];
	noise_constrs = noise_constrs + [ Pi_2_1 * [ acc_e.d * ones(2*dd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r1 ];

	noise_constrs = noise_constrs + [ Pi_1_2 * [ acc_e.d * ones(2*dd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r2 ];
	noise_constrs = noise_constrs + [ Pi_2_2 * [ acc_e.d * ones(2*dd*T,1) ; acc_e.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r2 ];


	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; acc_e.A^i];
	end

	G1 = [ (eye(n*(T+1))+S0*Q1*Cm0_1)*S0*E_big S0*Q1 (eye(n*(T+1))+S0*Q1*Cm0_1)*pre_xi ];
	G2 = [ (eye(n*(T+1))+S0*Q2*Cm0_2)*S0*E_big S0*Q2 (eye(n*(T+1))+S0*Q2*Cm0_2)*pre_xi ];

	bounded_disturb_matrix = [ [ eye(dd*T) ; -eye(dd*T) ] zeros(2*dd*T,p*T+n) ;
								zeros(2*p*T,dd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
								zeros(2*n,(p+dd)*T) [ eye(n) ; -eye(n) ] ];

	dual_equal_constrs = [ Pi_1_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G1 ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2_1 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G1];
	dual_equal_constrs = dual_equal_constrs + [Pi_1_2 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G2 ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G2];

	%Lower Diagonal Constraint
	l_diag_constr = [];
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q1(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
												[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
	end
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q2(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
												[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
	end

	%Add constraint for the shared entries of Q
	overlap_Q_constr = [Q1([1:m*2],[1:p*2]) == Q2([1:m*2],[1:p*2])];

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim2 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+overlap_Q_constr, ...
			alpha_2, ...
			ops)

	if optim2.problem ~= 0
		error(['The design problem was not completely solved.' optim2.info ])
	end

	% Save Feedback Matrices
	% ++++++++++++++++++++++
	hybrid_obsv.Q1 = value(Q1);
	hybrid_obsv.Q2 = value(Q2);
	hybrid_obsv.r1 = value(r1);
	hybrid_obsv.r2 = value(r2);
	hybrid_obsv.F1 = value( (inv(value(eye(size(Q1,1)) + Q1*Cm0_1*S0)) ) * Q1);
	hybrid_obsv.F2 = value( (inv(value(eye(size(Q2,1)) + Q2*Cm0_2*S0)) ) * Q2);
	hybrid_obsv.u0_1 = value( inv(value(eye(size(Q1,1)) + Q1*Cm0_1*S0)) * r1 );
	hybrid_obsv.u0_2 = value( inv(value(eye(size(Q2,1)) + Q2*Cm0_2*S0)) * r2 );

	hybrid_obsv.opt_obj = value(alpha_2);

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc;
	results.experim_params.T = T;
	results.experim_params.M1 = M1;
	results.L = L;
	results.L_star = L_star;

	results.hybrid_obsv = hybrid_obsv;
	results.worst_case_L = worst_case_L;

end