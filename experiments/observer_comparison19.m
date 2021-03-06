function [ results ] = observer_comparison19( varargin )
%	observer_comparison19.m
%		Description:
%			The objective of this experiment is to design a finite horizon, affine
%			estimator (fhae) that is robust against the possibility of 1 observation
%			missing in the entire sequence of length T.
%
%			What makes this experiment different from 13,14,16 and 18 is that:
%				- It uses the recent optimization presented by Prof. Yong.
%				- We make explicit use of the "E" matrix in the ACC example.
%				- A reduced order observer is incorporated to try to achieve Equalized Performance if necessary.
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
	%Using ACC System
	load('data/system_examples/acc_p.mat');

	%Create Reduced Order Observer's Error Dynamics
	acc_roo.A = acc.A(3,3);
	acc_roo.B = eye(1);
	acc_roo.C = acc.A([1:2],3);
	acc_roo.E = [acc.E(3,1) acc.A(3,[1:2]) zeros(1,3) ];
	% acc_roo.G1 = [acc.E([1:2],1) acc.A([1:2],[1:2])];
	% acc_roo.G2 = [acc.E([1:2],1) zeros(2) ];
	acc_roo.G = [acc.E([1:2],1) acc.A([1:2],[1:2]) zeros(2,1) eye(2) ];
	acc_roo.m = acc.m;
	acc_roo.d = acc.d; 

	%Create Error System
	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));

	test_sys = acc_e;

	%Create dimensions
	n = size(test_sys.A,1);
	p = size(test_sys.C,1);

	m = size(test_sys.B,2);
	dd = size(test_sys.E,2); %Delta Dimension

	b_dim = size(acc_roo.G,2);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%Calculate Big G and E Matrices
	% G_at_each_n = {}; E_bar = {};
	% G_big = zeros(T*size(acc_roo.G,1),(T+1)*b_dim);
	% E_big = zeros(T*size(acc_roo.E,1),(T+1)*b_dim);
	% for i = 1:T
	% 	G_big( [(i-1)*size(acc_roo.G,1) + 1 : i*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(i-1)*b_dim) acc_roo.G zeros(size(acc_roo.G,1),(T-i)*b_dim )];
	% 	E_at_each_n{i} = acc_roo.E; 
	% 	E_big( [(i-1)*size(acc_roo.E,1) + 1 : i*size(acc_roo.E,1) ] , : ) = [ zeros(size(acc_roo.E,1),(i-1)*b_dim) acc_roo.E zeros(size(acc_roo.E,1),(T-i)*b_dim )];
	% end
	% E_bar = [ blkdiag(E_at_each_n{:}) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Announce New Experiments
	disp('===========================================')
	disp('Experiment 1: Synthesis For No Data Missing')
	disp(' ')
	disp('Assuming that M1 and T are given (they are constants), find the minimal M2.')
	disp('Comparing the performance of the raw calculation with the function''s performance.')
	disp(' ')

	% Constants
	% +++++++++

	% E_big = [];
	% for i = 1:T
	% 	E_big = blkdiag(E_big,test_sys.E);
	% end

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

	[S0,H0,Cm0,xi0m,E_big] = create_skaf_n_boyd_matrices(test_sys,T);

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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Try to do the same design with a function %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[contr1_fcn,optim1_fcn] = create_fhae_w_eq_recovery('min_M2',M1,T,acc,'verbosity',2)

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

	results.none_missing.fcn = optim1_fcn;

end