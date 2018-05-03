function [ results ] = observer_comparison27( varargin )
%	observer_comparison27.m
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
			1,1,0,1,1,1;
			1,1,1,0,1,1;
			1,1,1,1,0,1 ];

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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Onto the new version of the optimization %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear Q;
	clear r;
	clear Pi_1
	clear Pi_2;

	sys0 = acc_e;

	% Optimization Variables
	% ++++++++++++++++++++++
	delta 		= sdpvar(dd*T,1,'full');
	mu 			= sdpvar(p*T,1,'full');
	sys0.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');

	alpha_l 	= sdpvar(T+1,1,'full');

	for pattern_ind = 1 : size(L,1)

		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T,p*T,'full');
		r{pattern_ind} = sdpvar(m*T,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(2*n*T,2*(dd+p)*T+2*n,'full');
		Pi_2{pattern_ind} = sdpvar(2*n,2*(dd+p)*T+2*n,'full');
	end

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];

	for pattern_ind = 1 : size(L,1)
		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,E_big] = create_skaf_n_boyd_matrices(sys0,T,'missing',find(L(pattern_ind,:) == 0)-1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
		end

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ sys0.d * ones(2*dd*T,1) ; sys0.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*kron(eye(T),sys0.B)*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ sys0.d * ones(2*dd*T,1) ; sys0.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*kron(eye(T),sys0.B)*r{pattern_ind} ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T
			pre_xi = [ pre_xi ; sys0.A^i];
		end

		G = [ (eye(n*(T+1))+S0*kron(eye(T),sys0.B)*Q{pattern_ind}*Cm0)*S0*E_big S0*kron(eye(T),sys0.B)*Q{pattern_ind} (eye(n*(T+1))+S0*kron(eye(T),sys0.B)*Q{pattern_ind}*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(dd*T) ; -eye(dd*T) ] zeros(2*dd*T,p*T+n) ;
									zeros(2*p*T,dd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
									zeros(2*n,(p+dd)*T) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

		%Lower Diagonal Constraint
		for bl_row_num = 1 : T-1
			l_diag_constr = l_diag_constr + [ Q{pattern_ind}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
																[bl_row_num*p+1:end] ) == 0 ];
		end

		%Add joint constraints for all 
		disp(size(L,1))
		disp(pattern_ind)
		for patt_i = pattern_ind+1:size(L,1)
			%Match
			p1 = L(pattern_ind,:);
			p2 = L(patt_i,:);
			% Add xor
			p_overlap = bitxor(p1,p2);
			%Bit at which things end up being different
			ind_identical = find(p_overlap,1) - 1;
			%Add constraints
			shared_Q_constrs = shared_Q_constrs +  [Q{pattern_ind}( [1:ind_identical*m] , [1:ind_identical*p] ) == Q{patt_i}( [1:ind_identical*m] , [1:ind_identical*p] )];
			shared_r_constrs = shared_r_constrs +  [r{pattern_ind}( [1:ind_identical*m] ) == r{patt_i}([1:ind_identical*m]) ];
		end

	end

	% OPTIMIZATION
	% ++++++++++++
	ops = sdpsettings('verbose',verbosity);
	optim2 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr + shared_Q_constrs + shared_r_constrs, ...
			alpha_2, ...
			ops)

	if optim2.problem ~= 0
		error(['The design problem was not completely solved.' optim2.info ])
	end

	% 

	% Save Feedback Matrices
	% ++++++++++++++++++++++s
	for pattern_ind = 1 : size(L,1)
		%Get Parameters
		[S0,H0,Cm0,~,E_big] = create_skaf_n_boyd_matrices(sys0,T,'missing',find(L(pattern_ind,:) == 0)-1);

		state_based.Q{pattern_ind} = value(Q{pattern_ind});
		state_based.r{pattern_ind} = value(r{pattern_ind});
		state_based.F{pattern_ind} = value( (inv(value(eye(size(S0,2)) + kron(eye(T),sys0.B)*Q{pattern_ind}*Cm0*S0)) ) *kron(eye(T),sys0.B)* Q{pattern_ind});
		state_based.u0{pattern_ind} = value( inv(value(eye(size(S0,2)) + kron(eye(T),sys0.B)*Q{pattern_ind}*Cm0*S0)) *kron(eye(T),sys0.B)* r{pattern_ind} );
	end

	state_based.opt_obj = value(alpha_2);


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %% Plot the Performance of this Designed Controller %%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% %Controller Definition
	% contr2.F = F2;
	% contr2.u0 = u0_2;

	% num_plots = 4;

	% %Plot Constants
	% [xi_t0,xi_mag_t0] = apply_controller_to_rollouts(acc_e,contr2,T,num_rollouts,M1);
	% xi{1} = xi_t0;
	% xi_mag{1} = xi_mag_t0;

	% figure('units','normalized','outerposition',[0 0 0.5 0.5]);
	% for i = 1 : num_plots
	% 	if i > 1
	% 		[xi{i},xi_mag{i}] = apply_controller_to_rollouts(acc_e,contr2,T,num_rollouts,M1,'missing',i-1);
	% 	end

	% 	subplot(2,2,i)
	% 	hold on;

	% 	bar([0:T]*T_s,[ M1 value(alpha_2)*ones(1,T-1) M1 ],'w')
	% 	for r_num = 1:rollouts_per_plot
	% 		plot([0:T]*T_s,xi_mag{i}(:,r_num))
	% 	end
	% 	xlabel('Time [sec]')
	% 	ylabel('$||x(t)-\hat{x}(t)||_{\infty}$','Interpreter','latex')
	% 	% if i > 1
	% 	% 	title(['Norm of the Estimation Error (Missing Data Occurs at t=' num2str(i-1) ')'])
	% 	% else
	% 	% 	title(['Norm of the Estimation Error (Data Always Available)'])
	% 	% end
	% 	legend('Guarantees')
	% 	xt = get(gca, 'XTick');
	% 	set(gca, 'FontSize', 16)
	% end


	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc;
	results.experim_params.T = T;
	results.experim_params.M1 = M1;
	results.L = L;
	results.L_star = L_star;

	results.state_based = state_based;
	results.worst_case_L = worst_case_L;

end

