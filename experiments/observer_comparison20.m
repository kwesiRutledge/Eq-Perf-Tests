function [ results ] = observer_comparison20( varargin )
%	observer_comparison16.m
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

	T_s = 0.5;

	%Default Values for observer_comparison
	if nargin < 2
		T = 6;
		M1 = 0.7;
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

	test_sys = acc_roo;
	test_sys.B = eye(size(test_sys.A,1));

	%Create dimensions
	n = size(test_sys.A,1);
	p = size(test_sys.C,1);

	m = size(test_sys.B,2);
	wd = size(acc.E,2); %Delta Dimension
	vd = size(acc.C,1); %measurement disturbance shifting

	b_dim = size(acc_roo.G,2);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Announce New Experiments
	disp('========================================================')
	disp('Experiment 1: Synthesis For No Data Missing, ROO for ACC')
	disp(' ')
	disp('Assuming that M1 and T are given (they are constants), find the minimal M2.')
	disp('Comparing the performance of the raw calculation with the function''s performance.')
	disp(' ')

	% Constants
	% +++++++++

	%Calculate Big C Matrix
	G_at_each_n = {}; E_bar = {};
	G_big = zeros(T*size(acc_roo.G,1),(T+1)*(wd+vd));
	E_big = zeros(T*size(acc_roo.E,1),(T+1)*(wd+vd));

	for i = 1:T
		% G1_at_each_n{i} = acc_roo.G1;
		% G2_at_each_n{i} = acc_roo.G2;
		G_big( [(i-1)*size(acc_roo.G,1) + 1 : i*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(i-1)*(wd+vd)) acc_roo.G zeros(size(acc_roo.G,1),(T-i)*(wd+vd) )];
		E_at_each_n{i} = acc_roo.E; 
		E_big( [(i-1)*size(acc_roo.E,1) + 1 : i*size(acc_roo.E,1) ] , : ) = [ zeros(size(acc_roo.E,1),(i-1)*(wd+vd)) acc_roo.E zeros(size(acc_roo.E,1),(T-i)*(wd+vd) )];
	end
	E_bar = [ blkdiag(E_at_each_n{:}) ];

	% Optimization Variables
	% ++++++++++++++++++++++
	% delta 		= sdpvar(wd*T,1,'full');
	% mu 			= sdpvar(p*T,1,'full');
	% b replaces process and measurement noises
	b = sdpvar((wd+vd)*(T+1),1,'full');
	%Last process noise unused

	test_sys.x0 = sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');

	alpha_l 	= sdpvar(T+1,1,'full');

	% Feedback Variables
	Q = sdpvar(m*T,p*T,'full');
	r = sdpvar(m*T,1,'full');

	% Dual Variables
	Pi_1 = sdpvar(2*n*T,2*(wd+p)*(T+1)+2*n,'full');
	Pi_2 = sdpvar(2*n,2*(wd+p)*(T+1)+2*n,'full');

	% Creating Constraints
	% ++++++++++++++++++++

	[S0,H0,Cm0,xi0m] = create_skaf_n_boyd_matrices(test_sys,T);

	positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	% Constraints that bound the noise
	sel_noise_mat = [];
	for i = 1 : T+1
		sel_noise_mat = [ sel_noise_mat; test_sys.d * ones(wd,1) ; test_sys.m * ones(vd,1) ];
	end
	sel_noise_mat = [ sel_noise_mat ; sel_noise_mat ]; %Repeat the same bounds (for the negative constraint)
	sel_noise_mat = [sel_noise_mat; M1*ones(2*n,1)];

	noise_constrs = [ Pi_1 * sel_noise_mat <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
	noise_constrs = noise_constrs + [ Pi_2 * sel_noise_mat <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; test_sys.A^i];
	end

	G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*E_big+S0*Q*G_big (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

	% bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,p*T+n) ;
	% 							zeros(2*p*T,wd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
	% 							zeros(2*n,(p+wd)*T) [ eye(n) ; -eye(n) ] ];

	bounded_disturb_matrix = [];
	bounded_disturb_matrix = [ bounded_disturb_matrix ; [eye((wd+vd)*(T+1)) ; -eye((wd+vd)*(T+1))] zeros(2*(wd+vd)*(T+1),n) ];
	bounded_disturb_matrix = [ bounded_disturb_matrix ; zeros(2*n,(wd+vd)*(T+1)) [ eye(n) ; -eye(n) ] ];

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

	opt_obj1 = value(alpha_2);
	Q1 = value(Q);
	r1 = value(r);
	F1 = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	u0_1 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Find the Minimum Level M_star such that the ROO can achieve Equalized Performance %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	perf_level = 1;
	T = 1; %this is for equalized performance

	delta_pl = perf_level/2;
	search_finished = false

	while ~search_finished
		% Constants
		G_big3 = zeros(T*size(acc_roo.G,1),(T+1)*(wd+vd));
		E_big3 = zeros(T*size(acc_roo.E,1),(T+1)*(wd+vd));

		for i = 1:T
			G_big3( [(i-1)*size(acc_roo.G,1) + 1 : i*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(i-1)*(wd+vd)) acc_roo.G zeros(size(acc_roo.G,1),(T-i)*(wd+vd) )];
			E_big3( [(i-1)*size(acc_roo.E,1) + 1 : i*size(acc_roo.E,1) ] , : ) = [ zeros(size(acc_roo.E,1),(i-1)*(wd+vd)) acc_roo.E zeros(size(acc_roo.E,1),(T-i)*(wd+vd) )];
		end

		b = sdpvar((wd+vd)*(T+1),1,'full');
		%Last process noise unused

		test_sys.x0 = sdpvar(n,1,'full');

		alpha_2 	= sdpvar(1,1,'full');

		% Feedback Variables
		Q = sdpvar(m*T,p*T,'full');
		r = sdpvar(m*T,1,'full');

		% Dual Variables
		Pi_1 = sdpvar(2*n*T,2*(wd+p)*(T+1)+2*n,'full');
		Pi_2 = sdpvar(2*n,2*(wd+p)*(T+1)+2*n,'full');

		% Creating Constraints
		% ++++++++++++++++++++

		[S3,H3,Cm3,xi0m_3] = create_skaf_n_boyd_matrices(test_sys,T);

		positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
		end

		% Constraints that bound the noise
		sel_noise_mat = [];
		for i = 1 : T+1
			sel_noise_mat = [ sel_noise_mat; test_sys.d * ones(wd,1) ; test_sys.m * ones(vd,1) ];
		end
		sel_noise_mat = [ sel_noise_mat ; sel_noise_mat ]; %Repeat the same bounds (for the negative constraint)
		sel_noise_mat = [sel_noise_mat; perf_level*ones(2*n,1)];

		noise_constrs = [ Pi_1 * sel_noise_mat <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S3*r ];
		noise_constrs = noise_constrs + [ Pi_2 * sel_noise_mat <= perf_level * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S3*r ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T
			pre_xi = [ pre_xi ; test_sys.A^i];
		end

		G = [ (eye(n*(T+1))+S3*Q*Cm3)*S3*E_big3+S3*Q*G_big3 (eye(n*(T+1))+S3*Q*Cm3)*pre_xi ];

		% bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,p*T+n) ;
		% 							zeros(2*p*T,wd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
		% 							zeros(2*n,(p+wd)*T) [ eye(n) ; -eye(n) ] ];

		bounded_disturb_matrix = [];
		bounded_disturb_matrix = [ bounded_disturb_matrix ; [eye((wd+vd)*(T+1)) ; -eye((wd+vd)*(T+1))] zeros(2*(wd+vd)*(T+1),n) ];
		bounded_disturb_matrix = [ bounded_disturb_matrix ; zeros(2*n,(wd+vd)*(T+1)) [ eye(n) ; -eye(n) ] ];

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
		optim3 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr, ...
				alpha_2, ...
				ops);

		% if optim3.problem ~= 0
		% 	error(['The design problem was not completely solved.' optim1.info ])
		% end

		%Check feasibility + delta_pl
		if (optim3.problem == 0) & (delta_pl<=0.01)
			search_finished = true;
		end

		if optim3.problem == 0
			%If feasible, reduce the performance level
			perf_level = perf_level - delta_pl;
		else
			perf_level = perf_level + delta_pl;
			delta_pl = delta_pl/2;
			disp(['delta_pl = ' num2str(delta_pl)])
		end
	end

	opt_obj3 = value(alpha_2);
	Q3 = value(Q);
	r3 = value(r);
	F3 = value( (inv(value(eye(size(Q,1)) + Q*Cm3*S3)) ) * Q);
	u0_3 = value( inv(value(eye(size(Q,1)) + Q*Cm3*S3)) * r );

	disp('Finished search for optimal performance level.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller for 1 missing Data Case %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	clear delta
	clear mu
	clear alpha_l

	%Announce new experiment
	disp('====================================================')
	disp('Experiment 2: Synthesis For "At Most 1" Data Missing')

	T = 6;
	T_allowed = 3;
	M1 = 1.5*perf_level;
	num_potential_missing = 2;

	T_not_feasible = true;

	while T_not_feasible

		% Constants
		G_big2 = zeros(T*size(acc_roo.G,1),(T+1)*(wd+vd));
		E_big2 = zeros(T*size(acc_roo.E,1),(T+1)*(wd+vd));

		for ind0 = 1:T
			G_big2( [(ind0-1)*size(acc_roo.G,1) + 1 : ind0*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(ind0-1)*(wd+vd)) acc_roo.G zeros(size(acc_roo.G,1),(T-ind0)*(wd+vd) )];
			E_big2( [(ind0-1)*size(acc_roo.E,1) + 1 : ind0*size(acc_roo.E,1) ] , : ) = [ zeros(size(acc_roo.E,1),(ind0-1)*(wd+vd)) acc_roo.E zeros(size(acc_roo.E,1),(T-ind0)*(wd+vd) )];
		end

		%Optimization Variables
		test_sys.x0 = sdpvar(n,1,'full');
		alpha_2 	= sdpvar(1,1,'full');
		alpha_l 	= sdpvar(T+1,1,'full');

		% Feedback Variables
		Q = sdpvar(m*T,p*T,'full');
		r = sdpvar(m*T,1,'full');

		% Dual Variables
		Pi_1 = sdpvar(2*n*T,2*(wd+p)*(T+1)+2*n,'full');
		Pi_2 = sdpvar(2*n,2*(wd+p)*(T+1)+2*n,'full');

		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m] = create_skaf_n_boyd_matrices(test_sys,T);

		positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for ind0 = 1 : T
			sel_influenced_states = [ sel_influenced_states ; select_m(ind0,T) ];
		end

		% Constraints that bound the noise
		sel_noise_mat = [];
		for t = 0 : T
			sel_noise_mat = [ sel_noise_mat; test_sys.d * ones(wd,1) ; test_sys.m * ones(vd,1) ];
		end
		sel_noise_mat = [ sel_noise_mat ; sel_noise_mat ]; %Repeat the same bounds (for the negative constraint)
		sel_noise_mat = [sel_noise_mat; M1*ones(2*n,1)];

		noise_constrs = [ Pi_1 * sel_noise_mat <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
		noise_constrs = noise_constrs + [ Pi_2 * sel_noise_mat <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

		%Dual relationship to design variables
		pre_xi = [];
		for ind0 = 0:T
			pre_xi = [ pre_xi ; test_sys.A^ind0];
		end

		G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*E_big2+S0*Q*G_big2 (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

		% bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,p*T+n) ;
		% 							zeros(2*p*T,wd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
		% 							zeros(2*n,(p+wd)*T) [ eye(n) ; -eye(n) ] ];

		bounded_disturb_matrix = [];
		bounded_disturb_matrix = [ bounded_disturb_matrix ; [eye((wd+vd)*(T+1)) ; -eye((wd+vd)*(T+1))] zeros(2*(wd+vd)*(T+1),n) ];
		bounded_disturb_matrix = [ bounded_disturb_matrix ; zeros(2*n,(wd+vd)*(T+1)) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

		%Lower Diagonal Constraint
		l_diag_constr = [];
		for bl_row_num = 1 : T-1
			l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(test_sys.B,2)+1:bl_row_num*size(test_sys.B,2)], ...
													[bl_row_num*size(test_sys.C,1)+1:end] ) == 0 ];
		end

		% Can we simply try to add further constraints?
		% +++++++++++++++++++++++++++++++++++++++++++++

		for missing_num = 1:num_potential_missing
			
			all_poss_comb = nchoosek([1:T_allowed],missing_num);

			for comb_num = 1:size(all_poss_comb,1)

				missing_locs = unique([ all_poss_comb(comb_num,:) (all_poss_comb(comb_num,:)+1) ]);

				%Calculate Big C Matrix
				[~,~,Cm,~] = create_skaf_n_boyd_matrices(test_sys,T,'missing',missing_locs);

				%Create Special selection matrices for selecting the proper variables
				R = [ zeros(n,n*T) eye(n) ];
				% mu_select = [];
				% for i = 1:T+1
				% 	%Takes out the row representing the observation at time t=i+1
				% 	if any(missing_loc+1+1 == i)
				% 		mu_select = [ mu_select ; zeros(wd+vd,(wd+vd)*(T+1)) ];
				% 	else
				% 		mu_select = [ mu_select ; [ zeros((wd+vd),vd*(i-1)) eye(wd+vd) zeros(wd+vd,p*(T-i)) ] ];
				% 	end
				% end

				% Pxd = [ Pxd ; (eye(n*(T+1))+S*Q*Cm)*S ];%E_bar*[ eye(b_dim*T) zeros(b_dim*T,b_dim) ];
				% Pxm = [ Pxm ; S*Q*mu_select ];
				% xi_tilde = [ xi_tilde ; (eye(n*(T+1)) + S*Q*Cm)*xi0m + S*r];

				G_big_temp = zeros(T*size(acc_roo.G,1),(T+1)*(wd+vd));
				for i = 1:T-1
					if any( i == (missing_locs+1)) || (any(i == (missing_locs+2) ))
						G_big_temp( [(i-1)*size(acc_roo.G,1) + 1 : i*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(T+1)*(wd+vd)) ];
					else
						G_big_temp( [(i-1)*size(acc_roo.G,1) + 1 : i*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(i-1)*(wd+vd)) acc_roo.G zeros(size(acc_roo.G,1),(T-i)*(wd+vd) )];
					end
				end

				G = [ (eye(n*(T+1))+S0*Q*Cm)*S0*E_big2+S0*Q*G_big_temp (eye(n*(T+1))+S0*Q*Cm)*pre_xi ];

				%Awd to the constraint set
				dual_equal_constrs = dual_equal_constrs + [Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G];
				dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

			end

		end

		% OPTIMIZE
		% ++++++++

		% ops = sdpsettings('verbose',verbosity);
		optim2 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr, ...
				alpha_2, ...
				ops)

		optim2.info

		if optim2.problem == 0
			T_not_feasible = false
		else
			T = T+1
		end

	end

	if optim2.problem ~= 0
		error(['The design problem was not completely solved.' optim2.info ])
	end

	% Save Feedback Matrices
	% ++++++++++++++++++++++
	Q2 = value(Q);
	r2 = value(r);
	F2 = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	u0_2 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );
	T_opt = T;

	opt_obj2 = value(alpha_2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot the Performance of these 2 Controllers %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	num_rollouts = 1000;

	T1 = 10; %First part of the trajectory (no missing data)
	T2 = T_opt;

	% Create the first part of the trajectory.
	% ++++++++++++++++++++++++++++++++++++++++
	x0_t = unifrnd(-M1,M1,n,num_rollouts);
	w1_t = unifrnd(-test_sys.d,test_sys.d,wd*(T1+1),num_rollouts);
	v1_t = unifrnd(-test_sys.m,test_sys.m,vd*(T1+1),num_rollouts);

	%Create b1
	b1_t = [];
	for t = 0:T1
		b1_t = [b1_t; w1_t([t*wd+1:(t+1)*wd],:) ; v1_t([t*vd+1:(t+1)*vd],:)];
	end

	x1_t = x0_t;
	for t = 1:T1
		x1_t = [ x1_t ; (test_sys.A+F3*test_sys.C)*x1_t([(t-1)*n+1:t*n],:) + (F3*test_sys.G+test_sys.E)*b1_t([(t-1)*(wd+vd)+1:(t+1)*(wd+vd)],:) ];
	end

	%Collect Magnitudes of x1_t
	x1_mag_t = [];
	for rollout_ind = 1:num_rollouts
		for t = 0:T1
			x1_mag_t(t+1,rollout_ind) = norm(x1_t([t*n+1:(t+1)*n],rollout_ind),Inf);
		end
	end

	% Create the second part of the trajectory.
	% +++++++++++++++++++++++++++++++++++++++++

	w2_t = unifrnd(-test_sys.d,test_sys.d,wd*(T2+1),num_rollouts);
	v2_t = unifrnd(-test_sys.m,test_sys.m,vd*(T2+1),num_rollouts);

	%Create b2
	b2_t = [];
	for t = 0:T2
		b2_t = [b2_t; w2_t([t*wd+1:(t+1)*wd],:) ; v2_t([t*vd+1:(t+1)*vd],:)];
	end

	poss_miss_indices = [1:T_allowed];

	x2_t0 = x1_t([T1*n+1:(T1+n)],:);
	for rollout_ind = 1 : num_rollouts
		rnd_missing_data_pattern = logical([]);
		for t = 1:T_allowed
			if t == 1
				rnd_missing_data_pattern(t) = ((unidrnd(2)-1) == 1);
			elseif sum(rnd_missing_data_pattern) >= num_potential_missing
				rnd_missing_data_pattern(t) = false;
			else
				rnd_missing_data_pattern(t) = ((unidrnd(2)-1) == 1);
			end
		end

		missing_locs2 = [];
		if any(rnd_missing_data_pattern)
			missing_locs2 = poss_miss_indices(rnd_missing_data_pattern);
		end

		% Constants
		G_big_t = zeros(T2*size(acc_roo.G,1),(T2+1)*(wd+vd));
		E_big_t = zeros(T2*size(acc_roo.E,1),(T2+1)*(wd+vd));

		for t = 1:T2
			% G_big_t( [(i-1)*size(acc_roo.G,1) + 1 : i*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(i-1)*(wd+vd)) acc_roo.G zeros(size(acc_roo.G,1),(T2-i)*(wd+vd) )];
			E_big_t( [(t-1)*size(acc_roo.E,1) + 1 : t*size(acc_roo.E,1) ] , : ) = [ zeros(size(acc_roo.E,1),(t-1)*(wd+vd)) acc_roo.E zeros(size(acc_roo.E,1),(T2-t)*(wd+vd) )];

			if any( t == (missing_locs2+1)) || (any(t == (missing_locs2+2) ))
				G_big_t( [(t-1)*size(acc_roo.G,1) + 1 : t*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(T+1)*(wd+vd)) ];
			else
				G_big_t( [(t-1)*size(acc_roo.G,1) + 1 : t*size(acc_roo.G,1) ] , : ) = [ zeros(size(acc_roo.G,1),(t-1)*(wd+vd)) acc_roo.G zeros(size(acc_roo.G,1),(T-t)*(wd+vd) )];
			end
		end

		%Create initial conditions for the next part
		test_sys.x0 = x2_t0(:,rollout_ind);
		if ~any(rnd_missing_data_pattern)
			[S0,H0,Cm0,xi0m] = create_skaf_n_boyd_matrices(test_sys,T2);
		else
			[S0,H0,Cm0,xi0m] = create_skaf_n_boyd_matrices(test_sys,T2,'missing',unique([missing_locs2 (missing_locs2+1) ]) );
		end

		%Calculate the next entries for x2_t
		x2_t(:,rollout_ind) = ((S0 + S0*F2*inv(eye(p*T2) - Cm0 * S0 * F2) * Cm0 * S0) * E_big_t + S0*F2*inv(eye(p*T2) - Cm0 * S0 * F2) * G_big_t )*b2_t(:,rollout_ind) + (eye(n*(T2+1)) + S0*F2*inv(eye(p*T2) - Cm0 * S0 * F2) * Cm0 ) * ( xi0m + S0 * u0_2 );

		for t = 0 : T2
			x2_mag_t(t+1,rollout_ind) = norm(select_m(t,T2)*x2_t(:,rollout_ind),Inf);
		end

	end


	%Combine x1 and x2
	x_mag_t = [ x1_mag_t ; x2_mag_t([2:end],:) ];

	figure;
	hold on;

	bar([0:T1+T2]*T_s,[ M1 M1*ones(1,T1-1) M1 opt_obj2*ones(1,T2-1) M1 ],'w')
	for rollout_ind = 1 : num_rollouts
		plot([0:T1+T2]*T_s,x_mag_t(:,rollout_ind))
	end

	xlabel('Time [sec]')
	ylabel('$||v_L(t)-\hat{v}_L(t)||_{\infty}$','Interpreter','latex')
	legend('Guarantees')
	axis([-0.5*T_s (T1+T2+0.5)*T_s 0 0.6])
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

	results.missing_n.optimization = optim2;
	results.missing_n.opt_obj = opt_obj2;
	results.missing_n.Q = Q2;
	results.missing_n.r = r2;
	results.missing_n.F = F2;
	results.missing_n.u0 = u0_2;
	results.missing_n.T_opt = T_opt;

	results.opt_pl.perf_level = perf_level;

end