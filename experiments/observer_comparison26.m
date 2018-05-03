function [ results ] = observer_comparison26( varargin )
%observer_comparison26.m
%Description:
%	The objective of this experiment is to test the choice of v == 0 when data is missing instead of
%	the current, incorrect method of not doing anything. 
%
%	Qualities of this Problem:
%			- Robust Optimization proposed by Sze Zheng Yong
%			- Worst Case Language L^* used.
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
	T = 6;

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

	paper_format = true;

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
	disp('Experiment 1: Synthesis For "At Most 1" Data Missing')
	disp('- Only can be missing in first 3 instants')
	disp('- Legacy method, mathematically incorrect')
	disp(' ')

	% Constants
	% +++++++++

	M1 = 1;
	perf_level = M1;
	L = [	1,0,1,1,1,1;
			1,1,0,1,1,1;
			1,1,1,0,1,1;
			1,1,1,1,0,1;
			1,1,1,1,1,1 ];

	L_star = ones(1,6);
	for sig_i = 1:size(L,1)
		L_star = bitand(L_star,L(sig_i,:))
	end



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

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Solving for Worst Case Word while overconstraining Q %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Creating Constraints
	% ++++++++++++++++++++

	[S0,H0,Cm0,xi0m,E_big] = create_skaf_n_boyd_matrices(acc_e,T,'missing',find(L_star==0)-1);

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
	Q_val = value(Q);
	r_val = value(r);
	F = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	u0 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

	opt_obj2 = value(alpha_2);

	v_free.optim = optim2;
	v_free.Q = Q_val;
	v_free.r = r_val;
	v_free.F = F;
	v_free.u0 = u0;
	v_free.M2 = opt_obj2;

	clear positive_constr
	clear noise_constrs
	clear dual_equal_constrs
	clear l_diag_constr

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Solving for Worst Case Word while properly Q %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Creating Constraints
	% ++++++++++++++++++++

	[S0,H0,Cm0,xi0m,E_big] = create_skaf_n_boyd_matrices(acc_e,T,'missing',find(L_star==0)-1);

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

	% Disturbance v=0
	v_constr = [];
	for v_ind = find(L_star==0)
		v_constr = v_constr + [ v([(v_ind-1)*p+1:v_ind*p])==0 ];
	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim2 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr + v_constr, ...
			alpha_2, ...
			ops);

	if optim2.problem ~= 0
		error(['The design problem was not completely solved.' optim2.info ])
	end

	% Save Feedback Matrices
	% ++++++++++++++++++++++
	Q_val = value(Q);
	r_val = value(r);
	F = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	u0 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

	opt_obj2 = value(alpha_2);

	v_constrained.optim = optim2;
	v_constrained.Q = Q_val;
	v_constrained.r = r_val;
	v_constrained.F = F;
	v_constrained.u0 = u0;
	v_constrained.M2 = opt_obj2;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot the Performance of this Designed Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Constants
	T_s = 0.5;

	%Controller Definition
	contr1.F = v_free.F;
	contr1.u0 = v_free.u0;

	contr2.F = v_constrained.F;
	contr2.u0 = v_constrained.u0;

	controllers = [contr1,contr2];

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

	for contr_num = 1:length(controllers)
		
		contr_t = controllers(contr_num);

		%Plot Constants
		[xi_t0,xi_mag_t0] = apply_controller_to_rollouts(acc_e,contr_t,T,num_rollouts,M1);
		xi = xi_t0;
		xi_mag = xi_mag_t0;

		for missing_ob_num = (find(L_star == 0)-1) %1:3
			[S,H,Cm,~] = create_skaf_n_boyd_matrices(acc_e,T,'missing',missing_ob_num);

			% Create some trajectory matrices
			Pxd = (eye(n*(T+1))+H*contr_t.F*inv(eye(p*T)-Cm*H*contr_t.F)*Cm)*S*E_big;
			Pxm = H*contr_t.F*inv(eye(p*T)-Cm*H*contr_t.F);
			xi_factor = H*contr_t.F*inv(eye(p*T)-Cm*H*contr_t.F )*Cm;

			delta 	= unifrnd(-acc_e.d,acc_e.d,wd*T,num_rollouts);
			mu 		= unifrnd(-acc_e.m,acc_e.m,vd*T,num_rollouts);
			xi_0 	= xi(end-n+1:end,:);

			mu([missing_ob_num*p+1:(missing_ob_num+1)*p],:) = 0;

			xi_temp =  A_col*xi_0 + H*repmat(contr_t.u0,1,num_rollouts) + ...
		        xi_factor*(A_col*xi_0 + H*repmat(contr_t.u0,1,num_rollouts)) + ...
		        Pxd * delta + Pxm * mu;

		    xi = [xi;xi_temp([n+1:end],:)];
		    for rollout_ind = 1 : num_rollouts 
				for t = 1:T
					xi_mag(1+(missing_ob_num-find(L_star == 0,1)+1)*T+t,rollout_ind) = norm(select_m(t,T)*xi_temp(:,rollout_ind),Inf);
				end
			end
		end

		size(xi_mag)

		%Create bars
		t_test = [0:T*(length(find(L_star == 0)))]*T_s;
		bar_heights = [ M1 ];
		for i = 1:length(find(L_star==0))
			bar_heights = [ bar_heights opt_obj2*ones(1,T-1) M1];
		end

		figure;
		hold on;
		bar(t_test,bar_heights,'w')
		for rollout_ind = 1 : num_rollouts
			plot(t_test,xi_mag(:,rollout_ind)' )
		end

		legend('Guarantees')
		axis([-0.5*T_s (T*(3+1)+0.5)*T_s 0 opt_obj2+0.5])

		xlabel('Time [sec]')
		ylabel('$||x(t)-\hat{x}(t)||_{\infty}$','Interpreter','latex')
		% title(['Estimator Performance when $M_1$=' num2str(M1_list(controller_num))],'Interpreter','latex')
		if ~paper_format
			title(['Estimator #' num2str(contr_num)])
		end
		
		xt = get(gca, 'XTick');
		set(gca, 'FontSize', 16)
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc;
	results.experim_params.T = T;
	results.experim_params.M1s = M1_list;

	results.v_free = v_free;
	results.v_constrained = v_constrained;

end