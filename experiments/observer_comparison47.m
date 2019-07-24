function [results] = observer_comparison47(varargin)
	%observer_comparison47.m
	%
	%Description:
	%	The objective of this is to prototype a form of "leader-follower" dynamics. Version 4.
	%
	%Assumptions:
	%	- I am currently assuming that the value of the initial estimation error and the state are the same.
	%	  This only works if we assume that the initial estimate is the zero vector.

	%% Constants
	%Model Parameters
	h = 0.2;
	num_sys = 5;
	tau_i = [0.1, 0.5,0.5,2,2];

	n_i = 3;
	p_i = n_i;

	eta_w_i = [0.06,0.0015,0.0015,0.0005,0.0005];
	eta_v = 0.002;

	%Meta Model Parameters
	r = 2; %Dimension for the follower cube.

	% ad_arr = [];
	for sys_idx = 1:num_sys
		temp_A = [1,h,0;0,1,h; 0,0,1-h/tau_i(sys_idx)];
		temp_B = [0;0;h/tau_i(sys_idx)];
		temp_B_w = [0;0;1];

		ad_arr(sys_idx) = Aff_Dyn(temp_A,temp_B,zeros(n_i,1),eye(n_i), eta_w_i(sys_idx), eta_v, temp_B_w , eye(p_i) );
		% ad_arr(sys_idx) = Aff_Dyn(kron(eye(2),temp_A),kron(eye(2),temp_B),zeros(2*n_i,1),eye(2*n_i), eta_w_i(sys_idx), eta_v, kron(eye(2),temp_B_w) , eye(2*p_i) );

	end

	results.ad_arr = ad_arr;

	%Define Meta Model
	w_poly = 1;
	for sys_idx = 1:num_sys
		w_poly = w_poly * Polyhedron('lb',-eta_w_i(sys_idx),'ub',eta_w_i(sys_idx));
	end

	ad0 = 		Aff_Dyn(	blkdiag(ad_arr(:).A),blkdiag(ad_arr(:).B),zeros(n_i*num_sys,1),eye(n_i*num_sys), ...
							w_poly, Polyhedron('lb',-eta_v*ones(1,p_i*num_sys),'ub',eta_v*ones(1,p_i*num_sys)), ...
							blkdiag(ad_arr(:).B_w),blkdiag(ad_arr(:).C_v) );

	results.ad = ad0;

	%% Generate Reference Trajectories

	contr_durs = [4,4,4];
	contr_vals = [1,0.5,0.75];

	%x_r = zeros(n_i,1,num_sys);

	x0_bar(:,1,1) = [0.5;0.5;0]; x0_bar(:,1,2) = [-0.5;0;0]; x0_bar(:,1,3) = [-2;1;0]; x0_bar(:,1,4) = [-4;0;0]; x0_bar(:,1,5) = [-5;1;0];
	x_r = x0_bar;
	for sys_idx = 1:num_sys	
		for time_idx = 1:sum(contr_durs)
			for dur_idx = 1:length(contr_durs)
				if dur_idx == 1
					if (time_idx >= 1) && (time_idx <= sum(contr_durs(1:dur_idx)))
						%
						x_k = x_r(:,time_idx,sys_idx);
						x_kp1 = ad_arr(sys_idx).A*x_k + ad_arr(sys_idx).B*contr_vals(dur_idx);
						x_r(:,time_idx+1,sys_idx) = x_kp1;
					end
				else
					if (time_idx > sum(contr_durs(1:dur_idx-1))) && (time_idx <= sum(contr_durs(1:dur_idx)))
						x_k = x_r(:,time_idx,sys_idx);
						x_kp1 = ad_arr(sys_idx).A*x_k + ad_arr(sys_idx).B*contr_vals(dur_idx);
						x_r(:,time_idx+1,sys_idx) = x_kp1;
					end
				end
			end
		end
	end

	results.ref_trajs = x_r;

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize a Filter %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	clear r

	unit_box = Polyhedron('lb',-ones(1,5*n_i),'ub',ones(1,5*n_i));

	M1 = 0.3*unit_box;
	M2 = 4*unit_box;
	M3 = 0.5*M1;

	mu2 = sdpvar(1,1,'full');

	%L = {[1,0,1,1,1,1],[1,1,0,1,1,1],[1,1,1,0,1,1],[1,1,1,1,0,1],[1,0,0,1,1,1]};
	%L = {ones(1,6)};
	L = {[1,0,1,1,1,1],[1,1,0,1,1,1]};

	% Constants
	n = size(ad0.A,1);
	m = size(ad0.B,2);
	p = size(ad0.C,1);
	wd = size(ad0.B_w,2);
	vd = size(ad0.C_v,2);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	results.params.L = L;
	results.params.M1 = M1;
	results.params.M2 = M2;
	results.params.M3 = M3;

	%%Perform Optimization %%

	% Optimization Variables
	% ++++++++++++++++++++++
	ad0.x0 	= sdpvar(n,1,'full');

	ops = sdpsettings('verbose',1);

	max_T_i = -1;
	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
		r{pattern_ind} = sdpvar(m*T_i,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(2*n*T_i,2*(wd+vd)*T_i+2*n,'full');
		Pi_2{pattern_ind} = sdpvar(2*n,2*(wd+vd)*T_i+2*n,'full');

		%Find the maximum T_i
		if T_i > max_T_i
			max_T_i = T_i;
		end
	end
	w	= sdpvar(wd*max_T_i,1,'full');

	obj_constrs = [];
	obj_fcn = [];

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_constrs = [];
	positive_constr = [];
	l_diag_constr = [];

	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Creating Constraints
		% ++++++++++++++++++++

		[H0,S0,Cm0,J0,f_bar,B_w_big,~] = get_mpc_matrices(ad_arr,'word',L{pattern_ind}+1);
		P_wT = 1; P_vT = 1;
		for t_idx = 1:T_i
			P_wT = P_wT*ad0.P_w;
			P_vT = P_vT*ad0.P_v;
		end

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		sel_influenced_states = [];
		prod_M2 = 1;
		for i = 1 : T_i
			%Select all influenced states
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
			%Create product for M2
			prod_M2 = prod_M2*M2;
		end

		dual_constrs = dual_constrs + [ Pi_1{pattern_ind} * [ P_wT.b ; P_vT.b ; M1.b ] <= mu2*ones(size(prod_M2.A,1),1) - prod_M2.A*sel_influenced_states*(S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar) ];
		dual_constrs = dual_constrs + [ Pi_2{pattern_ind} * [ P_wT.b ; P_vT.b ; M1.b ] <= M3.b - M3.A*select_m(T_i,T_i)*(S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar) ];

		G{pattern_ind} = [ 	(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*B_w_big ...
									S0*Q{pattern_ind}*C_v_big ...
									(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*J0 ];

		bounded_disturb_matrix = [ 	P_wT.A zeros(size(P_wT.A,1),vd*T_i+n) ;
									zeros(size(P_vT.A,1),size(P_wT.A,2)) P_vT.A zeros(size(P_vT.A,1),n) ;
									zeros(size(M1.A,1),(vd+wd)*T_i) M1.A ];

		dual_constrs = dual_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == kron(eye(T_i),unit_box.A)*sel_influenced_states*G{pattern_ind} ];
		dual_constrs = dual_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == M3.A*select_m(T_i,T_i)*G{pattern_ind}];

		%Lower Diagonal Constraint
		for bl_row_num = 1 : T_i-1
			l_diag_constr = l_diag_constr + [ Q{pattern_ind}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
																[bl_row_num*p+1:end] ) == 0 ];
		end

		%Awd joint constraints for all 
		for patt_i = pattern_ind+1:length(L)
			%Match
			p1 = L{pattern_ind};
			p2 = L{patt_i};
			%Truncate one if necessary.
			if length(p1) < length(p2)
				p2 = p2(1:length(p1));
			elseif length(p1) > length(p2)
				p1 = p1(1:length(p2));
			end
			% Add xor
			p_overlap = bitxor(p1,p2);
			%Bit at which things end up being different
			ind_identical = find(p_overlap,1) - 1;
			%Awd constraints
			shared_Q_constrs = shared_Q_constrs +  [Q{pattern_ind}( [1:ind_identical*m] , [1:ind_identical*p] ) == Q{patt_i}( [1:ind_identical*m] , [1:ind_identical*p] )];
			shared_r_constrs = shared_r_constrs +  [r{pattern_ind}( [1:ind_identical*m] ) == r{patt_i}([1:ind_identical*m]) ];
		end

	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim0 = optimize(positive_constr+dual_constrs+l_diag_constr+shared_Q_constrs+shared_r_constrs+obj_constrs, ...
			mu2, ...
			ops);

	opt_out = optim0;
	if opt_out.problem ~= 0
		contr = [];
	else
		% Save Feedback Matrices
		% ++++++++++++++++++++++
		Q_set = {}; r_set = {};
		F_set = {}; u0_set = {};
		for pattern_ind = 1 : length(L)
			T_i = length(L{pattern_ind});
			%Get Parameters
			[H0,S0,Cm0,J0,f_bar,B_w_big,~] = get_mpc_matrices(ad_arr,'word',L{pattern_ind}+1);

			Q_set{pattern_ind} = value(Q{pattern_ind});
			r_set{pattern_ind} = value(r{pattern_ind});
			F_set{pattern_ind} = value( (pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*S0)) ) * Q{pattern_ind});
			u0_set{pattern_ind} = value( pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*S0)) * r{pattern_ind} );
			
			%Fix up F and u0 to avoid NaN
			F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
			% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

			%Create Function Outputs
			opt_out.Q_set = Q_set;
			opt_out.r_set = r_set;

			contr = FHAE_pb(L,F_set,u0_set);

		end
	end

	results.prefix.opt_out = opt_out;
	results.prefix.contr = contr;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize a Time-Based Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear r Q Pi_1 Pi_2 T_i opt_out contr

	mu2 = sdpvar(1,1,'full');

	%%Perform Optimization %%

	% Optimization Variables
	% ++++++++++++++++++++++
	ad0.x0 	= sdpvar(n,1,'full');

	T = length(L{1});
	L_star = ones(1,T);
    for sig_i = 1:length(L)
		L_star = bitand(L_star,L{sig_i});
	end

	Pi_1 = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
	Pi_2 = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');

	% Feedback Variables
	Q = sdpvar(m*T,p*T,'full');
	r = sdpvar(m*T,1,'full');

	obj_constrs = [];
	obj_fcn = [];

	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];

	
	% Creating Constraints
	% ++++++++++++++++++++

	[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(ad0,T,'missing',find(L_star == 0)-1);
	P_wT = 1; P_vT = 1; prod_M2 = 1;
	for t_idx = 1:T
		P_wT = P_wT*ad0.P_w;
		P_vT = P_vT*ad0.P_v;
		prod_M2 = prod_M2*M2;
	end

	positive_constr = positive_constr + [ Pi_1 >= 0, Pi_2 >= 0 ];

	sel_influenced_states = [];
	for i = 1 : T
		%Select all influenced states
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	noise_constrs = noise_constrs + [ Pi_1 * [ P_wT.b ; P_vT.b ; M1.b ] <= mu2*ones(size(prod_M2.A,1),1) - prod_M2.A*sel_influenced_states*H0*r ];
	noise_constrs = noise_constrs + [ Pi_2 * [ P_wT.b ; P_vT.b ; M1.b ] <= M3.b - M3.A*select_m(T,T)*H0*r ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; ad0.A^i];
	end

	G = [ 	(eye(n*(T+1))+H0*Q*Cm0)*S0*B_w_big ...
						H0*Q*C_v_big ...
						(eye(n*(T+1))+H0*Q*Cm0)*pre_xi ];

	bounded_disturb_matrix = [ 	P_wT.A zeros(size(P_wT.A,1),vd*T+n) ;
								zeros(size(P_vT.A,1),size(P_wT.A,2)) P_vT.A zeros(size(P_vT.A,1),n) ;
								zeros(size(M1.A,1),(vd+wd)*T) M1.A ];

	dual_equal_constrs = dual_equal_constrs + [Pi_1 * bounded_disturb_matrix == kron(eye(T),unit_box.A)*sel_influenced_states*G ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == M3.A*select_m(T,T)*G];

	%Lower Diagonal Constraint
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*m+1:bl_row_num*m], [bl_row_num*p+1:end] ) == 0 ];
	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim1 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+obj_constrs, ...
			mu2, ...
			ops);

	opt_out = optim1;
	if opt_out.problem ~= 0
		contr = [];
	else
		% Save Feedback Matrices
		% ++++++++++++++++++++++
		
		[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(ad0,T,'missing',find(L{pattern_ind} == 0)-1);

		opt_out.Q = value(Q);
		opt_out.Q(isnan(value(Q))) = 0;
		opt_out.r = value(r);
		contr.F   = value( (pinv(value(eye(size(Q,1)) + opt_out.Q*Cm0*H0)) ) * opt_out.Q);
		contr.u0  = value( pinv(value(eye(size(Q,1)) + opt_out.Q*Cm0*H0)) * r );
		
		%Fix up F and u0 to avoid NaN
		F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
		% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

 		% contr = FHAE_pb(L,F_set,u0_set);

	end

	results.time.opt_out = opt_out;
	results.time.contr = contr;

	%%%%%%%%%%%%%%%%%%%%%
	%% Function Method %%
	%%%%%%%%%%%%%%%%%%%%%
	%[tempA,tempB] = ad0.rec_synthesis('free' , 'prefix' , 'Min_M2' , M1, M3 , L')

end