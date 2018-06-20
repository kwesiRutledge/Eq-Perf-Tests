function [results] = observer_comparison30(varargin)
%observer_comparison30.m	
%	Description:
%		The objective of this test is to see if ACC has this property that we can use to synthesize "modular" filters in the following sense.
%		A filter will be modular if there exists a time t_1 such that if we use the first portion of the filter from [t_0,t_1], then
%		the second part of the filter can be used from [t_1,T].

	%%%%%%%%%%%%%%%%%%%
	%% Manage Inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin == 0
		verbosity	= 0;
	elseif nargin == 1
		verbosity	= varargin{1};
	else
		error('Unacceptable number of arguments.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	L = {};
	L{1} = [1,1,0,1];

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	%Create Aff_Dyn object with the data from acc_e
	% acc_e = acc;
	% acc_e.B = eye(size(acc.A,1));

	% acc_ad = Aff_Dyn(acc_e.A,acc_e.B,zeros(size(acc_e.A,1),1), acc_e.C, acc_e.d , acc_e.m, acc_e.E , eye(size(acc.C,1)) );

	sys0 = Aff_Dyn(1,1,0,1,0.25,0,-1,1)

	n = size(sys0.A,1);
	m = size(sys0.B,2);
	p = size(sys0.C,1);
	wd = size(sys0.B_w,2);
	vd = size(sys0.C_v,2);

	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	M1 = 1;
	% T = length(L,2);
	ops = sdpsettings('verbose',verbosity);

	%Save Initial Constants to Results
	results.con.sys = sys0;
	results.con.M1 = M1;
	results.con.L = L;

	disp('Constants created and saved to results.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Implement Prefix-Based Filter Design %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('1. There will be 2 problems:')
	disp('		- L = { [0,0] }, an infeasible problem.')
	disp('		- L = { [1,0,0] }, a feasible ER problem.')
	disp(' ')

	[opt_out,contr] = eq_rec_design_pb( sys0 , 'Feasible Set' , M1 , 3 , {[0,0]} );

	results.infeasible_prob.opt_out = opt_out;
	results.infeasible_prob.contr = contr;

	[opt_out,contr] = eq_rec_design_pb( sys0 , 'Feasible Set' , M1 , 3 , {[1,0,0]} );

	results.feasible_prob.opt_out = opt_out;
	results.feasible_prob.contr = contr;

	[opt_out,contr] = free_rec_design_pb( sys0 , 'Min_M3' , M1 , 3 , {[1]} );
	
	results.half1.opt_out = opt_out;
	results.half1.contr = contr;

	[opt_out,contr] = free_rec_design_pb( sys0 , 'Feasible Set' , results.half1.opt_out.M3 , 3 , M1 , {[0,0]} );
	
	results.half2.opt_out = opt_out;
	results.half2.contr = contr;

	disp('============================================================================')
	disp('2. Simulate how one could easily generate a ''modular'' filter in this case.')
	disp(' ')

	num_runs = 1000;

	filter_bad = FHAE_pb( {[0,0]} , { results.feasible_prob.contr.F_set{1}(2:end,2:end) } , { results.feasible_prob.contr.u0_set{1}(2:end) } );
	[run1_raw,run1_norms] = filter_bad.simulate_n_runs( sys0 , results.half1.opt_out.M3 , num_runs );

	run1_norms_mod = [];
	for ind = 1:length(run1_norms)
		run1_norms_mod(:,:,ind) = run1_norms{ind};
	end

	T_temp = length(filter_bad.L{1});

	figure;
	plot([0:T_temp],reshape(run1_norms_mod,T_temp+1,num_runs));
	axis([-0.5 3.5 -0.5 3.5])

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesis of Filter with Sufficient Condition Added %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ad = sys0;
	M2 = 3;
	M3 = M1;

	% Optimization Variables
	% ++++++++++++++++++++++
	ad.x0 	= sdpvar(n,1,'full');

	max_T_i = -1;
	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
		r{pattern_ind} = sdpvar(m*T_i,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(2*n*T_i,2*(wd+vd)*T_i+2*n,'full');
		Pi_2{pattern_ind} = sdpvar(2*n,2*(wd+vd)*T_i+2*n,'full');
		Pi_3{pattern_ind} = sdpvar(2*n*T_i,2*(p+n),'full');

		%Find the maximum T_i
		if T_i > max_T_i
			max_T_i = T_i;
		end
	end
	w	= sdpvar(wd*max_T_i,1,'full');

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];
	obj_constrs = [];

	obj_fcn = [];

	% for pattern_ind = 1 : length(L)
	pattern_ind = 1;
		T_i = length(L{pattern_ind});
		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 , Pi_3{pattern_ind} >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T_i
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
		end

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ ad.eta_w * ones(2*wd*T_i,1) ; ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M2 * ones(2*n*T_i,1) - [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*H0*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ ad.eta_w * ones(2*wd*T_i,1) ; ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M3 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T_i,T_i)*H0*r{pattern_ind} ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T_i
			pre_xi = [ pre_xi ; ad.A^i];
		end

		G{pattern_ind} = [ 	(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*S0*B_w_big ...
							H0*Q{pattern_ind}*C_v_big ...
							(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T_i) ; -eye(wd*T_i) ] zeros(2*wd*T_i,vd*T_i+n) ;
									zeros(2*vd*T_i,wd*T_i) [ eye(vd*T_i) ; -eye(vd*T_i) ] zeros(2*vd*T_i,n) ;
									zeros(2*n,(vd+wd)*T_i) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == [eye(n*T_i); -eye(n*T_i)]*sel_influenced_states*G{pattern_ind} ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T_i,T_i)*G{pattern_ind}];

		%Extra Filter Constraints

		size( Pi_3{pattern_ind}*[ eye(vd) zeros(vd,n) ; -eye(vd) zeros(vd,n) ; zeros(n,vd) eye(n) ; zeros(vd,n) -eye(n)  ] )
		size( sel_influenced_states*[ zeros(n,p) eye(n) ; H0([1:n*(T_i+1-1)],[1:m*(T_i-1)]) * Q{pattern_ind}([m+1:end],[p+1:end])*Cm0([p+1:end],[n+1:end])*pre_xi(1:T_i*n,:)*Q{pattern_ind}([1:m],[1:p]) - H0([1:n*(T_i+1-1)],[1:m*(T_i-1)])*Q{pattern_ind}([m+1:end],[1:p]) ...
																								H0([1:n*(T_i+1-1)],[1:m*(T_i-1)]) * Q{pattern_ind}([m+1:end],[p+1:end])*Cm0([p+1:end],[n+1:end])*pre_xi(1:T_i*n,:)*Q{pattern_ind}([1:m],[1:p])*L{1}(1)*ad.C - H0([1:n*(T_i+1-1)],[1:m*(T_i-1)])*Q{pattern_ind}([m+1:end],[1:p])*L{1}(1)*ad.C	] )

		noise_constrs = noise_constrs + [ Pi_3{pattern_ind} * [ ad.eta_v * ones(2*p,1) ; M1 * ones(2*n,1) ] <= [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*[ zeros(n,m*(T_i)) ; eye(m*(T_i)) + H0([1:n*(T_i+1-1)],[1:m*(T_i-1)]) * Q{pattern_ind}([m+1:end],[p+1:end])*Cm0([p+1:end],[n+1:end]) ]*pre_xi(1:T_i*n,:)*r{pattern_ind}(1:m) ];
		dual_equal_constrs = dual_equal_constrs + [Pi_3{pattern_ind}*[ eye(vd) zeros(vd,n) ; -eye(vd) zeros(vd,n) ; zeros(n,vd) eye(n) ; zeros(vd,n) -eye(n)  ] == [eye(n*T_i);-eye(n*T_i)]* ...
													sel_influenced_states*[ zeros(n,p) eye(n) ; H0([1:n*(T_i+1-1)],[1:m*(T_i-1)]) * Q{pattern_ind}([m+1:end],[p+1:end])*Cm0([p+1:end],[n+1:end])*pre_xi(1:T_i*n,:)*Q{pattern_ind}([1:m],[1:p]) - H0([1:n*(T_i+1-1)],[1:m*(T_i-1)])*Q{pattern_ind}([m+1:end],[1:p]) ...
																								H0([1:n*(T_i+1-1)],[1:m*(T_i-1)]) * Q{pattern_ind}([m+1:end],[p+1:end])*Cm0([p+1:end],[n+1:end])*pre_xi(1:T_i*n,:)*Q{pattern_ind}([1:m],[1:p])*L{1}(1)*ad.C - H0([1:n*(T_i+1-1)],[1:m*(T_i-1)])*Q{pattern_ind}([m+1:end],[1:p])*L{1}(1)*ad.C	] ...
																								 ]

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

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim0 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+shared_Q_constrs+shared_r_constrs+obj_constrs, ...
			obj_fcn, ...
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
			[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

			Q_set{pattern_ind} = value(Q{pattern_ind});
			r_set{pattern_ind} = value(r{pattern_ind});
			F_set{pattern_ind} = value( (inv(value(eye(size(S0,2)) + Q{pattern_ind}*Cm0*H0)) ) * Q{pattern_ind});
			u0_set{pattern_ind} = value( inv(value(eye(size(S0,2)) + Q{pattern_ind}*Cm0*H0)) * r{pattern_ind} );
			
			%Fix up F and u0 to avoid NaN
			F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
			% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

			%Create Function Outputs
			opt_out.Q_set = Q_set;
			opt_out.r_set = r_set;

			contr = FHAE_pb(L,F_set,u0_set);

		end
	end

	results.extra_cond.opt_out = opt_out;
	results.extra_cond.contr = contr;

end