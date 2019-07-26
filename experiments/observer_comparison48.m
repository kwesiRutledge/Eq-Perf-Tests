function [results] = observer_comparison48(varargin)
	%Description:
	%	Meant to test the new constraint generator which should allow for the equalized recovery function
	%	to accept Polyhedra and eventually zonotopes.
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Trying to show the Lane-Keeping Example %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear r

	dt = 0.1;
	[~,lk_ad] = get_lk_aff_dyn(dt);

	L = {ones(1,6)};
	mu1 = 1;
	mu2 = 10;
	mu3 = 3;

	M1 = Polyhedron('lb',-1*[0.5,0.5,0.05,0.05],'ub',[0.5,0.5,0.05,0.05]);
	M2 = 3*M1;
	M3 = M1;

	% Constants
	n = size(lk_ad.A,1);
	m = size(lk_ad.B,2);
	p = size(lk_ad.C,1);
	wd = size(lk_ad.B_w,2);
	vd = size(lk_ad.C_v,2);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%%Perform Optimization %%

	% Optimization Variables
	% ++++++++++++++++++++++
	lk_ad.x0 	= sdpvar(n,1,'full');

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
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];

	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(lk_ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		sel_influenced_states = [];
		prod_M2 = 1;
		for i = 1 : T_i
			%Select all influenced states
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
			%Create product for M2
			prod_M2 = prod_M2*M2;
		end

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ lk_ad.eta_w * ones(2*wd*T_i,1) ; lk_ad.eta_v * ones(2*p*T_i,1) ; M1.b ] <= prod_M2.b - prod_M2.A*sel_influenced_states*H0*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ lk_ad.eta_w * ones(2*wd*T_i,1) ; lk_ad.eta_v * ones(2*p*T_i,1) ; M1.b ] <= M3.b - M3.A*select_m(T_i,T_i)*H0*r{pattern_ind} ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T_i
			pre_xi = [ pre_xi ; lk_ad.A^i];
		end

		G{pattern_ind} = [ 	(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*S0*B_w_big ...
							H0*Q{pattern_ind}*C_v_big ...
							(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T_i) ; -eye(wd*T_i) ] zeros(2*wd*T_i,vd*T_i+n) ;
									zeros(2*vd*T_i,wd*T_i) [ eye(vd*T_i) ; -eye(vd*T_i) ] zeros(2*vd*T_i,n) ;
									zeros(size(M1.A,1),(vd+wd)*T_i) M1.A ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == prod_M2.A*sel_influenced_states*G{pattern_ind} ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == M3.A*select_m(T_i,T_i)*G{pattern_ind}];

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
			[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(lk_ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

			Q_set{pattern_ind} = value(Q{pattern_ind});
			r_set{pattern_ind} = value(r{pattern_ind});
			F_set{pattern_ind} = value( (pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*H0)) ) * Q{pattern_ind});
			u0_set{pattern_ind} = value( pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*H0)) * r{pattern_ind} );
			
			%Fix up F and u0 to avoid NaN
			F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
			% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

			%Create Function Outputs
			opt_out.Q_set = Q_set;
			opt_out.r_set = r_set;

			contr = FHAE_pb(L,F_set,u0_set);

		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Try to Reproduce with Function %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[opt_out2,contr2] = lk_ad.rec_synthesis('Equalized','prefix','Feasible Set', M1 , M2 , L );

	results.fcn_result.opt_out = opt_out2;
	results.fcn_result.contr = contr2;

end