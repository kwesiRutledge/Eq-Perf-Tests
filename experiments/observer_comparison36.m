function [results] = observer_comparison36(varargin)
	%observer_comparison36.m
	%	Verifying the behavior of a simple Lane Keeping system function generator.

	%% Constnats
	dt = 0.01;

	[~,lk_dyn] = get_lk_aff_dyn(dt);

	M1 = 0.1;

	extra_time = 5;
	% L = {ones(1,extra_time)}
	L = {[1,0,ones(1,extra_time)],[1,1,0,ones(1,extra_time-1)]};

	num_runs = 100;

	%% Sanity Checking of the Model

	%% Synthesis

	[ oc36_opt1 , oc36_contr1 ] = eq_rec_design_pb( lk_dyn , 'Min_M2' , M1 , L );

	%% Simulate
	disp('Ideally, I would''ve liked to put this into the NAHS publication. It seems as if the controller can not guarantee very low values for M2.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Double Integrator Formulation %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear M1 n_x m p dt L

	L = {	[1,0,1,0,1,0,1],...
			[1,1,0,1,0,1,0], ...
			[1,0,1,1,1,0,1], ...
			[1,0,0,0,1,1,1]};

	n_x = 4;
	A = [ zeros(n_x,1) [1;-20;0;0] zeros(n_x,1) [0;0;1;-20] ];
	B = [ 	0, 0;
			1, 0;
			0, 0;
			0, 1];
	F = zeros(4,1);
	C = [1,0,0,0;0,0,1,0];
	n_y = size(C,1);

	temp_sys = ss(A,B,C,0);
	dt = 0.05;
	temp_dsys = c2d(temp_sys,dt);

	eta_w = 0.1; eta_v = 0.2;
	M1 = 0.3;

	simp_int_dyn = Aff_Dyn(	temp_dsys.A,temp_dsys.B,F,C,...
							eta_w,eta_v, ...
							temp_dsys.B, eye(n_y))

	[ oc36_opt2 , oc36_contr2 ] = simp_int_dyn.eq_rec_design_pb( 'Min_M2' , M1 , L );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Double Integrator "Steered Invariance" %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	clear L
	extra_time = 10;
	L = {	[1,0,ones(1,extra_time)],...
			[1,1,0,ones(1,extra_time-1)],...
			[1,1,1,0,ones(1,extra_time-2)],...
			[1,1,1,1,0,ones(1,extra_time-3)]};

	%target sets
	M2 = 2;
	% M1 = 0.3
	poly1 = Polyhedron('lb',-M1*ones(1,n_x),'ub',M1*ones(1,n_x)) + [1;0;1;0];
	poly2 = Polyhedron('lb',-M2*ones(1,n_x),'ub',M2*ones(1,n_x));

	%constants
	n = n_x;
	m = size(simp_int_dyn.B,2);
	p = size(simp_int_dyn.C,1);
	wd = size(simp_int_dyn.B_w,2);
	vd = size(simp_int_dyn.C_v,2);

	max_T_i = -1;

	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	verbosity = 1;
	ops = sdpsettings('verbose',verbosity);

	%Begin Synthesis

	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
		r{pattern_ind} = sdpvar(m*T_i,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(size(poly2.b,1)*T_i,2*(wd+vd)*T_i+2*n,'full');
		Pi_2{pattern_ind} = sdpvar(size(poly1.b,1),2*(wd+vd)*T_i+2*n,'full');

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

	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(simp_int_dyn,T_i,'missing',find(L{pattern_ind} == 0)-1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T_i
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
		end

		%Create special diagnoalization of poly1.A
		temp_larg_diag = [];
		for i = 1:T_i
			temp_larg_diag = [ temp_larg_diag ; zeros(size(poly2.A,1),(i-1)*n) poly2.A zeros(size(poly2.A,1),(T_i-i)*n) ];
		end

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ simp_int_dyn.eta_w * ones(2*wd*T_i,1) ; simp_int_dyn.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= repmat(poly2.b,T_i,1) - temp_larg_diag*sel_influenced_states*H0*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ simp_int_dyn.eta_w * ones(2*wd*T_i,1) ; simp_int_dyn.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= poly1.b - poly1.A*select_m(T_i,T_i)*H0*r{pattern_ind} ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T_i
			pre_xi = [ pre_xi ; simp_int_dyn.A^i];
		end

		G{pattern_ind} = [ 	(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*S0*B_w_big ...
							H0*Q{pattern_ind}*C_v_big ...
							(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T_i) ; -eye(wd*T_i) ] zeros(2*wd*T_i,vd*T_i+n) ;
									zeros(2*vd*T_i,wd*T_i) [ eye(vd*T_i) ; -eye(vd*T_i) ] zeros(2*vd*T_i,n) ;
									zeros(2*n,(vd+wd)*T_i) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == temp_larg_diag*sel_influenced_states*G{pattern_ind} ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == poly1.A*select_m(T_i,T_i)*G{pattern_ind}];

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
	optim0 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+shared_Q_constrs+shared_r_constrs, ...
			[], ...
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
			[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(simp_int_dyn,T_i,'missing',find(L{pattern_ind} == 0)-1);

			Q_set{pattern_ind} = value(Q{pattern_ind});
			r_set{pattern_ind} = value(r{pattern_ind});
			F_set{pattern_ind} = value( (inv(value(eye(size(H0,2)) + Q{pattern_ind}*Cm0*H0)) ) * Q{pattern_ind});
			u0_set{pattern_ind} = value( inv(value(eye(size(H0,2)) + Q{pattern_ind}*Cm0*H0)) * r{pattern_ind} );
			
			%Fix up F and u0 to avoid NaN
			F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
			% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

		end

		%Create Function Outputs
		opt_out.Q_set = Q_set;
		opt_out.r_set = r_set;

		contr = FHAE_pb(L,F_set,u0_set);
	end

	%Plot Target Set and intermediate set
	start_poly = poly1 - [1;0;1;0];
	figure;
	hold on;
	plot(poly2.slice([2,4],[0,0]))
	plot(poly1.slice([2,4],[0,0]),'color','g')

	%Simulate
	[ ctrl_sim1, ~ ] = contr.simulate_n_runs( simp_int_dyn , M1 , num_runs )

	ctrl_sim1_mod = [];
	for ind = 1:length(ctrl_sim1)
		ctrl_sim1_mod(:,:,ind) = simp_int_dyn.C*ctrl_sim1{ind};
	end

	figure;
	hold on;
	plot(poly2.slice([2,4],[0,0]),'color','w')
	plot(poly1.slice([2,4],[0,0]),'color','g')
	plot(start_poly.slice([2,4],[0,0]),'color','y')

	for sim_num = 1:num_runs
		plot(ctrl_sim1_mod(1,:,sim_num),ctrl_sim1_mod(2,:,sim_num))
	end

	%Results
	results.exp1.sys = lk_dyn;
	results.exp1.opt = oc36_opt1;
	results.exp1.contr = oc36_contr1;

	results.exp2.sys = simp_int_dyn;
	results.exp2.opt = oc36_opt2;
	results.exp2.contr = oc36_contr2;

	results.exp3.sys = simp_int_dyn;
	results.exp3.p1 = poly1;
	results.exp3.p2 = poly2;
	results.exp3.opt_out = optim0;
end